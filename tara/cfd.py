"""Steady axisymmetric Euler finite volumes for a calorically perfect gas.

The domain ends at the base plane. Wall forces include only pressure on the
resolved meridian, at zero incidence. No viscous shear, base wake, or heat
transfer is present. Local pseudo-time stepping is a steady solver accelerator.
"""
from __future__ import annotations

from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from .geometry import RingMesh

Array = NDArray[np.float64]
Flux = Literal["rusanov", "roe"]
Reconstruction = Literal["first_order", "muscl_minmod"]


def primitive(U: Array, gamma: float = 1.4) -> tuple[Array, Array, Array, Array]:
    """Return density, axial/radial velocities and thermodynamic pressure in SI."""
    rho = U[..., 0]
    u, v = U[..., 1] / rho, U[..., 2] / rho
    return rho, u, v, (gamma - 1) * (U[..., 3] - 0.5 * rho * (u*u + v*v))


def state(rho: Any, u: Any, v: Any, p: Any, gamma: float = 1.4) -> Array:
    """Build conservative state vectors ordered rho, rho*u, rho*v, rho*E."""
    rho, u, v, p = np.broadcast_arrays(rho, u, v, p)
    return np.stack((rho, rho*u, rho*v, p/(gamma-1) + 0.5*rho*(u*u+v*v)), axis=-1)


def flux(
    UL: Array, UR: Array, S: Array, gamma: float, scheme: Flux = "rusanov"
) -> tuple[Array, Array]:
    """Return an oriented, area-integrated inviscid flux and CFL spectral bound."""
    area = np.linalg.norm(S, axis=-1)
    nx, nr = np.moveaxis(S / np.maximum(area[..., None], 1e-100), -1, 0)
    rl, ul, vl, pl = primitive(UL, gamma)
    rr, ur, vr, pr = primitive(UR, gamma)
    ql, qr = ul*nx + vl*nr, ur*nx + vr*nr
    hl, hr = (UL[..., 3] + pl) / rl, (UR[..., 3] + pr) / rr
    FL = np.stack((rl*ql, rl*ul*ql + pl*nx, rl*vl*ql + pl*nr, rl*hl*ql), axis=-1)
    FR = np.stack((rr*qr, rr*ur*qr + pr*nx, rr*vr*qr + pr*nr, rr*hr*qr), axis=-1)
    speed = np.maximum(abs(ql) + np.sqrt(gamma*pl/rl), abs(qr) + np.sqrt(gamma*pr/rr))
    dissipation = speed[..., None] * (UR - UL)
    if scheme == "roe":
        left, right = np.sqrt(rl), np.sqrt(rr)
        weight = left + right
        u = (left*ul + right*ur) / weight
        v = (left*vl + right*vr) / weight
        h = (left*hl + right*hr) / weight
        kinetic = 0.5 * (u*u + v*v)
        a = np.sqrt((gamma - 1) * (h - kinetic))
        q, tangent = u*nx + v*nr, -u*nr + v*nx
        dp, dq, density = pr - pl, qr - ql, left * right
        minus = (dp - density*a*dq) / (2*a*a)
        plus = (dp + density*a*dq) / (2*a*a)
        contact = rr - rl - dp/(a*a)
        shear = density * (-(ur-ul)*nr + (vr-vl)*nx)

        def entropy_fix(speed: Array) -> Array:
            delta = 0.1 * a
            return np.where(abs(speed) < delta, 0.5*(speed*speed/delta + delta), abs(speed))

        acoustic_minus = np.stack((np.ones_like(q), u-a*nx, v-a*nr, h-a*q), axis=-1)
        acoustic_plus = np.stack((np.ones_like(q), u+a*nx, v+a*nr, h+a*q), axis=-1)
        entropy_wave = np.stack((np.ones_like(q), u, v, kinetic), axis=-1)
        shear_wave = np.stack((np.zeros_like(q), -nr, nx, tangent), axis=-1)
        roe = ((minus*entropy_fix(q-a))[..., None]*acoustic_minus
               + (plus*entropy_fix(q+a))[..., None]*acoustic_plus
               + abs(q)[..., None]*(contact[..., None]*entropy_wave + shear[..., None]*shear_wave))
        # Switch to the more dissipative flux if Roe's first or last acoustic
        # intermediate state has nonpositive density or internal energy.
        valid = _physical(UL + minus[..., None]*acoustic_minus)
        valid &= _physical(UR - plus[..., None]*acoustic_plus)
        dissipation = np.where(valid[..., None], roe, dissipation)
    elif scheme != "rusanov":
        raise ValueError(f"Unknown flux scheme {scheme!r}")
    return 0.5*(FL + FR - dissipation)*area[..., None], speed*area


def _physical(U: Array) -> NDArray[np.bool_]:
    """Positivity without dividing by density or clipping conserved states."""
    rho = U[..., 0]
    return (np.isfinite(U).all(axis=-1) & (rho > 0)
            & (2*rho*U[..., 3] > U[..., 1]**2 + U[..., 2]**2))


def _faces(U: Array, axis: int, gamma: float, reconstruction: Reconstruction) -> tuple[Array, Array]:
    """Reconstruct one-sided face states with minmod slopes on the logical mesh."""
    if reconstruction == "first_order":
        return U, U
    if reconstruction != "muscl_minmod":
        raise ValueError(f"Unknown reconstruction {reconstruction!r}")
    W = np.moveaxis(np.stack(primitive(U, gamma), axis=-1), axis, 0)
    delta = np.diff(W, axis=0)
    slope = np.zeros_like(W)
    slope[1:-1] = np.where(
        delta[:-1]*delta[1:] > 0,
        np.sign(delta[:-1])*np.minimum(abs(delta[:-1]), abs(delta[1:])),
        0,
    )
    left = state(*np.moveaxis(W - 0.5*slope, -1, 0), gamma)
    right = state(*np.moveaxis(W + 0.5*slope, -1, 0), gamma)
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


def rhs(
    U: Array, mesh: RingMesh, inf: Array, gamma: float,
    scheme: Flux = "rusanov", reconstruction: Reconstruction = "first_order",
) -> tuple[Array, Array]:
    """Evaluate the semi-discrete Euler operator including the ring metric source."""
    x_left, x_right = _faces(U, 0, gamma, reconstruction)
    r_left, r_right = _faces(U, 1, gamma, reconstruction)
    inflow = np.broadcast_to(inf, (1, mesh.nr, 4))
    pin = float(primitive(inf, gamma)[3])
    outlet = U[-1:].copy()
    rho, u, v, p = primitive(outlet, gamma)
    subsonic = u < np.sqrt(gamma*p/rho)
    outlet[..., 3] = np.where(subsonic, pin/(gamma-1) + 0.5*rho*(u*u+v*v), outlet[..., 3])
    Fx, sx = flux(np.concatenate((inflow, x_right)), np.concatenate((x_left, outlet)), mesh.Sx, gamma, scheme)

    outer = U[:, -1:].copy()
    rho, u, v, p = primitive(outer, gamma)
    un = u*mesh.Sr[:, -1:, 0] + v*mesh.Sr[:, -1:, 1]
    normal_speed = un / np.linalg.norm(mesh.Sr[:, -1:], axis=-1)
    normal_mach = normal_speed / np.sqrt(gamma*p/rho)
    # Rusanov diffuses contact variables even at zero normal speed. A hard
    # inflow/outflow switch at a grazing boundary makes the residual discontinuous.
    inflow_weight = np.clip(0.5 - normal_mach/0.04, 0, 1)[..., None]
    outer = inflow_weight*inf + (1 - inflow_weight)*outer
    rho, u, v, p = primitive(outer, gamma)
    outer[..., 3] = np.where(
        normal_speed < np.sqrt(gamma*p/rho),
        pin/(gamma-1) + 0.5*rho*(u*u+v*v), outer[..., 3],
    )
    Fr, sr = flux(
        np.concatenate((r_left[:, :1], r_right), axis=1),
        np.concatenate((r_left, outer), axis=1), mesh.Sr, gamma, scheme,
    )
    p = primitive(U, gamma)[3]
    Fr[:, 0] = 0
    Fr[:, 0, 1:3] = p[:, 0, None]*mesh.Sr[:, 0]
    rate = -(Fx[1:] - Fx[:-1] + Fr[:, 1:] - Fr[:, :-1])/mesh.V[..., None]
    rate[..., 2] += 2*np.pi*p*mesh.area2/mesh.V
    return rate, sx[1:] + sx[:-1] + sr[:, 1:] + sr[:, :-1]


def admissibility(row: dict, c: dict) -> tuple[list[str], float]:
    """Screen continuum, low-Mach and perfect-gas model limits at a flight state."""
    settings = c["cfd"]
    mach = row["Mach"]
    total_temp = row["T_K"] * (1 + 0.5*(settings["gamma"]-1)*mach*mach)
    reasons = []
    if mach < settings["minimum_mach"]:
        reasons.append("low_Mach_no_preconditioning")
    if mach > settings["maximum_mach"]:
        reasons.append("beyond_Mach_range")
    if total_temp > settings["maximum_total_temperature_K"]:
        reasons.append("calorically_perfect_gas_limit")
    if row["Kn_D"] > settings["maximum_Knudsen"]:
        reasons.append("rarefied_flow_requires_kinetic_model")
    if not row["atmosphere_valid"]:
        reasons.append("atmosphere_extrapolated")
    return reasons, total_temp


def solve(
    c: dict, row: dict, nx: int = 140, nr: int = 32, steps: int = 3000,
    perturb: tuple[Array, Array] | None = None, initial: Array | None = None,
) -> dict:
    """Advance to steady state or the iteration cap and report the actual outcome."""
    settings = c["cfd"]
    g = settings["gamma"]
    scheme = settings.get("flux_scheme", "rusanov")
    reconstruction = settings.get("reconstruction", "first_order")
    local_cfl = settings.get("local_cfl", settings["cfl"])
    required = int(settings.get("consecutive_checkpoints", 3))
    reference = np.array([row[k] for k in ("rho_kgm3", "v_mps", "p_Pa", "a_mps", "q_Pa")])
    if (not np.isfinite(reference).all() or min(reference) <= 0
            or not 0 < local_cfl <= 1 or not 1 < g < 2
            or required < 1 or steps < 1 or nx < 12 or nr < 2):
        raise ValueError("Invalid CFD state, mesh, iteration budget, or CFL")
    if min(settings["residual_tolerance"], settings["drag_change_tolerance"]) <= 0:
        raise ValueError("CFD convergence tolerances must be positive")
    mesh = RingMesh(c, int(row["stage"]), nx, nr, perturb)
    rho, speed, pressure, sound, q = reference
    inf = state(rho, speed, 0, pressure, g)
    U = np.broadcast_to(inf, (nx, nr, 4)).copy() if initial is None else np.array(initial, dtype=float)
    if U.shape != (nx, nr, 4) or not _physical(U).all():
        raise ValueError("Initial CFD state must match mesh with positive density and pressure")
    area = np.pi*c["diameter_m"]**2/4
    velocity_scale = max(speed, sound)
    rate_scale = np.array([rho, rho*velocity_scale, rho*velocity_scale, max(pressure, q)])
    rate_scale *= velocity_scale / (mesh.x[-1] - mesh.x[0])
    history = []
    tau = 0.0
    drag_old = None
    streak = retries = 0
    for it in range(steps + 1):
        rate, wave = rhs(U, mesh, inf, g, scheme, reconstruction)
        if it % 50 == 0 or it == steps:
            wallp = primitive(U[:, 0], g)[3]
            drag = float(-np.sum((wallp - pressure)*mesh.Sr[:, 0, 0]))
            scaled = rate/rate_scale
            residual = float(np.sqrt(np.mean(scaled**2)))
            change = abs(drag-drag_old)/max(abs(drag), 0.01*q*area) if drag_old is not None else None
            mass_error = abs(np.sum(rate[..., 0]*mesh.V))/(rho*speed*area)
            history.append(dict(
                iteration=it, pseudo_time_s=tau, residual=residual,
                residual_Linf=float(np.max(abs(scaled))), pressure_drag_N=drag,
                CD_pressure_only=drag/(q*area), drag_relative_change=change,
                mass_flux_imbalance=mass_error,
            ))
            within_limits = (it >= 500 and residual < settings["residual_tolerance"]
                             and change is not None and change < settings["drag_change_tolerance"])
            streak = streak + 1 if within_limits else 0
            drag_old = drag
            if streak >= required:
                break
        if it == steps:
            break
        dt = local_cfl*mesh.V/np.maximum(wave, 1e-100)
        for _ in range(15):
            trial = U + dt[..., None]*rate
            if _physical(trial).all():
                break
            dt *= 0.5
            retries += 1
        else:
            raise RuntimeError("Euler positivity retries exhausted")
        U = trial
        tau += float(np.min(dt))

    density, u, v, p = primitive(U, g)
    mach = np.hypot(u, v)/np.sqrt(g*p/density)
    surface = [
        dict(x_from_nose_m=float(x), z_from_original_aft_m=float(c["tip_z_m"]-x),
             radius_m=float(r), pressure_Pa=float(pw), Cp=float((pw-pressure)/q))
        for x, r, pw in zip(mesh.xc[:, 0], 0.5*(mesh.R[:-1, 0]+mesh.R[1:, 0]), p[:, 0])
    ]
    final = history[-1]
    summary = dict(
        converged=streak >= required, iterations=it, pseudo_time_s=tau, residual=final["residual"],
        pressure_drag_N=final["pressure_drag_N"], CD_pressure_only=final["CD_pressure_only"],
        drag_relative_change=final["drag_relative_change"], residual_Linf=final["residual_Linf"],
        mass_flux_imbalance=final["mass_flux_imbalance"], positivity_retries=retries,
        ring_metric_error=mesh.check(), Nx=nx, Nr=nr, flux_scheme=scheme,
        reconstruction=reconstruction, residual_tolerance=settings["residual_tolerance"],
        drag_change_tolerance=settings["drag_change_tolerance"], convergence_streak=streak,
        local_cfl=local_cfl, time_coordinate="sum_of_minimum_cell_steps_not_physical_time",
    )
    return dict(mesh=mesh, U=U, p=p, Mach=mach, surface=surface, history=history, summary=summary)
