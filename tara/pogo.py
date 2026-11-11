"""Delay-aware scalar longitudinal feedback screening for a supplied loop model."""
from __future__ import annotations

from math import factorial

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import root


def _polynomial(values: ArrayLike) -> np.ndarray:
    result = np.trim_zeros(np.asarray(values, dtype=float), "f")
    if result.ndim != 1 or not len(result) or not np.isfinite(result).all():
        raise ValueError("Transfer polynomials must be nonzero finite coefficient vectors")
    return result


def assess(
    structural_num: ArrayLike, structural_den: ArrayLike,
    feed_num: ArrayLike, feed_den: ArrayLike, delay_s: float,
) -> dict:
    """Evaluate 1 + Gs Zf exp(-delay*s) for a strictly proper rational Gs Zf.

    Coefficients descend in powers of s in 1/s. The negative feedback sign and
    dimensionless loop gain must be encoded consistently in the supplied model.
    Nyquist counts right-half-plane roots. Reported eigenvalues are a subset of
    the delay spectrum, checked against the exact characteristic equation.
    """
    num = np.polymul(_polynomial(structural_num), _polynomial(feed_num))
    den = np.polymul(_polynomial(structural_den), _polynomial(feed_den))
    if not np.isfinite(delay_s) or delay_s < 0 or len(num) >= len(den):
        raise ValueError("Pogo analysis needs a nonnegative delay and a strictly proper loop")
    poles, zeros = np.roots(den), np.roots(num)
    scale = max(1., float(np.max(abs(poles), initial=0)), float(np.max(abs(zeros), initial=0)))
    if np.any(abs(poles.real) < 1e-9*scale):
        raise ValueError("Pogo analysis requires open-loop poles off the imaginary axis")
    if any(np.min(abs(z - poles), initial=np.inf) < 1e-8*scale for z in zeros):
        raise ValueError("Supply a minimal loop with explicit pole-zero cancellations removed")

    def characteristic(s: np.ndarray | complex) -> np.ndarray | complex:
        return 1 + np.polyval(num, s)/np.polyval(den, s)*np.exp(-delay_s*s)

    radius = 2*scale
    for _ in range(30):
        bound = (abs(num[0]/den[0])*radius**(len(num)-len(den))
                 *np.prod(1 + abs(zeros)/radius)/np.prod(1 - abs(poles)/radius))
        if bound < 0.25:
            break
        radius *= 2
    else:
        raise RuntimeError("Could not bound the high-frequency Nyquist closure")

    # Close in the right half-plane. Since |L|<1/4 on the arc, the
    # corresponding 1+L path can be contracted to the direct chord.
    previous = None
    rhp_poles = int(np.count_nonzero(poles.real > 0))
    for n in (2049, 4097, 8193, 16385, 32769, 65537, 131073):
        if 2*radius*delay_s/(n-1) > 0.15:
            continue
        values = characteristic(1j*np.linspace(-radius, radius, n))
        if np.min(abs(values)) < 1e-6:
            break
        angle = np.angle(np.roll(values, -1)/values)
        winding = int(np.rint(angle.sum()/(2*np.pi)))
        count = rhp_poles - winding
        # Exclude the contractible chord from the local sampling diagnostic.
        if count == previous and count >= 0 and np.max(abs(angle[:-1])) < 0.25:
            break
        previous = count
    else:
        n = 0
    resolved = (n > 0 and previous is not None and np.min(abs(values)) >= 1e-6
                and count == previous and count >= 0 and np.max(abs(angle[:-1])) < 0.25)

    # Seed an exact-root check using a scaled low-order Padé representation.
    n_order = 6 if delay_s else 0
    coeff = np.array([
        factorial(2*n_order-k)*factorial(n_order)
        / (factorial(2*n_order)*factorial(n_order-k)*factorial(k))*(delay_s*scale)**k
        for k in range(n_order+1)
    ])
    scaled_den = den*scale**np.arange(len(den)-1, -1, -1)
    scaled_num = num*scale**np.arange(len(num)-1, -1, -1)
    approximate = np.polyadd(np.polymul(scaled_den, coeff[::-1]),
                             np.polymul(scaled_num, (coeff*(-1.)**np.arange(len(coeff)))[::-1]))
    roots = []
    for seed in sorted(np.roots(approximate)*scale, key=lambda x: -x.real)[:18]:
        if not delay_s:
            candidate = seed
        else:
            def residual(v: np.ndarray) -> list[float]:
                try:
                    with np.errstate(over="raise", divide="raise", invalid="raise"):
                        value = characteristic(complex(*v))
                except FloatingPointError:
                    return [1e100, 1e100]
                return [float(value.real), float(value.imag)]

            fit = root(residual, [seed.real, seed.imag])
            if np.linalg.norm(residual(fit.x)) > 1e-7:
                continue
            candidate = complex(*fit.x)
        error = float(abs(characteristic(candidate)))
        if error <= 1e-7 and all(abs(candidate - z[0]) > 1e-7*scale for z in roots):
            roots.append((candidate, error))

    return dict(
        status=("UNRESOLVED" if not resolved else ("RHP_INSTABILITY" if count else "NO_RHP_ROOTS_DETECTED")),
        open_loop_RHP_poles=rhp_poles, closed_loop_RHP_roots=count if resolved else None,
        nyquist_radius_radps=radius, nyquist_samples=n, high_frequency_bound=bound,
        minimum_sampled_distance_to_minus_one=float(np.min(abs(values))) if n else None,
        eigenvalue_scope="subset_checked_against_exact_delay_characteristic",
        eigenvalues=[dict(real_per_s=float(z.real), imag_radps=float(z.imag), characteristic_residual=e)
                     for z, e in sorted(roots, key=lambda pair: (-pair[0].real, pair[0].imag))],
    )
