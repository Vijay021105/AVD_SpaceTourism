"""Free-free, planar Timoshenko beam model for supplied homogenized sections."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh

Array = NDArray[np.float64]


@dataclass(frozen=True)
class BeamSection:
    """One uniform section, with all properties in SI.

    EI and kGA describe one bending plane. Rotary inertias are per unit length.
    Wet properties represent rigidly attached liquid. A separate slosh model
    is required for moving liquid relative to the tank.
    """

    length_m: float
    EA_N: float
    EI_Nm2: float
    kGA_N: float
    dry_mass_per_m: float
    wet_mass_per_m: float
    dry_rotary_inertia_per_m: float
    wet_rotary_inertia_per_m: float
    elements: int


def modes(sections: list[BeamSection], wet: bool, count: int = 8) -> dict:
    """Assemble free-free mass/stiffness and return elastic modes and residuals.

    DOFs per node are axial displacement, lateral displacement and rotation.
    Reduced integration of the linear element's shear energy avoids locking.
    Structural interfaces are rigid and mesh convergence must be checked.
    """
    if not sections or count < 1:
        raise ValueError("A positive mode count and at least one beam section are required")
    for s in sections:
        properties = [s.length_m, s.EA_N, s.EI_Nm2, s.kGA_N,
                      s.dry_mass_per_m, s.wet_mass_per_m,
                      s.dry_rotary_inertia_per_m, s.wet_rotary_inertia_per_m]
        if (not np.isfinite(properties).all() or min(properties) <= 0
                or not isinstance(s.elements, int) or s.elements < 2
                or s.wet_mass_per_m < s.dry_mass_per_m
                or s.wet_rotary_inertia_per_m < s.dry_rotary_inertia_per_m):
            raise ValueError("Invalid section properties, wet mass or element count")
    total = sum(s.elements for s in sections)
    if total > 300:
        raise ValueError("Dense beam model supports at most 300 elements")
    z = np.zeros(total+1)
    K = np.zeros((3*(total+1),)*2)
    M = np.zeros_like(K)
    i = 0
    for s in sections:
        L = s.length_m/s.elements
        mass = s.wet_mass_per_m if wet else s.dry_mass_per_m
        inertia = s.wet_rotary_inertia_per_m if wet else s.dry_rotary_inertia_per_m
        extension = np.array([-1/L, 0, 0, 1/L, 0, 0])
        curvature = np.array([0, 0, -1/L, 0, 0, 1/L])
        shear = np.array([0, -1/L, -0.5, 0, 1/L, -0.5])
        ke = L*(s.EA_N*np.outer(extension, extension)
                + s.EI_Nm2*np.outer(curvature, curvature)
                + s.kGA_N*np.outer(shear, shear))
        me = L/6*np.kron([[2., 1.], [1., 2.]], np.diag([mass, mass, inertia]))
        for _ in range(s.elements):
            dofs = slice(3*i, 3*i+6)
            K[dofs, dofs] += ke
            M[dofs, dofs] += me
            z[i+1] = z[i] + L
            i += 1

    # Explicit rigid basis avoids classifying the first low elastic root as a
    # rigid root in a beam with a very broad stiffness or mass spectrum.
    from scipy.linalg import cholesky, null_space, solve_triangular
    rigid = np.zeros((len(K), 3))
    rigid[0::3, 0] = rigid[1::3, 1] = rigid[2::3, 2] = 1
    rigid[1::3, 2] = z
    chol = cholesky(M, lower=True)
    basis = null_space((chol.T @ rigid).T)
    transform = solve_triangular(chol.T, basis, lower=False)
    vals, shapes = eigh(transform.T @ K @ transform,
                        subset_by_index=(0, min(count, len(K)-3)-1))
    if np.min(vals) <= 0:
        raise RuntimeError("Beam elastic eigenspectrum contains a nonpositive root")
    phi = transform @ shapes
    left = K @ phi
    right = (M @ phi)*vals
    error = np.linalg.norm(left-right, axis=0)/np.maximum(
        np.linalg.norm(left, axis=0) + np.linalg.norm(right, axis=0), 1e-300)
    if not np.isfinite(error).all() or np.max(error) > 1e-6:
        raise RuntimeError("Beam eigenpair residual exceeds tolerance")
    axial_fraction = np.sum(phi[0::3]*(M @ phi)[0::3], axis=0)
    return dict(z_m=z, M=M, K=K, modes=phi.reshape(-1, 3, len(vals)),
                frequency_Hz=np.sqrt(vals)/(2*np.pi), axial_kinetic_fraction=axial_fraction,
                eigen_residual=error)
