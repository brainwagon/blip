"""Support radius sweep optimizer — assembles K/f once, re-solves per radius."""

import numpy as np
from plate_fem import (
    create_circular_mesh,
    assemble_plate_system,
    find_support_dofs,
    solve_plate,
)
from rms_calc import compute_rms


# Material properties for Pyrex (borosilicate glass)
PYREX_E = 63e9          # Young's modulus (Pa)
PYREX_NU = 0.20         # Poisson's ratio
PYREX_RHO = 2230.0      # Density (kg/m^3)
GRAVITY = 9.81          # m/s^2


def optimize_support_radius(
    radius_m,
    thickness_m,
    E=PYREX_E,
    nu=PYREX_NU,
    rho=PYREX_RHO,
    nrefs=5,
    n_points=25,
    r_min_frac=0.2,
    r_max_frac=0.8,
    obstruction_radius=0.0,
):
    """Sweep support radius to find optimal placement minimizing RMS deformation.

    Assembles the stiffness matrix and load vector once, then re-solves for
    each support radius by only changing the constraint DOFs.

    Args:
        radius_m: Mirror radius in meters.
        thickness_m: Mirror thickness in meters.
        E: Young's modulus (Pa).
        nu: Poisson's ratio.
        rho: Density (kg/m^3).
        nrefs: Mesh refinement level.
        n_points: Number of sweep points.
        r_min_frac: Minimum support radius as fraction of mirror radius.
        r_max_frac: Maximum support radius as fraction of mirror radius.

    Returns:
        Tuple of (radii_frac, rms_values, optimal_frac, min_rms):
            radii_frac: Array of support radius fractions tested.
            rms_values: Array of RMS values in nanometers.
            optimal_frac: Optimal support radius as fraction of mirror radius.
            min_rms: Minimum RMS in nanometers.
    """
    # Flexural rigidity
    D = E * thickness_m**3 / (12 * (1 - nu**2))
    # Distributed load (pressure from self-weight)
    q = rho * GRAVITY * thickness_m

    # Create mesh and assemble system once
    mesh = create_circular_mesh(radius_m, nrefs=nrefs)
    K, f, basis = assemble_plate_system(mesh, D, q, nu)

    print(f"Mesh: {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} elements, {K.shape[0]} DOFs")

    # Sweep support radius
    radii_frac = np.linspace(r_min_frac, r_max_frac, n_points)
    rms_values = np.zeros(n_points)

    for i, frac in enumerate(radii_frac):
        support_r = frac * radius_m
        support_dofs = find_support_dofs(basis, mesh, support_r)
        w = solve_plate(K, f, support_dofs)
        rms_values[i] = compute_rms(mesh, basis, w, obstruction_radius)

    # Find optimum
    best_idx = np.argmin(rms_values)
    optimal_frac = radii_frac[best_idx]
    min_rms = rms_values[best_idx]

    return radii_frac, rms_values, optimal_frac, min_rms


def evaluate_single(radius_m, thickness_m, support_frac, E=PYREX_E, nu=PYREX_NU, rho=PYREX_RHO, nrefs=5, obstruction_radius=0.0):
    """Evaluate deformation for a single support radius.

    Returns:
        Tuple of (mesh, basis, w, rms_nm, support_points):
            mesh: The FEM mesh.
            basis: The FEM basis.
            w: Full DOF solution vector.
            rms_nm: RMS deflection in nanometers.
            support_points: (N, 2) array of support point coordinates.
    """
    D = E * thickness_m**3 / (12 * (1 - nu**2))
    q = rho * GRAVITY * thickness_m

    mesh = create_circular_mesh(radius_m, nrefs=nrefs)
    K, f, basis = assemble_plate_system(mesh, D, q, nu)

    print(f"Mesh: {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} elements, {K.shape[0]} DOFs")

    support_r = support_frac * radius_m
    support_dofs = find_support_dofs(basis, mesh, support_r)
    w = solve_plate(K, f, support_dofs)
    rms_nm = compute_rms(mesh, basis, w, obstruction_radius)

    # Compute support point coordinates for visualization
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    support_points = np.column_stack([
        support_r * np.cos(angles),
        support_r * np.sin(angles),
    ])

    return mesh, basis, w, rms_nm, support_points
