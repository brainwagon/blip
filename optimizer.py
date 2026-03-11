"""Support radius sweep optimizer — assembles K/f once, re-solves per radius."""

import numpy as np
from plate_fem import (
    create_circular_mesh,
    assemble_plate_system,
    find_support_dofs,
    solve_plate,
)
from rms_calc import compute_rms, compute_pv, compute_wavefront_rms, compute_wavefront_pv


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
    metric='rms',
    mode='standard',
):
    """Sweep support radius to find optimal placement.

    Assembles the stiffness matrix and load vector once, then re-solves for
    each support radius by only changing the constraint DOFs.

    In standard mode, both RMS and PV surface deformation (nm) are computed.
    In PLOP mode, RMS and PV wavefront error (waves at 550nm) are computed
    after defocus removal (simulating refocusing).

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
        obstruction_radius: Central obstruction radius in meters.
        metric: 'rms' or 'pv' — which metric to minimize.
        mode: 'standard' or 'plop'.

    Returns:
        Dict with keys:
            radii_frac: Array of support radius fractions tested.
            rms_values: Array of RMS values (nm or waves depending on mode).
            pv_values: Array of PV values (nm or waves depending on mode).
            optimal_frac: Optimal support radius as fraction of mirror radius.
            min_rms: RMS at optimal radius.
            min_pv: PV at optimal radius.
            metric: Which metric was optimized ('rms' or 'pv').
            mode: 'standard' or 'plop'.
    """
    # Flexural rigidity
    D = E * thickness_m**3 / (12 * (1 - nu**2))
    # Distributed load (pressure from self-weight)
    q = rho * GRAVITY * thickness_m

    # Create mesh and assemble system once
    mesh = create_circular_mesh(radius_m, nrefs=nrefs)
    K, f, basis = assemble_plate_system(mesh, D, q, nu)

    print(f"Mesh: {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} elements, {K.shape[0]} DOFs")

    # Select metric functions based on mode
    if mode == 'plop':
        rms_func = compute_wavefront_rms
        pv_func = compute_wavefront_pv
    else:
        rms_func = compute_rms
        pv_func = compute_pv

    # Sweep support radius
    radii_frac = np.linspace(r_min_frac, r_max_frac, n_points)
    rms_values = np.zeros(n_points)
    pv_values = np.zeros(n_points)

    for i, frac in enumerate(radii_frac):
        support_r = frac * radius_m
        support_dofs = find_support_dofs(basis, mesh, support_r)
        w = solve_plate(K, f, support_dofs)
        rms_values[i] = rms_func(mesh, basis, w, obstruction_radius)
        pv_values[i] = pv_func(mesh, basis, w, obstruction_radius)

    # Find optimum based on chosen metric
    target = rms_values if metric == 'rms' else pv_values
    best_idx = np.argmin(target)
    optimal_frac = radii_frac[best_idx]

    return {
        'radii_frac': radii_frac,
        'rms_values': rms_values,
        'pv_values': pv_values,
        'optimal_frac': optimal_frac,
        'min_rms': rms_values[best_idx],
        'min_pv': pv_values[best_idx],
        'metric': metric,
        'mode': mode,
    }


def evaluate_single(radius_m, thickness_m, support_frac, E=PYREX_E, nu=PYREX_NU, rho=PYREX_RHO, nrefs=5, obstruction_radius=0.0, mode='standard'):
    """Evaluate deformation for a single support radius.

    Args:
        radius_m: Mirror radius in meters.
        thickness_m: Mirror thickness in meters.
        support_frac: Support radius as fraction of mirror radius.
        E: Young's modulus (Pa).
        nu: Poisson's ratio.
        rho: Density (kg/m^3).
        nrefs: Mesh refinement level.
        obstruction_radius: Central obstruction radius in meters.
        mode: 'standard' or 'plop'.

    Returns:
        Dict with keys:
            mesh: The FEM mesh.
            basis: The FEM basis.
            w: Full DOF solution vector.
            rms_nm: RMS deflection in nanometers (standard mode).
            pv_nm: Peak-to-valley deflection in nanometers (standard mode).
            wf_rms_waves: RMS wavefront error in waves (PLOP mode, or None).
            wf_pv_waves: PV wavefront error in waves (PLOP mode, or None).
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

    # Always compute surface deformation metrics
    rms_nm = compute_rms(mesh, basis, w, obstruction_radius)
    pv_nm = compute_pv(mesh, basis, w, obstruction_radius)

    # Compute wavefront metrics if in PLOP mode
    wf_rms_waves = None
    wf_pv_waves = None
    if mode == 'plop':
        wf_rms_waves = compute_wavefront_rms(mesh, basis, w, obstruction_radius)
        wf_pv_waves = compute_wavefront_pv(mesh, basis, w, obstruction_radius)

    # Compute support point coordinates for visualization
    angles = np.linspace(0, 2 * np.pi, 3, endpoint=False)
    support_points = np.column_stack([
        support_r * np.cos(angles),
        support_r * np.sin(angles),
    ])

    return {
        'mesh': mesh,
        'basis': basis,
        'w': w,
        'rms_nm': rms_nm,
        'pv_nm': pv_nm,
        'wf_rms_waves': wf_rms_waves,
        'wf_pv_waves': wf_pv_waves,
        'support_points': support_points,
    }
