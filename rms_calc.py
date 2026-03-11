"""Surface error metrics with piston/tilt removal and obstruction masking.

Includes both surface deformation metrics (Standard Mode) and wavefront
error metrics (PLOP Mode) with defocus removal for refocusing simulation.
"""

import numpy as np

# Reference wavelength for wavefront error in waves
WAVELENGTH_NM = 550.0  # nm (visible light)


def compute_nodal_areas(mesh):
    """Compute area associated with each node (1/3 of each adjacent triangle).

    Args:
        mesh: Triangular mesh with p (nodes) and t (elements) arrays.

    Returns:
        Array of nodal areas, shape (n_nodes,).
    """
    n_nodes = mesh.p.shape[1]
    areas = np.zeros(n_nodes)

    # Compute triangle areas and distribute 1/3 to each vertex
    for i in range(mesh.t.shape[1]):
        tri_nodes = mesh.t[:, i]
        coords = mesh.p[:, tri_nodes].T  # (3, 2)
        # Signed area via cross product
        v1 = coords[1] - coords[0]
        v2 = coords[2] - coords[0]
        tri_area = 0.5 * abs(v1[0] * v2[1] - v1[1] * v2[0])
        for node in tri_nodes:
            areas[node] += tri_area / 3.0

    return areas


def remove_piston_tilt(w, x, y, areas):
    """Remove best-fit plane (piston + tilt) using area-weighted least squares.

    Fits w = a + b*x + c*y and subtracts it.

    Args:
        w: Deflection values at nodes.
        x, y: Node coordinates.
        areas: Nodal area weights.

    Returns:
        Deflection with piston and tilt removed.
    """
    # Weighted least squares: minimize sum(areas * (w - a - b*x - c*y)^2)
    A = np.column_stack([np.ones_like(x), x, y])
    W = np.diag(areas)
    # Normal equations: (A^T W A) c = A^T W w
    AtW = A.T @ W
    coeffs = np.linalg.solve(AtW @ A, AtW @ w)
    plane = A @ coeffs
    return w - plane


def remove_piston_tilt_defocus(w, x, y, areas):
    """Remove best-fit piston, tilt, and defocus using area-weighted least squares.

    Fits w = a + b*x + c*y + d*(x^2 + y^2) and subtracts it.
    The defocus term (r^2) simulates the observer refocusing the telescope
    to minimize wavefront error at the focal plane.

    Args:
        w: Deflection values at nodes.
        x, y: Node coordinates.
        areas: Nodal area weights.

    Returns:
        Deflection with piston, tilt, and defocus removed.
    """
    r2 = x**2 + y**2
    A = np.column_stack([np.ones_like(x), x, y, r2])
    W = np.diag(areas)
    AtW = A.T @ W
    coeffs = np.linalg.solve(AtW @ A, AtW @ w)
    fit = A @ coeffs
    return w - fit


def _corrected_deflection(mesh, basis, w, obstruction_radius=0.0):
    """Extract masked, piston/tilt-corrected nodal deflection.

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.

    Returns:
        Tuple of (w_corrected, areas_masked) — corrected deflections and
        corresponding nodal areas for unobstructed nodes. Returns (None, None)
        if all nodes are masked.
    """
    w_nodal = w[basis.nodal_dofs[0]]
    x = mesh.p[0]
    y = mesh.p[1]

    areas = compute_nodal_areas(mesh)

    # Mask out nodes inside the central obstruction
    r = np.sqrt(x**2 + y**2)
    mask = r >= obstruction_radius
    if not np.any(mask):
        return None, None

    w_masked = w_nodal[mask]
    x_masked = x[mask]
    y_masked = y[mask]
    areas_masked = areas[mask]

    w_corrected = remove_piston_tilt(w_masked, x_masked, y_masked, areas_masked)
    return w_corrected, areas_masked


def _wavefront_corrected(mesh, basis, w, obstruction_radius=0.0):
    """Extract masked, piston/tilt/defocus-corrected wavefront error.

    Converts surface deformation to wavefront error (2x for reflection)
    and removes piston, tilt, and defocus to simulate refocusing.

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.

    Returns:
        Tuple of (wf_corrected, areas_masked) — corrected wavefront error
        in meters and corresponding nodal areas. Returns (None, None)
        if all nodes are masked.
    """
    w_nodal = w[basis.nodal_dofs[0]]
    x = mesh.p[0]
    y = mesh.p[1]

    areas = compute_nodal_areas(mesh)

    # Mask out nodes inside the central obstruction
    r = np.sqrt(x**2 + y**2)
    mask = r >= obstruction_radius
    if not np.any(mask):
        return None, None

    # Wavefront error = 2 * surface deformation (reflection doubles OPD)
    wf_masked = 2.0 * w_nodal[mask]
    x_masked = x[mask]
    y_masked = y[mask]
    areas_masked = areas[mask]

    # Remove piston, tilt, and defocus (simulates refocusing)
    wf_corrected = remove_piston_tilt_defocus(wf_masked, x_masked, y_masked, areas_masked)
    return wf_corrected, areas_masked


def compute_rms(mesh, basis, w, obstruction_radius=0.0):
    """Compute area-weighted RMS of deflection after piston/tilt removal.

    Nodes inside the central obstruction (secondary mirror shadow) are
    excluded from both the piston/tilt fit and the RMS computation.

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.
            Nodes with r < obstruction_radius are excluded. Default 0 (none).

    Returns:
        RMS deflection in nanometers.
    """
    w_corrected, areas_masked = _corrected_deflection(mesh, basis, w, obstruction_radius)
    if w_corrected is None:
        return 0.0

    total_area = np.sum(areas_masked)
    rms = np.sqrt(np.sum(areas_masked * w_corrected**2) / total_area)

    # Convert from meters to nanometers
    return rms * 1e9


def compute_pv(mesh, basis, w, obstruction_radius=0.0):
    """Compute peak-to-valley deflection after piston/tilt removal.

    Nodes inside the central obstruction (secondary mirror shadow) are
    excluded from both the piston/tilt fit and the PV computation.

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.
            Nodes with r < obstruction_radius are excluded. Default 0 (none).

    Returns:
        Peak-to-valley deflection in nanometers.
    """
    w_corrected, areas_masked = _corrected_deflection(mesh, basis, w, obstruction_radius)
    if w_corrected is None:
        return 0.0

    pv = np.max(w_corrected) - np.min(w_corrected)

    # Convert from meters to nanometers
    return pv * 1e9


def compute_wavefront_rms(mesh, basis, w, obstruction_radius=0.0):
    """Compute area-weighted RMS wavefront error after piston/tilt/defocus removal.

    Converts surface deformation to wavefront error (2x for reflection),
    removes piston, tilt, and defocus (simulating refocusing), then computes
    area-weighted RMS. Result is in waves (lambda = 550nm).

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.

    Returns:
        RMS wavefront error in waves (lambda = 550nm).
    """
    wf_corrected, areas_masked = _wavefront_corrected(mesh, basis, w, obstruction_radius)
    if wf_corrected is None:
        return 0.0

    total_area = np.sum(areas_masked)
    rms_m = np.sqrt(np.sum(areas_masked * wf_corrected**2) / total_area)

    # Convert from meters to waves
    return rms_m / (WAVELENGTH_NM * 1e-9)


def compute_wavefront_pv(mesh, basis, w, obstruction_radius=0.0):
    """Compute peak-to-valley wavefront error after piston/tilt/defocus removal.

    Converts surface deformation to wavefront error (2x for reflection),
    removes piston, tilt, and defocus (simulating refocusing), then computes
    P-V. Result is in waves (lambda = 550nm).

    Args:
        mesh: Triangular mesh.
        basis: FEM basis (Morley element).
        w: Full DOF solution vector.
        obstruction_radius: Radius of central obstruction in meters.

    Returns:
        Peak-to-valley wavefront error in waves (lambda = 550nm).
    """
    wf_corrected, areas_masked = _wavefront_corrected(mesh, basis, w, obstruction_radius)
    if wf_corrected is None:
        return 0.0

    pv_m = np.max(wf_corrected) - np.min(wf_corrected)

    # Convert from meters to waves
    return pv_m / (WAVELENGTH_NM * 1e-9)
