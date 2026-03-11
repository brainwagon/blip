"""Area-weighted RMS computation with piston/tilt removal."""

import numpy as np


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
    # Extract nodal deflection DOFs (first row of nodal_dofs)
    w_nodal = w[basis.nodal_dofs[0]]
    x = mesh.p[0]
    y = mesh.p[1]

    areas = compute_nodal_areas(mesh)

    # Mask out nodes inside the central obstruction
    r = np.sqrt(x**2 + y**2)
    mask = r >= obstruction_radius
    if not np.any(mask):
        return 0.0

    w_masked = w_nodal[mask]
    x_masked = x[mask]
    y_masked = y[mask]
    areas_masked = areas[mask]

    w_corrected = remove_piston_tilt(w_masked, x_masked, y_masked, areas_masked)

    total_area = np.sum(areas_masked)
    rms = np.sqrt(np.sum(areas_masked * w_corrected**2) / total_area)

    # Convert from meters to nanometers
    return rms * 1e9
