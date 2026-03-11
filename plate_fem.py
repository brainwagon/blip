"""FEA engine: mesh generation, assembly, support constraints, and solve."""

import numpy as np
from skfem import *
from skfem.models import laplace
from skfem.helpers import ddot, dd


def create_circular_mesh(radius, nrefs=5):
    """Create a circular triangular mesh scaled to the given radius.

    Args:
        radius: Mirror radius in meters.
        nrefs: Number of uniform refinements (5 gives ~3000 elements).

    Returns:
        MeshTri instance.
    """
    mesh = MeshTri.init_circle(nrefs=nrefs)
    # init_circle creates a unit circle; scale to physical radius
    mesh = mesh.scaled(radius)
    return mesh


def assemble_plate_system(mesh, D, q, nu):
    """Assemble the Kirchhoff-Love plate bending system using Morley elements.

    Uses the biharmonic equation: D * nabla^4(w) = q
    where D = E*t^3 / (12*(1-nu^2)), q = rho*g*t.

    The weak form uses ddot(C : dd(u), dd(v)) where C is the bending
    stiffness tensor for an isotropic plate.

    Args:
        mesh: Triangular mesh.
        D: Flexural rigidity (Pa * m^3).
        q: Distributed load (Pa), typically rho * g * t.
        nu: Poisson's ratio.

    Returns:
        Tuple of (K, f, basis) — stiffness matrix, load vector, basis object.
    """
    element = ElementTriMorley()
    basis = CellBasis(mesh, element)

    # Bending stiffness bilinear form
    # For isotropic plate: M = D * [(1-nu)*kappa + nu*tr(kappa)*I]
    # Energy = integral of M:kappa = D * [(1-nu)*ddot(dd(u),dd(v)) + nu*div(grad(u))*div(grad(v))]
    @BilinearForm
    def bilinear(u, v, w):
        # dd(u) is the Hessian matrix of u
        # For isotropic plate bending:
        # a(u,v) = D * integral[ (1-nu)*ddot(dd(u), dd(v)) + nu*trace(dd(u))*trace(dd(v)) ]
        d2u = dd(u)
        d2v = dd(v)
        return D * ((1 - nu) * ddot(d2u, d2v) + nu * (d2u[0, 0] + d2u[1, 1]) * (d2v[0, 0] + d2v[1, 1]))

    @LinearForm
    def load(v, w):
        return q * v

    K = bilinear.assemble(basis)
    f = load.assemble(basis)

    return K, f, basis


def find_support_dofs(basis, mesh, support_radius, num_supports=3):
    """Find DOF indices for support points at equal angular intervals.

    Places num_supports points equally spaced in angle at the given radius,
    starting at angle=0.

    Args:
        basis: FEM basis object.
        mesh: Triangular mesh.
        support_radius: Radial distance of supports from center (meters).
        num_supports: Number of support points (default 3).

    Returns:
        Array of DOF indices corresponding to the deflection DOF at each support.
    """
    nodes = mesh.p.T  # shape (n_nodes, 2)
    angles = np.linspace(0, 2 * np.pi, num_supports, endpoint=False)
    support_points = np.column_stack([
        support_radius * np.cos(angles),
        support_radius * np.sin(angles),
    ])

    dof_indices = []
    for pt in support_points:
        dists = np.linalg.norm(nodes - pt, axis=1)
        nearest_node = np.argmin(dists)
        # Morley element: nodal_dofs[0] are the deflection DOFs
        dof_idx = basis.nodal_dofs[0, nearest_node]
        dof_indices.append(dof_idx)

    return np.array(dof_indices, dtype=int)


def solve_plate(K, f, support_dofs):
    """Solve the plate bending problem with zero-deflection constraints at supports.

    Args:
        K: Stiffness matrix.
        f: Load vector.
        support_dofs: DOF indices where deflection is constrained to zero.

    Returns:
        Full deflection vector (all DOFs).
    """
    return solve(*condense(K, f, D=support_dofs))
