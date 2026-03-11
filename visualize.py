"""Visualization: deformation map and RMS-vs-radius curve."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def plot_deformation(mesh, w_nodal, support_points, title=None, obstruction_radius=0.0):
    """Plot the mirror surface deformation as a filled contour map.

    Args:
        mesh: Triangular mesh.
        w_nodal: Nodal deflection values (meters).
        support_points: (N, 2) array of support point coordinates.
        title: Optional plot title.
        obstruction_radius: Central obstruction radius in meters (0 = none).
    """
    x = mesh.p[0]
    y = mesh.p[1]
    triangles = mesh.t.T  # (n_elements, 3)

    tri = mtri.Triangulation(x, y, triangles)

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    # Convert to micrometers for display
    w_um = w_nodal * 1e6
    tc = ax.tripcolor(tri, w_um, shading='gouraud', cmap='RdBu_r')
    cb = fig.colorbar(tc, ax=ax, label='Deflection (µm)')

    # Mark support points
    ax.plot(support_points[:, 0], support_points[:, 1], 'k^', markersize=10,
            markeredgewidth=2, label='Support points')

    # Draw mirror outline
    theta = np.linspace(0, 2 * np.pi, 100)
    r_max = np.max(np.sqrt(x**2 + y**2))
    ax.plot(r_max * np.cos(theta), r_max * np.sin(theta), 'k-', linewidth=1)

    # Draw central obstruction
    if obstruction_radius > 0:
        ax.plot(obstruction_radius * np.cos(theta), obstruction_radius * np.sin(theta),
                'k--', linewidth=1.5, label='Secondary obstruction')

    ax.set_aspect('equal')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.legend(loc='upper right')
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Mirror Surface Deformation')

    plt.tight_layout()
    return fig


def plot_rms_vs_radius(radii_frac, rms_values, optimal_frac, min_rms):
    """Plot RMS deformation vs support radius fraction.

    Args:
        radii_frac: Array of support radius fractions (r/R).
        rms_values: Array of RMS values in nanometers.
        optimal_frac: Optimal support radius fraction.
        min_rms: Minimum RMS in nanometers.
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(radii_frac, rms_values, 'b-o', markersize=4, linewidth=1.5)
    ax.axvline(optimal_frac, color='r', linestyle='--', alpha=0.7,
               label=f'Optimum: {optimal_frac:.4f}R')
    ax.plot(optimal_frac, min_rms, 'r*', markersize=15, zorder=5)

    # Annotate the optimum
    ax.annotate(
        f'{optimal_frac:.4f}R\n{min_rms:.1f} nm RMS',
        xy=(optimal_frac, min_rms),
        xytext=(optimal_frac + 0.05, min_rms + (max(rms_values) - min(rms_values)) * 0.15),
        arrowprops=dict(arrowstyle='->', color='red'),
        fontsize=10,
        color='red',
    )

    # Mark the classical 0.6789R reference (PV-optimal)
    ax.axvline(0.6789, color='green', linestyle=':', alpha=0.5,
               label='Classical 0.6789R (PV-optimal)')

    ax.set_xlabel('Support Radius (fraction of mirror radius)')
    ax.set_ylabel('RMS Surface Deformation (nm)')
    ax.set_title('RMS Deformation vs Support Radius')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig
