"""Visualization: deformation map and metric-vs-radius curve."""

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


def plot_metric_vs_radius(results):
    """Plot surface error metric vs support radius fraction.

    Plots both RMS and PV curves, highlighting the optimized metric's optimum.
    Adapts labels and units based on mode (standard: nm, PLOP: waves).

    Args:
        results: Dict from optimize_support_radius() containing radii_frac,
                 rms_values, pv_values, optimal_frac, min_rms, min_pv, metric, mode.
    """
    radii_frac = results['radii_frac']
    rms_values = results['rms_values']
    pv_values = results['pv_values']
    optimal_frac = results['optimal_frac']
    metric = results['metric']
    mode = results.get('mode', 'standard')

    # Units and labels depend on mode
    if mode == 'plop':
        unit = 'waves'
        rms_label = 'RMS Wavefront Error (waves)'
        pv_label = 'P-V Wavefront Error (waves)'
        title_prefix = 'Wavefront Error'
    else:
        unit = 'nm'
        rms_label = 'RMS Surface Deformation (nm)'
        pv_label = 'Peak-to-Valley Deformation (nm)'
        title_prefix = 'Surface Deformation'

    # Compute Strehl ratio for each radius point
    if mode == 'plop':
        rms_waves = np.array(rms_values)
    else:
        rms_waves = np.array(rms_values) / 550.0  # nm to waves
    strehl_values = np.exp(-(2 * np.pi * rms_waves)**2)

    fig, (ax1, ax_strehl) = plt.subplots(2, 1, figsize=(9, 9),
                                          gridspec_kw={'height_ratios': [3, 2]})

    # Plot RMS on left axis
    color_rms = 'tab:blue'
    line_rms, = ax1.plot(radii_frac, rms_values, '-o', color=color_rms,
                         markersize=4, linewidth=1.5, label='RMS')
    ax1.set_xlabel('Support Radius (fraction of mirror radius)')
    ax1.set_ylabel(rms_label, color=color_rms)
    ax1.tick_params(axis='y', labelcolor=color_rms)

    # Plot PV on right axis
    ax2 = ax1.twinx()
    color_pv = 'tab:orange'
    line_pv, = ax2.plot(radii_frac, pv_values, '-s', color=color_pv,
                        markersize=4, linewidth=1.5, label='PV')
    ax2.set_ylabel(pv_label, color=color_pv)
    ax2.tick_params(axis='y', labelcolor=color_pv)

    # Format annotation based on mode
    if mode == 'plop':
        fmt = '.4f'  # more decimals for waves (small numbers)
    else:
        fmt = '.1f'

    # Mark the optimum on the optimized metric's axis
    if metric == 'rms':
        opt_val = results['min_rms']
        ax1.plot(optimal_frac, opt_val, 'r*', markersize=15, zorder=5)
        ax1.axvline(optimal_frac, color='r', linestyle='--', alpha=0.7)
        ann_ax = ax1
        ann_label = f'{optimal_frac:.4f}R\n{opt_val:{fmt}} {unit} RMS'
        ann_y = opt_val
        y_range = max(rms_values) - min(rms_values)
    else:
        opt_val = results['min_pv']
        ax2.plot(optimal_frac, opt_val, 'r*', markersize=15, zorder=5)
        ax2.axvline(optimal_frac, color='r', linestyle='--', alpha=0.7)
        ann_ax = ax2
        ann_label = f'{optimal_frac:.4f}R\n{opt_val:{fmt}} {unit} PV'
        ann_y = opt_val
        y_range = max(pv_values) - min(pv_values)

    ann_ax.annotate(
        ann_label,
        xy=(optimal_frac, ann_y),
        xytext=(optimal_frac + 0.05, ann_y + y_range * 0.15),
        arrowprops=dict(arrowstyle='->', color='red'),
        fontsize=10,
        color='red',
    )

    # Mark the classical 0.6789R reference
    ax1.axvline(0.6789, color='green', linestyle=':', alpha=0.5)

    # Combined legend
    lines = [line_rms, line_pv]
    labels = [l.get_label() for l in lines]
    labels.append('Classical 0.6789R')
    lines.append(plt.Line2D([0], [0], color='green', linestyle=':', alpha=0.5))
    ax1.legend(lines, labels, loc='upper left')

    metric_name = 'RMS' if metric == 'rms' else 'PV'
    ax1.set_title(f'{title_prefix} vs Support Radius (optimizing {metric_name})')
    ax1.grid(True, alpha=0.3)

    # --- Strehl ratio subplot ---
    ax_strehl.plot(radii_frac, strehl_values, '-o', color='tab:purple',
                   markersize=4, linewidth=1.5, label='Strehl ratio')

    # Mark optimum
    opt_idx = np.argmin(np.abs(np.array(radii_frac) - optimal_frac))
    opt_strehl = strehl_values[opt_idx]
    ax_strehl.plot(optimal_frac, opt_strehl, 'r*', markersize=15, zorder=5)
    ax_strehl.axvline(optimal_frac, color='r', linestyle='--', alpha=0.7)
    ax_strehl.annotate(
        f'{optimal_frac:.4f}R\nStrehl {opt_strehl:.4f}',
        xy=(optimal_frac, opt_strehl),
        xytext=(optimal_frac + 0.05, opt_strehl - 0.05),
        arrowprops=dict(arrowstyle='->', color='red'),
        fontsize=10, color='red',
    )

    # Mark classical reference
    ax_strehl.axvline(0.6789, color='green', linestyle=':', alpha=0.5, label='Classical 0.6789R')

    ax_strehl.set_xlabel('Support Radius (fraction of mirror radius)')
    ax_strehl.set_ylabel('Strehl Ratio')
    ax_strehl.set_title('Strehl Ratio vs Support Radius')
    ax_strehl.set_ylim(0, 1.05)
    ax_strehl.legend(loc='lower right')
    ax_strehl.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig
