#!/usr/bin/env python3
"""Experiment: Trefoil vs Defocus Dominance Across Support Radii.

Sweeps support radius and decomposes mirror deformation into defocus (r²)
and trefoil (r³cos3θ, r³sin3θ) components to show why refocusing shifts
the optimal support radius inward from the Grubb radius.

Usage:
    python experiment_trefoil.py --diameter 200 --thickness 25
    python experiment_trefoil.py --diameter 200 --thickness 25 -o trefoil_experiment.png
"""

import argparse
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from plate_fem import (
    create_circular_mesh,
    assemble_plate_system,
    find_support_dofs,
    solve_plate,
)
from rms_calc import compute_nodal_areas, remove_piston_tilt

# Material properties (Pyrex)
PYREX_E = 63e9
PYREX_NU = 0.20
PYREX_RHO = 2230.0
GRAVITY = 9.81


def zernike_decomposition(w_pt, x, y, areas):
    """Decompose piston/tilt-removed deflection into defocus and trefoil.

    Fits: w = a + b*r² + c*r³*cos(3θ) + d*r³*sin(3θ)

    Args:
        w_pt: Deflection with piston/tilt already removed.
        x, y: Node coordinates.
        areas: Nodal area weights.

    Returns:
        Dict with defocus_coeff (b), trefoil_cos (c), trefoil_sin (d),
        trefoil_amplitude (sqrt(c²+d²)), and the defocus/trefoil
        contributions in meters RMS.
    """
    r2 = x**2 + y**2
    r = np.sqrt(r2)
    theta = np.arctan2(y, x)

    r3_cos3 = r**3 * np.cos(3 * theta)
    r3_sin3 = r**3 * np.sin(3 * theta)

    # Weighted least squares fit
    A = np.column_stack([np.ones_like(x), r2, r3_cos3, r3_sin3])
    W = np.diag(areas)
    AtW = A.T @ W
    coeffs = np.linalg.solve(AtW @ A, AtW @ w_pt)

    a, b, c, d = coeffs

    # Compute RMS contribution of each component
    total_area = np.sum(areas)
    defocus_signal = b * r2
    trefoil_signal = c * r3_cos3 + d * r3_sin3

    defocus_rms = np.sqrt(np.sum(areas * defocus_signal**2) / total_area)
    trefoil_rms = np.sqrt(np.sum(areas * trefoil_signal**2) / total_area)

    return {
        'defocus_coeff': b,
        'trefoil_cos': c,
        'trefoil_sin': d,
        'trefoil_amplitude': np.sqrt(c**2 + d**2),
        'defocus_rms_m': defocus_rms,
        'trefoil_rms_m': trefoil_rms,
    }


def run_experiment(radius_m, thickness_m, nrefs=5, n_points=50):
    """Sweep support radius and compute all metrics.

    Returns:
        Dict with arrays of results indexed by sweep point.
    """
    D = PYREX_E * thickness_m**3 / (12 * (1 - PYREX_NU**2))
    q = PYREX_RHO * GRAVITY * thickness_m

    mesh = create_circular_mesh(radius_m, nrefs=nrefs)
    K, f, basis = assemble_plate_system(mesh, D, q, PYREX_NU)

    print(f"Mesh: {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} elements, {K.shape[0]} DOFs")

    x = mesh.p[0]
    y = mesh.p[1]
    areas = compute_nodal_areas(mesh)

    radii_frac = np.linspace(0.2, 0.8, n_points)
    standard_rms = np.zeros(n_points)
    refocused_rms = np.zeros(n_points)
    defocus_rms_nm = np.zeros(n_points)
    trefoil_rms_nm = np.zeros(n_points)
    defocus_fraction = np.zeros(n_points)

    # Store solutions for deformation maps at selected radii
    map_fracs = [0.35, 0.55, 0.75]
    map_solutions = {}

    for i, frac in enumerate(radii_frac):
        support_r = frac * radius_m
        support_dofs = find_support_dofs(basis, mesh, support_r)
        w = solve_plate(K, f, support_dofs)
        w_nodal = w[basis.nodal_dofs[0]]

        # Standard RMS (piston/tilt removed)
        w_pt = remove_piston_tilt(w_nodal, x, y, areas)
        std_rms_m = np.sqrt(np.sum(areas * w_pt**2) / np.sum(areas))
        standard_rms[i] = std_rms_m * 1e9  # nm

        # Refocused RMS (piston/tilt/defocus removed)
        r2 = x**2 + y**2
        A_def = np.column_stack([np.ones_like(x), x, y, r2])
        W = np.diag(areas)
        AtW = A_def.T @ W
        coeffs = np.linalg.solve(AtW @ A_def, AtW @ w_nodal)
        w_refoc = w_nodal - A_def @ coeffs
        ref_rms_m = np.sqrt(np.sum(areas * w_refoc**2) / np.sum(areas))
        refocused_rms[i] = ref_rms_m * 1e9  # nm

        # Zernike decomposition (on piston/tilt-removed data)
        decomp = zernike_decomposition(w_pt, x, y, areas)
        defocus_rms_nm[i] = decomp['defocus_rms_m'] * 1e9
        trefoil_rms_nm[i] = decomp['trefoil_rms_m'] * 1e9

        # Defocus fraction of variance
        if std_rms_m > 0:
            defocus_fraction[i] = 1.0 - (ref_rms_m**2 / std_rms_m**2)
        else:
            defocus_fraction[i] = 0.0

        # Save solutions for deformation maps
        for mf in map_fracs:
            if abs(frac - mf) < (0.6 / n_points) / 2 and mf not in map_solutions:
                map_solutions[mf] = {
                    'w_nodal': w_nodal.copy(),
                    'w_pt': w_pt.copy(),
                    'frac': frac,
                }

    return {
        'radii_frac': radii_frac,
        'standard_rms': standard_rms,
        'refocused_rms': refocused_rms,
        'defocus_rms_nm': defocus_rms_nm,
        'trefoil_rms_nm': trefoil_rms_nm,
        'defocus_fraction': defocus_fraction,
        'map_solutions': map_solutions,
        'mesh': mesh,
        'radius_m': radius_m,
    }


def plot_results(results, output_path=None):
    """Create 4-panel figure."""
    radii = results['radii_frac']
    mesh = results['mesh']
    x = mesh.p[0]
    y = mesh.p[1]

    fig = plt.figure(figsize=(10, 14))
    gs = GridSpec(4, 3, figure=fig, height_ratios=[1, 1, 1, 1],
                  hspace=0.35, wspace=0.3)

    # --- Panel 1: Standard RMS vs Refocused RMS ---
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(radii, results['standard_rms'], 'b-', linewidth=2, label='Standard RMS')
    ax1.plot(radii, results['refocused_rms'], 'r-', linewidth=2, label='Refocused RMS')

    # Mark optima
    std_opt_idx = np.argmin(results['standard_rms'])
    ref_opt_idx = np.argmin(results['refocused_rms'])
    ax1.axvline(radii[std_opt_idx], color='b', linestyle='--', alpha=0.5)
    ax1.axvline(radii[ref_opt_idx], color='r', linestyle='--', alpha=0.5)
    ax1.plot(radii[std_opt_idx], results['standard_rms'][std_opt_idx],
             'bv', markersize=10, label=f'Standard opt: {radii[std_opt_idx]:.3f}R')
    ax1.plot(radii[ref_opt_idx], results['refocused_rms'][ref_opt_idx],
             'rv', markersize=10, label=f'Refocused opt: {radii[ref_opt_idx]:.3f}R')

    # Grubb reference
    ax1.axvline(0.6789, color='gray', linestyle=':', alpha=0.7, label='Grubb (0.6789R)')

    ax1.set_xlabel('Support radius (r/R)')
    ax1.set_ylabel('RMS surface error (nm)')
    ax1.set_title('Standard vs Refocused RMS')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Defocus fraction ---
    ax2 = fig.add_subplot(gs[1, :])
    ax2.plot(radii, results['defocus_fraction'], 'g-', linewidth=2)
    ax2.axhline(0.5, color='gray', linestyle=':', alpha=0.5, label='50% line')
    ax2.fill_between(radii, results['defocus_fraction'], alpha=0.15, color='g')
    ax2.set_xlabel('Support radius (r/R)')
    ax2.set_ylabel('Defocus fraction of variance')
    ax2.set_title('Defocus Fraction (removable by refocusing)')
    ax2.set_ylim(-0.05, 1.05)
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # --- Panel 3: Component amplitudes ---
    ax3 = fig.add_subplot(gs[2, :])
    ax3.plot(radii, results['defocus_rms_nm'], 'b-', linewidth=2, label='Defocus RMS')
    ax3.plot(radii, results['trefoil_rms_nm'], 'r-', linewidth=2, label='Trefoil RMS')

    # Find crossover
    diff = results['defocus_rms_nm'] - results['trefoil_rms_nm']
    sign_changes = np.where(np.diff(np.sign(diff)))[0]
    if len(sign_changes) > 0:
        # Linear interpolation for crossover point
        idx = sign_changes[0]
        x_cross = radii[idx] + (radii[idx+1] - radii[idx]) * abs(diff[idx]) / (abs(diff[idx]) + abs(diff[idx+1]))
        y_cross = np.interp(x_cross, radii, results['defocus_rms_nm'])
        ax3.plot(x_cross, y_cross, 'ko', markersize=8, label=f'Crossover: {x_cross:.3f}R')

    ax3.set_xlabel('Support radius (r/R)')
    ax3.set_ylabel('Component RMS (nm)')
    ax3.set_title('Zernike Component Amplitudes')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # --- Panel 4: Deformation maps ---
    map_fracs = sorted(results['map_solutions'].keys())
    triangulation = matplotlib.tri.Triangulation(x, y, mesh.t.T)

    for j, mf in enumerate(map_fracs):
        ax = fig.add_subplot(gs[3, j])
        sol = results['map_solutions'][mf]
        w_nm = sol['w_pt'] * 1e9  # piston/tilt removed, in nm
        vmax = max(abs(w_nm.max()), abs(w_nm.min()))
        tc = ax.tripcolor(triangulation, w_nm, shading='gouraud',
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        ax.set_aspect('equal')
        ax.set_title(f'r/R = {sol["frac"]:.3f}', fontsize=10)
        ax.set_xlim(-results['radius_m'] * 1.05, results['radius_m'] * 1.05)
        ax.set_ylim(-results['radius_m'] * 1.05, results['radius_m'] * 1.05)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(tc, ax=ax, shrink=0.8, label='nm')

    fig.suptitle('Trefoil vs Defocus Dominance Across Support Radii',
                 fontsize=14, fontweight='bold', y=0.98)

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    else:
        plt.show()

    return fig


def print_summary(results):
    """Print key results."""
    radii = results['radii_frac']
    std_opt_idx = np.argmin(results['standard_rms'])
    ref_opt_idx = np.argmin(results['refocused_rms'])

    print()
    print("Results")
    print("=" * 55)
    print(f"Standard RMS optimal radius:  {radii[std_opt_idx]:.4f}R "
          f"({results['standard_rms'][std_opt_idx]:.2f} nm)")
    print(f"Refocused RMS optimal radius: {radii[ref_opt_idx]:.4f}R "
          f"({results['refocused_rms'][ref_opt_idx]:.2f} nm)")
    print(f"Classical Grubb radius:       0.6789R")
    print()

    # Print table
    print(f"{'r/R':>8s}  {'Std RMS':>10s}  {'Ref RMS':>10s}  "
          f"{'Defoc RMS':>10s}  {'Trefoil':>10s}  {'Defoc frac':>10s}")
    print(f"{'':>8s}  {'(nm)':>10s}  {'(nm)':>10s}  "
          f"{'(nm)':>10s}  {'(nm)':>10s}  {'':>10s}")
    print("-" * 65)
    # Print every 5th point plus the optima
    for i in range(len(radii)):
        if i % 5 == 0 or i == std_opt_idx or i == ref_opt_idx:
            marker = ''
            if i == std_opt_idx:
                marker = ' <-- std opt'
            elif i == ref_opt_idx:
                marker = ' <-- ref opt'
            print(f"{radii[i]:8.4f}  {results['standard_rms'][i]:10.2f}  "
                  f"{results['refocused_rms'][i]:10.2f}  "
                  f"{results['defocus_rms_nm'][i]:10.2f}  "
                  f"{results['trefoil_rms_nm'][i]:10.2f}  "
                  f"{results['defocus_fraction'][i]:10.4f}{marker}")


def main():
    parser = argparse.ArgumentParser(
        description='Trefoil vs Defocus experiment')
    parser.add_argument('--diameter', type=float, required=True,
                        help='Mirror diameter in mm')
    parser.add_argument('--thickness', type=float, required=True,
                        help='Mirror thickness in mm')
    parser.add_argument('--nrefs', type=int, default=5,
                        help='Mesh refinement level (default: 5)')
    parser.add_argument('--n-points', type=int, default=50,
                        help='Number of sweep points (default: 50)')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Save plot to file (e.g. trefoil_experiment.png)')
    parser.add_argument('--no-plot', action='store_true', default=False,
                        help='Suppress plot display')
    args = parser.parse_args()

    if args.diameter <= 0 or args.thickness <= 0:
        print("Error: diameter and thickness must be positive", file=sys.stderr)
        sys.exit(1)

    radius_m = (args.diameter / 2.0) / 1000.0
    thickness_m = args.thickness / 1000.0

    print(f"Trefoil vs Defocus Experiment")
    print(f"{'=' * 45}")
    print(f"Mirror diameter:  {args.diameter:.1f} mm")
    print(f"Mirror thickness: {args.thickness:.1f} mm")
    print()

    results = run_experiment(radius_m, thickness_m, nrefs=args.nrefs,
                             n_points=args.n_points)
    print_summary(results)

    if not args.no_plot or args.output:
        if args.output:
            matplotlib.use('Agg')
        plot_results(results, output_path=args.output)


if __name__ == '__main__':
    main()
