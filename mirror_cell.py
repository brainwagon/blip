#!/usr/bin/env python3
"""Mirror Cell Support Optimizer — CLI entry point.

Analyzes and optimizes the placement of three symmetric support points
for a Newtonian telescope primary mirror using finite element analysis.

Supports two modes:
  Standard Mode — minimizes surface deformation (RMS or PV in nm).
  PLOP Mode — minimizes wavefront error at focus after refocusing,
               expressed in waves (lambda = 550nm).

Usage:
    python mirror_cell.py --diameter 150 --thickness 25
    python mirror_cell.py --diameter 150 --thickness 25 --support-radius 0.68
    python mirror_cell.py --diameter 150 --thickness 25 --metric pv
    python mirror_cell.py --diameter 150 --thickness 25 --mode plop --focal-length 750
"""

import argparse
import sys

import numpy as np

from optimizer import optimize_support_radius, evaluate_single
from visualize import plot_deformation, plot_metric_vs_radius


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Mirror Cell Support Optimizer — find optimal 3-point support placement',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            'Examples:\n'
            '  python mirror_cell.py --diameter 150 --thickness 25\n'
            '  python mirror_cell.py --diameter 150 --thickness 25 --metric pv\n'
            '  python mirror_cell.py --diameter 200 --thickness 30 --support-radius 0.68\n'
            '  python mirror_cell.py --diameter 150 --thickness 25 --mode plop --focal-length 750\n'
        ),
    )
    parser.add_argument('--diameter', type=float, required=True,
                        help='Mirror diameter in mm')
    parser.add_argument('--thickness', type=float, required=True,
                        help='Mirror thickness in mm')
    parser.add_argument('--secondary', type=float, default=None,
                        help='Secondary mirror diameter in mm (central obstruction). '
                             'Deformations inside this radius are excluded from RMS/PV.')
    parser.add_argument('--support-radius', type=float, default=None,
                        help='Support radius as fraction of mirror radius (0.0-1.0). '
                             'If omitted, --optimize is assumed.')
    parser.add_argument('--optimize', action='store_true', default=False,
                        help='Sweep support radius to find optimum (default if no --support-radius)')
    parser.add_argument('--metric', choices=['rms', 'pv'], default='rms',
                        help='Optimization metric: rms (default) or pv (peak-to-valley)')
    parser.add_argument('--mode', choices=['standard', 'plop'], default='standard',
                        help='Analysis mode: standard (surface deformation in nm) or '
                             'plop (wavefront error in waves at 550nm after refocusing)')
    parser.add_argument('--focal-length', type=float, default=None,
                        help='Mirror focal length in mm (used for PLOP mode f-ratio display)')
    parser.add_argument('--no-plot', action='store_true', default=False,
                        help='Suppress plot display')
    parser.add_argument('--n-points', type=int, default=50,
                        help='Number of sweep points for optimization (default: 50)')
    parser.add_argument('--nrefs', type=int, default=5,
                        help='Mesh refinement level (default: 5, use 6 for higher accuracy)')
    return parser.parse_args(argv)


def validate_args(args):
    if args.diameter <= 0:
        print("Error: diameter must be positive", file=sys.stderr)
        sys.exit(1)
    if args.thickness <= 0:
        print("Error: thickness must be positive", file=sys.stderr)
        sys.exit(1)
    if args.thickness >= args.diameter:
        print("Error: thickness must be less than diameter", file=sys.stderr)
        sys.exit(1)
    if args.secondary is not None:
        if args.secondary <= 0:
            print("Error: secondary diameter must be positive", file=sys.stderr)
            sys.exit(1)
        if args.secondary >= args.diameter:
            print("Error: secondary diameter must be less than primary diameter", file=sys.stderr)
            sys.exit(1)
    if args.support_radius is not None:
        if not (0.0 < args.support_radius < 1.0):
            print("Error: support-radius must be between 0 and 1 (exclusive)", file=sys.stderr)
            sys.exit(1)
    if args.mode == 'plop' and args.focal_length is not None:
        if args.focal_length <= 0:
            print("Error: focal-length must be positive", file=sys.stderr)
            sys.exit(1)


def main(argv=None):
    args = parse_args(argv)
    validate_args(args)

    # Convert mm to meters
    radius_m = (args.diameter / 2.0) / 1000.0
    thickness_m = args.thickness / 1000.0
    obstruction_m = (args.secondary / 2.0) / 1000.0 if args.secondary else 0.0

    # Default to optimize mode if no support radius specified
    do_optimize = args.optimize or args.support_radius is None

    metric_name = 'RMS' if args.metric == 'rms' else 'Peak-to-Valley'
    is_plop = args.mode == 'plop'

    print(f"Mirror Cell Support Optimizer")
    print(f"{'='*45}")
    if is_plop:
        print(f"Mode:             PLOP (wavefront error)")
    else:
        print(f"Mode:             Standard (surface deformation)")
    print(f"Mirror diameter:  {args.diameter:.1f} mm")
    print(f"Mirror radius:    {args.diameter/2:.1f} mm")
    print(f"Mirror thickness: {args.thickness:.1f} mm")
    if args.focal_length:
        f_ratio = args.focal_length / args.diameter
        print(f"Focal length:     {args.focal_length:.1f} mm (f/{f_ratio:.1f})")
    if args.secondary:
        print(f"Secondary diam:   {args.secondary:.1f} mm (central obstruction)")
        print(f"  Obstruction:    {args.secondary/args.diameter*100:.1f}% by diameter")
    print(f"Material:         Pyrex (borosilicate glass)")
    print(f"  E  = 63 GPa")
    print(f"  nu = 0.20")
    print(f"  rho = 2230 kg/m^3")
    if is_plop:
        print(f"Wavelength:       550 nm (reference)")
        print(f"Refocusing:       enabled (defocus term removed)")
    if do_optimize:
        print(f"Optimizing:       {metric_name}")
    print()

    if do_optimize:
        print(f"Sweeping support radius from 0.2R to 0.8R ({args.n_points} points)...")
        print()

        results = optimize_support_radius(
            radius_m, thickness_m, nrefs=args.nrefs, n_points=args.n_points,
            obstruction_radius=obstruction_m, metric=args.metric,
            mode=args.mode,
        )

        radii_frac = results['radii_frac']
        rms_values = results['rms_values']
        pv_values = results['pv_values']
        optimal_frac = results['optimal_frac']
        min_rms = results['min_rms']
        min_pv = results['min_pv']

        # Print results table
        if is_plop:
            unit = 'waves'
            rms_hdr = 'RMS (waves)'
            pv_hdr = 'PV (waves)'
            fmt = '10.4f'
        else:
            unit = 'nm'
            rms_hdr = 'RMS (nm)'
            pv_hdr = 'PV (nm)'
            fmt = '10.2f'

        print()
        print(f"{'r/R':>8s}  {rms_hdr:>12s}  {pv_hdr:>12s}")
        print(f"{'-'*8}  {'-'*12}  {'-'*12}")
        for r, rms, pv in zip(radii_frac, rms_values, pv_values):
            marker = ' <-- optimum' if abs(r - optimal_frac) < 1e-10 else ''
            print(f"{r:8.4f}  {rms:{fmt}}  {pv:{fmt}}{marker}")

        print()
        print(f"Optimized metric: {metric_name}")
        print(f"Optimal support radius: {optimal_frac:.4f}R = {optimal_frac * args.diameter/2:.2f} mm")
        if is_plop:
            print(f"  RMS wavefront error at optimum: {min_rms:.4f} waves ({min_rms/4:.4f} Rayleigh)")
            print(f"  P-V wavefront error at optimum: {min_pv:.4f} waves")
            # Optical quality assessment
            if min_rms <= 1/14:
                quality = "Diffraction-limited (Marechal criterion: RMS <= lambda/14)"
            elif min_rms <= 1/10:
                quality = "Good (RMS <= lambda/10)"
            elif min_rms <= 1/4:
                quality = "Acceptable (RMS <= lambda/4)"
            else:
                quality = "Poor (RMS > lambda/4)"
            print(f"  Optical quality: {quality}")
        else:
            print(f"  RMS at optimum: {min_rms:.2f} nm")
            print(f"  PV  at optimum: {min_pv:.2f} nm")

        # Compare with classical reference
        classical_pv = 0.6789
        deviation = abs(optimal_frac - classical_pv) / classical_pv * 100
        print()
        print(f"Classical reference (Grubb, PV-optimal): 0.6789R")
        print(f"Deviation from classical PV-optimal:     {deviation:.1f}%")

        if not args.no_plot:
            import matplotlib.pyplot as plt

            # Also evaluate at optimum for deformation plot
            result = evaluate_single(
                radius_m, thickness_m, optimal_frac, nrefs=args.nrefs,
                obstruction_radius=obstruction_m, mode=args.mode,
            )
            w_nodal = result['w'][result['basis'].nodal_dofs[0]]

            plot_deformation(
                result['mesh'], w_nodal, result['support_points'],
                title=f'Deformation at optimal support ({optimal_frac:.4f}R, {metric_name})',
                obstruction_radius=obstruction_m,
            )
            plot_metric_vs_radius(results)
            plt.show()

    else:
        # Single evaluation mode
        support_frac = args.support_radius
        print(f"Evaluating support radius: {support_frac:.4f}R = {support_frac * args.diameter/2:.2f} mm")
        print()

        result = evaluate_single(
            radius_m, thickness_m, support_frac, nrefs=args.nrefs,
            obstruction_radius=obstruction_m, mode=args.mode,
        )

        w_nodal = result['w'][result['basis'].nodal_dofs[0]]
        print(f"RMS surface deformation:            {result['rms_nm']:.2f} nm")
        print(f"Peak-to-valley surface deformation: {result['pv_nm']:.2f} nm")

        if is_plop:
            print()
            print(f"Wavefront error (after refocusing at 550nm):")
            print(f"  RMS wavefront error: {result['wf_rms_waves']:.4f} waves")
            print(f"  P-V wavefront error: {result['wf_pv_waves']:.4f} waves")

        if not args.no_plot:
            import matplotlib.pyplot as plt
            plot_deformation(
                result['mesh'], w_nodal, result['support_points'],
                title=f'Deformation at {support_frac:.4f}R support radius',
                obstruction_radius=obstruction_m,
            )
            plt.show()


if __name__ == '__main__':
    main()
