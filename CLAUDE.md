# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BLIP is a mirror cell support optimizer for amateur telescope makers. It uses finite element analysis (FEA) to find the optimal placement of three symmetric support points for Newtonian telescope primary mirrors. Supports two modes: Standard (minimizes surface deformation in nm) and PLOP (minimizes wavefront error in waves at 550nm after refocusing/defocus removal).

## Running the Tool

```bash
# Basic optimization (finds optimal support radius)
python mirror_cell.py --diameter 200 --thickness 25

# Single radius evaluation
python mirror_cell.py --diameter 200 --thickness 25 --support-radius 0.67

# With central obstruction and PV metric
python mirror_cell.py --diameter 200 --thickness 25 --secondary 50 --metric pv

# Higher accuracy mesh (slower)
python mirror_cell.py --diameter 200 --thickness 25 --nrefs 6

# Suppress plot windows
python mirror_cell.py --diameter 200 --thickness 25 --no-plot

# PLOP mode (wavefront error in waves after refocusing)
python mirror_cell.py --diameter 200 --thickness 25 --mode plop --focal-length 1000

# PLOP mode with single radius evaluation
python mirror_cell.py --diameter 200 --thickness 25 --mode plop --support-radius 0.67
```

There are no tests, no linter, and no build system. This is a pure Python project.

## Dependencies

Install with: `pip install -r requirements.txt`

Core: numpy, scipy, matplotlib, scikit-fem (Morley triangle elements for plate bending FEA).

## Architecture

The pipeline flows: **CLI → FEA → Optimization → Visualization**

- **mirror_cell.py** — Entry point. Argument parsing, validation, orchestration, output formatting.
- **plate_fem.py** — FEA engine. Mesh generation (`create_circular_mesh`), stiffness/load assembly (`assemble_plate_system`), constraint solving (`solve_plate`). Key design: K matrix and load vector are assembled once; only constraint DOFs change per support radius, making sweeps fast.
- **optimizer.py** — Sweeps support radius from 0.2R to 0.8R, calling `solve_plate` at each point. `optimize_support_radius()` returns the optimal radius and metric curves.
- **rms_calc.py** — Surface error metrics. Standard mode: RMS and PV after piston/tilt removal. PLOP mode: wavefront RMS and PV after piston/tilt/defocus removal (2x surface deformation for reflection OPD). Supports central obstruction masking.
- **visualize.py** — Matplotlib plots for deformation maps and metric-vs-radius curves.

## Physics Model

Kirchhoff-Love thin plate bending (biharmonic equation D∇⁴w = q) with three point supports at 120° spacing. Default material is Pyrex/borosilicate glass. The tool compares results against the classical Grubb optimal support radius (0.6789R).
