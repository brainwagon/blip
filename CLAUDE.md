# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

BLIP is a mirror cell support optimizer for amateur telescope makers. It uses finite element analysis (FEA) to find the optimal placement of three symmetric support points for Newtonian telescope primary mirrors, minimizing gravitational surface deformation.

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
- **rms_calc.py** — Surface error metrics. Computes area-weighted RMS and peak-to-valley (PV) after removing piston/tilt. Supports central obstruction masking.
- **visualize.py** — Matplotlib plots for deformation maps and metric-vs-radius curves.

## Physics Model

Kirchhoff-Love thin plate bending (biharmonic equation D∇⁴w = q) with three point supports at 120° spacing. Default material is Pyrex/borosilicate glass. The tool compares results against the classical Grubb optimal support radius (0.6789R).
