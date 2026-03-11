# Mirror Cell Support Optimizer

A free, open-source tool for amateur telescope makers to analyze and optimize the placement of three symmetric support points for a Newtonian telescope primary mirror. It uses finite element analysis to compute gravitational surface deformation and sweeps the support radius to find the placement that minimizes either RMS or peak-to-valley (PV) wavefront error.

## How It Works

The mirror is modeled as a thin circular plate under uniform gravitational loading using Kirchhoff-Love plate bending theory. The governing equation is the biharmonic equation:

```
D * nabla^4(w) = q
```

where `D = E*t^3 / (12*(1-nu^2))` is the flexural rigidity and `q = rho*g*t` is the self-weight pressure.

The FEA uses the **Morley triangle element** via [scikit-fem](https://github.com/kinnala/scikit-fem), a pure-Python finite element library. Three point supports at equal angular spacing (120 degrees apart) are modeled as zero-displacement constraints at the nearest mesh node.

**Key optimization:** The mesh, stiffness matrix K, and load vector f are assembled once. For each candidate support radius in the sweep, only the constraint DOFs change — the system is re-condensed and re-solved. With ~8000 DOFs (nrefs=5), each solve takes milliseconds, making the full sweep fast.

### Metrics

Two surface error metrics are computed at every sweep point:

- **RMS** — Area-weighted root-mean-square of the surface deflection after removing piston and tilt (best-fit plane subtraction via weighted least squares). RMS captures the overall surface quality and correlates well with the Strehl ratio.
- **PV** (Peak-to-Valley) — Maximum minus minimum deflection after piston/tilt removal. PV is sensitive to the single worst point on the surface. The classical Grubb result of 0.6789R was derived to minimize this metric.

The `--metric` flag controls which one drives the optimization. Both are always reported in the output table.

### Central Obstruction

In a Newtonian telescope the secondary mirror casts a shadow on the primary. Deformations inside this shadow don't contribute to the image. The `--secondary` option excludes nodes inside the secondary's radius from both the piston/tilt fit and the RMS/PV calculations, giving a more realistic measure of optical quality.

## Installation

Requires Python 3.10+.

```bash
pip install -r requirements.txt
```

Dependencies: numpy, scipy, matplotlib, scikit-fem.

## Usage

```
python mirror_cell.py --diameter <mm> --thickness <mm> [options]
```

### Options

| Flag | Description |
|---|---|
| `--diameter` | Primary mirror diameter in mm (required) |
| `--thickness` | Mirror thickness in mm (required) |
| `--secondary` | Secondary mirror diameter in mm (central obstruction) |
| `--metric` | Optimization metric: `rms` (default) or `pv` (peak-to-valley) |
| `--support-radius` | Evaluate a single support radius as a fraction of R (0.0-1.0) |
| `--optimize` | Sweep support radius to find the optimum (default if no `--support-radius`) |
| `--n-points` | Number of sweep points (default: 50) |
| `--nrefs` | Mesh refinement level (default: 5, use 6 for higher accuracy) |
| `--no-plot` | Suppress plot windows |

### Examples

Optimize support placement for a 150mm mirror (RMS, default):

```bash
python mirror_cell.py --diameter 150 --thickness 25 --optimize
```

Optimize for peak-to-valley instead:

```bash
python mirror_cell.py --diameter 150 --thickness 25 --optimize --metric pv
```

With a 43mm secondary obstruction:

```bash
python mirror_cell.py --diameter 150 --thickness 25 --secondary 43 --optimize
```

Evaluate a specific support radius:

```bash
python mirror_cell.py --diameter 150 --thickness 25 --support-radius 0.68
```

## Example Results

All examples below use a 150mm diameter, 25mm thick Pyrex mirror with `--nrefs 6 --n-points 100` for higher accuracy.

---

### 1. No obstruction, optimizing RMS (default)

```
Optimizing:       RMS
Optimal support radius: 0.6485R = 48.64 mm
  RMS at optimum: 1.24 nm
  PV  at optimum: 5.29 nm
```

![Surface deformation, no obstruction, RMS-optimal](images/deformation_no_obstruction.png)

The deformation map shows the characteristic three-fold symmetric pattern. The mirror sags between the support points (red) and lifts near the supports and at the center (blue).

![RMS and PV vs support radius, no obstruction, optimizing RMS](images/metric_vs_radius_no_obstruction_rms.png)

Both RMS (blue, left axis) and PV (orange, right axis) are plotted simultaneously. The red dashed line marks the RMS optimum at 0.6485R. The green dotted line shows the classical Grubb reference at 0.6789R. The curves are broadly similar in shape but their minima occur at slightly different radii.

---

### 2. No obstruction, optimizing PV

```
Optimizing:       Peak-to-Valley
Optimal support radius: 0.6364R = 47.73 mm
  RMS at optimum: 1.25 nm
  PV  at optimum: 5.25 nm
```

![Surface deformation, no obstruction, PV-optimal](images/deformation_no_obstruction_pv.png)

![RMS and PV vs support radius, no obstruction, optimizing PV](images/metric_vs_radius_no_obstruction_pv.png)

The PV-optimal radius (0.6364R) is slightly smaller than the RMS-optimal (0.6485R). The PV curve has a broader, flatter minimum than RMS, making it less sensitive to exact placement near the optimum.

---

### 3. With 43mm obstruction, optimizing RMS

```
Secondary diam:   43.0 mm (central obstruction)
  Obstruction:    28.7% by diameter
Optimizing:       RMS
Optimal support radius: 0.6364R = 47.73 mm
  RMS at optimum: 1.29 nm
  PV  at optimum: 5.25 nm
```

![Surface deformation, 43mm obstruction](images/deformation_with_obstruction.png)

The deformation pattern is physically identical — the glass deforms the same way regardless of the secondary shadow. The dashed circle marks the 43mm obstruction boundary; deformations inside it are excluded from RMS/PV calculations.

![RMS and PV vs support radius, 43mm obstruction, optimizing RMS](images/metric_vs_radius_with_obstruction_rms.png)

---

### 4. With 43mm obstruction, optimizing PV

```
Optimizing:       Peak-to-Valley
Optimal support radius: 0.6364R = 47.73 mm
  RMS at optimum: 1.29 nm
  PV  at optimum: 5.25 nm
```

![Surface deformation, 43mm obstruction, PV-optimal](images/deformation_with_obstruction_pv.png)

![RMS and PV vs support radius, 43mm obstruction, optimizing PV](images/metric_vs_radius_with_obstruction_pv.png)

---

### Summary comparison

| Configuration | Metric | Optimal r/R | RMS (nm) | PV (nm) |
|---|---|---|---|---|
| No obstruction | RMS | 0.6485 | 1.24 | 5.29 |
| No obstruction | PV | 0.6364 | 1.25 | 5.25 |
| 43mm obstruction | RMS | 0.6364 | 1.29 | 5.25 |
| 43mm obstruction | PV | 0.6364 | 1.29 | 5.25 |

Key observations:

- **RMS vs PV optimum:** Without obstruction, the RMS-optimal radius (0.649R) is slightly larger than the PV-optimal (0.636R). The difference is small (~1% of mirror radius) and both are near the classical 0.6789R.
- **Effect of obstruction on RMS:** Excluding the central zone raises the minimum RMS slightly (1.24 to 1.29 nm) because the well-behaved central nodes that were pulling the average down are now excluded.
- **Effect of obstruction on PV:** PV is unchanged at 5.25 nm — the peak and valley both occur in the outer annulus, outside the obstruction.
- **Optimal radius shifts inward with obstruction:** The RMS optimum moves from 0.649R to 0.636R because the piston/tilt fit, now computed only over the outer annulus, changes the residual shape that RMS measures.
- **Flat minimum:** The minimum is broad in all cases, meaning placement errors of a few percent have minimal impact on optical quality.

## File Structure

```
mirror_cell.py    — CLI entry point, argument parsing, output formatting
plate_fem.py      — FEA engine: mesh, stiffness/load assembly, support constraints, solve
rms_calc.py       — Surface error metrics (RMS, PV) with piston/tilt removal and obstruction masking
optimizer.py      — Support radius sweep (reuses single K/f assembly)
visualize.py      — Matplotlib: deformation map, dual-axis metric-vs-radius curve
requirements.txt  — Python dependencies
```

## Material Properties

Default material is Pyrex (borosilicate glass):

| Property | Value |
|---|---|
| Young's modulus (E) | 63 GPa |
| Poisson's ratio (nu) | 0.20 |
| Density (rho) | 2230 kg/m^3 |

## Notes

- The mesh is generated from `MeshTri.init_circle()` and scaled to the mirror radius. Support points snap to the nearest mesh node, which introduces small quantization steps visible in the sweep curves. Use `--nrefs 6` for finer resolution (at the cost of ~4x more DOFs and longer solve times).
- The classical Grubb optimal radius of 0.6789R was derived analytically for a uniformly loaded thin plate. The FEA results (~0.64-0.65R) are close but not identical, partly due to mesh discretization and partly because the analytical derivation uses different assumptions about the piston/tilt reference.
- Deflections for typical amateur mirrors (150-300mm) are on the order of nanometers — far below the wavelength of visible light (~550 nm). Three-point support is adequate for mirrors up to roughly 300mm; larger mirrors typically need more support points.
