# Experiment: Trefoil vs Defocus Dominance Across Support Radii

## Hypothesis

With 3-point support at 120° spacing, mirror deformation is a superposition of
radially symmetric (defocus) and 3-fold symmetric (trefoil) components. Their
relative contribution changes with support radius:

- **Small radii** (r/R < 0.5): The mirror sags like a bowl around the supports.
  The dominant error is radially symmetric — a defocus term (Z₂⁰ = r²). This is
  removable by refocusing the telescope.

- **Large radii** (r/R > 0.6): The supports are near the edge and the mirror
  develops a characteristic 3-lobed deformation (trefoil). This error has 3-fold
  symmetry (Z₃³ = r³cos3θ, Z₃⁻³ = r³sin3θ) and cannot be removed by refocusing.

- **Crossover**: There is a transition region where defocus and trefoil contribute
  roughly equally to the total surface error.

## Consequence for Optimal Support Radius

In standard mode (no refocusing), the optimal radius balances defocus and trefoil
to minimize total RMS — this lands near the classical Grubb radius (0.6789R).

In PLOP mode (refocusing enabled), defocus is removed for free, so the optimizer
only needs to minimize trefoil. Since trefoil grows with radius, the optimal radius
shifts inward (toward ~0.4–0.5R) where trefoil is smaller, even though defocus is
larger there — it doesn't matter because refocusing removes it.

## Zernike Basis Functions

We use unnormalized Zernike-like polynomials for simplicity:

| Term    | Formula           | Physical meaning           |
|---------|-------------------|----------------------------|
| Piston  | 1                 | Constant offset (removed)  |
| Tilt-X  | x                 | Tip (removed)              |
| Tilt-Y  | y                 | Tilt (removed)             |
| Defocus | r²                | Bowl/focus shift (Z₂⁰)    |
| Trefoil | r³cos(3θ), r³sin(3θ) | 3-fold pattern (Z₃³, Z₃⁻³) |

After piston/tilt removal, we fit: `w = a + b·r² + c·r³cos(3θ) + d·r³sin(3θ)`

- Defocus amplitude = |b| (in physical units, scaled by r² range)
- Trefoil amplitude = √(c² + d²)

These are not orthonormalized Zernike coefficients — they're sufficient for showing
relative contributions but shouldn't be compared to Noll-indexed Zernike values.

## Expected Results

1. **Standard vs Refocused RMS**: Two curves with different optimal radii.
   Standard optimum near 0.68R, refocused optimum near 0.4–0.5R.

2. **Defocus fraction** (`1 - refocused_rms²/standard_rms²`): High (~0.7–0.9) at
   small radii, decreasing to low (~0.1–0.3) at large radii. This directly shows
   how much of the variance is removable by refocusing.

3. **Component amplitudes**: Defocus amplitude decreases with radius (monotonic).
   Trefoil amplitude increases with radius (monotonic). They cross somewhere in
   the 0.5–0.6R range.

4. **Deformation maps**: At 0.35R, a smooth bowl shape. At 0.55R, mixed. At 0.75R,
   clear 3-lobed trefoil pattern.

## Verification Criteria

- Standard RMS optimal near 0.68R (Grubb radius)
- Refocused RMS optimal shifts to smaller radius (~0.4–0.5R)
- Defocus fraction high at small radii, low at large radii
- Trefoil amplitude grows monotonically with radius
- Deformation maps show clear visual transition from bowl to trefoil
