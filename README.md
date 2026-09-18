# RCWA Convergence Study

A small, self-contained C++ program that simulates light passing through a
diffraction grating using a simplified version of **Rigorous Coupled-Wave
Analysis (RCWA)**, and shows how the answer becomes more accurate as you
give the simulation more "resolution."

No external libraries are required — just a C++17 compiler.

## What is RCWA, in plain terms?

Imagine shining a laser through a striped grating (alternating strips of
glass and air, repeating like a fence). Some light goes straight through,
some gets bent off at angles ("diffraction orders"). RCWA is a standard
method in optics for calculating exactly how much light ends up in each
of those directions.

To do the math, RCWA describes the grating's pattern as a sum of waves
(a **Fourier series**) instead of a sharp on/off pattern. The more waves
you include, the more faithfully the simulation represents the real
grating — but each extra wave also makes the underlying matrix math
bigger and slower. Finding the sweet spot between *accurate* and *fast*
is called a **convergence study**, and that's what this program does.

## What this code actually does

`RCWASolver` in [`src/rcwa_convergence.cpp`](src/rcwa_convergence.cpp):

1. Builds the Fourier coefficients of the grating's permittivity
   (its optical "material pattern") for a chosen number of harmonics.
2. Assembles that into a Toeplitz matrix, a standard step in RCWA that
   encodes how the different diffraction orders couple to each other.
3. Solves a linear system (via straightforward Gauss-Jordan elimination)
   to find the field of each diffraction order.
4. Reports the fraction of light transmitted straight through
   (the "zeroth order" transmittance, `T0`).

`main()` then repeats this for an increasing number of harmonics
(1 through 18) and prints how `T0` settles down to a stable value —
that's the convergence study.

**Note on scope:** this is a simplified, educational single-layer,
normal-incidence model built for demonstrating the convergence
methodology, not a full production RCWA solver. A complete RCWA
implementation additionally does an eigenmode decomposition inside each
layer and stitches layers together with a scattering-matrix (S-matrix)
algorithm for numerical stability. Those pieces are natural next steps
if you want to extend this into a general-purpose simulator (see
[Ideas for extending this](#ideas-for-extending-this) below).

## Sample output

```
=========================================================
      RCWA Convergence Study: Modes vs Transmittance
=========================================================
 Harmonics (m) | Total Modes (2m+1) | Transmittance (T0)
---------------------------------------------------------
            1 |                3 |         0.55748243
            2 |                5 |         0.88558651
            3 |                7 |         0.81592991
            ...
           17 |               35 |         0.80343400 <-- Converged
           18 |               37 |         0.80343368 <-- Converged
=========================================================
```

Each row uses more harmonics (and therefore a bigger matrix). Notice how
the transmittance jumps around at first, then settles down — once the
change between rows drops below `1e-5`, the program marks it
`<-- Converged`.

## Building and running

With `make`:

```bash
make        # builds ./rcwa_convergence
make run    # builds (if needed) and runs it
make clean  # removes the built binary
```

Or directly with a compiler:

```bash
g++ -std=c++17 -O2 -o rcwa_convergence src/rcwa_convergence.cpp
./rcwa_convergence
```

## Changing the simulated grating

The physical setup is defined in `main()`:

| Variable        | Meaning                              | Default |
|------------------|---------------------------------------|---------|
| `wavelength`     | Free-space wavelength of the light    | 1.0 µm  |
| `period`         | Repeat distance of the grating        | 1.5 µm  |
| `thickness`      | Grating layer depth                   | 0.5 µm  |
| `n_super`        | Refractive index above the grating    | 1.0 (air) |
| `n_sub`          | Refractive index below the grating    | 1.45 (glass) |
| `n_grating_hi`   | High-index stripe material            | 2.0 |
| `n_grating_lo`   | Low-index stripe material             | 1.0 (air) |
| `duty_cycle`     | Fraction of each period that is the high-index stripe | 0.5 |

Change any of these and re-run to see how the converged transmittance,
and the number of harmonics needed to reach it, change.

## Ideas for extending this

- Add TM polarization (this version approximates TE only).
- Add multiple stacked layers with a proper S-matrix.
- Compute the full set of diffraction orders, not just the 0th.
- Replace the hand-written Gauss-Jordan solver with a proper linear
  algebra library (e.g. Eigen) for speed and numerical robustness at
  large harmonic counts.

## License

MIT — see [`LICENSE`](LICENSE).
