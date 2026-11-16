# v2dileptons - documentation

## Physics and formulas

For a fluid cell with temperature `T` [MeV], density `n`, flow (`gamma`, `v_T`, `v_z`, direction `phi0`) and pion factor, the dilepton
emission rate at invariant mass `M`, transverse momentum `pT`, rapidity `y` and azimuth `theta` (relative to the cell's flow direction) is

```
mT   = sqrt(M^2 + pT^2),   p0 = mT cosh y
E*   = gamma (p0 - v_T pT cos(theta) - v_z mT sinh y)               (energy in the cell frame, dileptons::qCM)
q*   = sqrt(E*^2 - M^2)
dN/(d4x dM dpT ...) = -(2 pi M pT) alpha_EM^2 / pi^3 * L(M)/M^2 * f_BE(E*, T) * Im D(M, q*; n, T)        (dileptons::emmision, L = 1)
```

`Im D` is the in-medium vector-meson spectral function, interpolated from the tables (below) and multiplied by a coupling-type
normalisation `g = 829^4/(2.7*4 pi)` (`spectralFunction::g`); the sum of the `rho` column only is used (`cells[j][i] = line[2]`, the omega and phi columns are commented out).
Outside the triangulated `(n, T)` region (`Triangle() == -1`) the code returns the analytic continuum value `-3 M^2/(12 pi) * 5/9` instead.

The elliptic-flow weight is `v2 = (px^2 - py^2)/(px^2 + py^2)` with `px = pT cos(phi + phi0)`, `py = pT sin(phi + phi0)`.
An extra empirical factor `lambda_pion(n_pi, T) = (n_pi / (3 n_T0(T)))^1.3` (a quartic polynomial `n_T0(T)` in `main.cpp`) accounts for the pion chemical potential.

## Program flow (`main.cpp`)

1. Read `T, n, gamma, vx, vy, vz, pion` (binary `double`, `100^4` each) from `<data_dir>`.
2. One `dileptons` object per OpenMP thread; each loads the spectral-function tables in `dileptons::init()` (`./TriangFiles/`).
3. For each time step `t` in `[t_start, t_end]`: loop over the central `60 x 60 x 60` cells (indices 20..79 of 100, i.e. +-15 fm at `dx = 0.5 fm`),
   parallel over the `y` index. Cells enter only if `T > 50 MeV` and `n * 0.13868509 > 0.144` (with `Nemp = n * 0.13868509/0.16` passed on as the density).
4. For every selected cell `momenta_integration` sums the rate over `M` bins 1..88 (89 bin edges), `pT` bins 1..39 (40 edges), 40 rapidity bins in `[-0.5, 0.5]`
   and 36 azimuth bins, accumulating `yield` and `v2 * yield` per thread.
5. `output_PT` adds the threads and writes, per `(M, pT)`, the yield and the `v2`-weighted yield times `A * (2 y_max/40) * dx^4` with `A = 1000/197^4` (fm^-4 to GeV^-4 type conversion, factor 1000 for MeV).

Output: `<data_dir>/rho-V2-<t>.bin` (text, tab-separated: `M  pT  yield  v2*yield`). The v2 is then the ratio of the integrated columns.

## Spectral-function tables (`TriangFiles/`)

| File | Content |
|---|---|
| `ImDvec_T_<T>.out` | 29 files (`T = 40 ... 180 MeV` in steps of 5, 12 MB each at most). After a 197-line header, one block per baryon density (57 blocks) with `Im D_rho, Im D_omega, Im D_phi` on the `M x pT` grid (89 x 40) |
| `nTs.dat` | the 180 `(T, n)` points actually used (T = 40, 50, ..., 180, 12 densities each) |
| `triangles.dat` | 308 triangles (indices into `nTs.dat`) covering the `(n, T)` plane |
| `nt.dat` | 106 `(T, n)` pairs (not read by the code) |
| `[0-9]*-[0-9]*.dat` | 4274 files (`<a>-<b>.dat` with `a` = 0..104, `b` = 0..40, 90 lines each), not read by the code |

`spectralFunction::set_table` reads every second file (the array `Ts[15] = 40, 50, ..., 180`) and every fifth density block (`frequency_Ns = 5`, 57 blocks -> 12), which gives the 180 points of `nTs.dat`.
Interpolation: planar interpolation in `(n, T)` over the triangle that contains the point, then bilinear in `(M, pT)` between the surrounding bins.

## Cluster scripts

`dileptons.sh` is a SLURM batch script (20 tasks, 1 GB each, 40 h); it runs `./a.out $SCRATCHDIR/Dileptons <t_start> <t_end> 20` with `t_start`, `t_end` taken from its 2nd and 3rd argument.
`dileptons-script.bash` prints the `sbatch` commands for the equations of state `v2cmf, ptv2, pt2v2` and beam energies `123, 200, 400, 600, 1000` (the path it prints as first argument is not used by the job script; the data folder comes from `SCRATCHDIR`, which is a cluster path you need to change).

## Things to know (left unchanged)

* The program has to be started in the folder with `TriangFiles/`; the number of triangles/points (`308`, `180`) is hard coded in `dileptons::init`.
* `main` returns 1 on success; four unused variables give compiler warnings (`e0, n0, d4x, xmax` in `main`).
* The Mbin/Qbin tables are repeated in `main.cpp` (twice) and in `spectral_func.h`; change all copies together.
* Memory: all seven `100^4` arrays plus the per-thread yield arrays (about 5.6 GB in total).
* `dileptons.sh` ignores its first argument (see above); the job name and `--partition` are cluster specific.
* Rapp-type `ImDvec` headers mention `mu_B` and `rhoN_eff`; the code identifies the table density only via the order of the blocks (`nTs.dat`), not by reading these header lines.
