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

## Performance notes

* `dileptons::get_Triangle` / `spectralFunction::Triangle` is a linear scan over the 308 triangles, called **once per cell**. It now tries the previous cell's triangle first
  (neighbouring cells have similar `(n, T)`), but this is a minor effect: the time is dominated by the 4.9 million emission evaluations per cell.
* Originally every evaluation built several `std::vector`s, did linear scans for the `M` and `pT` bins and repeated the planar `(n, T)` interpolation four times, and
  `momenta_integration` took the `dileptons` object **by value**, i.e. copied the whole spectral-function table (about 6 MB) for every cell.
* Now: `dileptons::set_cell(n, T)` finds the triangle and combines the three tables with the barycentric weights into one `89 x 40` array
  (planar interpolation is linear in the table values, so this is exact); `emission_cell_theta` evaluates all 36 azimuth angles of one `(M, pT, y)` with the
  `cosh/sinh/sqrt` terms hoisted and a bilinear `(M, Q)` interpolation using `std::upper_bound` for the bins. The v2 weight `(px^2 - py^2)/(px^2 + py^2)` equals `cos(2(theta + phi0))`, so it is computed once per angle.
  The old per-point functions (`emmision`, `spectral_function`) are kept and used by `make test` as the reference.
* Result: relative difference to the original rates 4e-14 (random cells, including cells outside the triangulated region), 0.015 s instead of about 0.4 s per cell.
  Further speed-up is possible by caching `exp(E*/T)`, which is the largest remaining cost, or by a tabulated `Im D` on a `(M, Q)` grid shared by all cells with the same triangle.

## Square-grid interpolation (default)

The triangle search is replaced by a regular grid in `(T, n)`: `spectralFunction::build_grid(dT = 5, dn = 0.25)` samples the collapsed `M x Q` tables at
`29 x 84` nodes once (shared by all OpenMP threads, 69 MB), filling nodes outside the table region by planar extrapolation from the nearest triangle, and `grid_cell(n, T)`
does a bilinear interpolation of the four surrounding nodes. Whether a point is inside the table region is an O(1) test against the straight boundary lines between the row ends
(outside: the same vacuum value as with `Triangle() == -1`).

Accuracy against the exact triangulation (`make test`, random cells): the collapsed table differs by at most 3.9e-3 of the table maximum (median 2.4e-4), and the emission rate by
at most 1e-2 (mean 3e-4) in relative terms. The data exist only at 180 points, so a finer grid cannot add information; it only approaches the piecewise-linear result
(`dT = 2.5, dn = 0.1`: 338 MB, maximum 1.6e-3). If the spectral functions are ever recomputed on a regular `(T, n)` grid, this bilinear lookup becomes exact.
`dileptons::init(false)` uses the exact triangulation instead (about the same speed; the triangle search is not the bottleneck).
`isin` has a `1e-12` tolerance so that points exactly on a table row (e.g. `T = 120`) are not mistaken for points outside the region.

## Triangulation (`tools/triangulation.py`, `nTRapp.ipynb`)

`nTRapp.ipynb` reads Rapp's output `rho-spectra-function.txt`, extracts the `(T, n/n0)` points of the table blocks (`Tn.dat`), loads `nt.dat` (a `(T, n)` point list, filtered to `T > 50`,
`n > 0.9`) and builds a Delaunay triangulation of it with `matplotlib.tri`. Its point list (`nt.dat` has 93 points) is not the 180-point set of `TriangFiles/nTs.dat`, and the notebook does not
write `triangles.dat`, so it is not clear that it produced the shipped file.

`triangles.dat` has a strip structure: 15 temperature rows with 12 densities each; between neighbouring rows every one of the 11 quadrilaterals is cut into two
triangles along one diagonal (14 strips x 11 x 2 = 308). The choice of diagonals does not follow a
simple rule (shorter diagonal in any scaling, or Delaunay, reproduces it), so keep the shipped file. `python3 tools/triangulation.py` verifies it; `--write <file>` writes a new strip triangulation with shorter diagonals.

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
