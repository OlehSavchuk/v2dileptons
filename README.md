# v2dileptons

C++/OpenMP code that computes the **thermal dilepton yield and its elliptic flow v2** from a hydrodynamic /
transport evolution of a heavy-ion collision. For every space-time cell with temperature `T`, baryon density `n`,
flow velocity and pion chemical potential it evaluates the in-medium vector-meson (rho, omega, phi) spectral function,
folds it with a Bose distribution and the flow boost, integrates over the cell sizes and writes the yield and the
`v2`-weighted yield as a function of invariant mass `M` and transverse momentum `pT`.

The spectral functions are tabulated on a `(T, n)` grid (`TriangFiles/ImDvec_T_*.out`) and interpolated linearly on a
triangulation of the grid in the `(n, T)` plane (`triangles.dat`), then bilinearly in `(M, pT)`.

## Contents

| File | Purpose |
|---|---|
| `dilepton_yield.cpp` | Driver: reads the cell data, loops over cells and time steps (OpenMP), integrates over `M, pT, y, phi`, writes `rho-V2-<t>.bin` |
| `dileptons.h/.cpp` | `dileptons` class: emission rate per phase-space point (`emmision`), energy in the cell frame (`qCM`), Bose factor |
| `spectral_function.h/.cpp` | `spectralFunction` class: reads the `ImDvec_T_*.out` tables, finds the triangle containing `(n, T)`, interpolates |
| `TriangFiles/` | Spectral-function tables and the `(n, T)` triangulation (see below) |
| `run_dileptons.slurm` | SLURM job script (20 tasks, 40 h) |
| `submit_all.sh` | Prints the `sbatch` lines for 3 EoS variants times 5 beam energies |
| `Makefile` | `make` builds `dilepton_v2`; `make test` checks the fast path against the original rate (needs the `ImDvec` tables); `make check-triangulation`; `make clean` |
| `push.sh` | Creates the private GitHub repo and pushes (needs `gh auth login`) |
| `references.bib` | BibTeX for the papers and the spectral-function references |
| `tests/test_fast.cpp` | Compares `emmision` with the per-cell fast path and times both |
| `rapp_nT_points.ipynb`, `rho-spectra-function.txt`, `Tn.dat`, `nt.dat` | Notebook and Rapp spectral-function output it reads (extracts the `(T, n)` points, Delaunay triangulation) |
| `tools/triangulation.py` | Verifies `triangles.dat` (strip structure, coverage) and can write a new strip triangulation |
| `docs/index.md` | Detailed documentation |

## Build and run

```bash
make                                   # g++ -O2 -std=c++14 -fopenmp ... -> ./dilepton_v2
./dilepton_v2 <data_dir> <t_start> <t_end> <threads>
```

* `<data_dir>` holds the cell data as raw binary `double` arrays, 100 x 100 x 100 x 100 each (time, x, y, z; `dx = 0.5 fm`):
  `T.dat`, `n.dat`, `gamma.dat`, `vx.dat`, `vy.dat`, `vz.dat`, `pion.dat`. Memory use is about 5.6 GB (all seven arrays are held at once).
* the program must be started from the folder that contains `TriangFiles/` (path `./TriangFiles/` is hard coded in `dileptons::init`).
* output: `<data_dir>/rho-V2-<t>.bin` for every time step `t` (a text file despite the extension): columns `M`, `pT`, yield, `v2`-weighted yield.

## Speed

The triangle search itself (once per cell, 308 triangles) was never the bottleneck: the cost is the 4.9 million emission evaluations per cell
(88 `M` x 39 `pT` x 40 `y` x 36 azimuth bins). The code now does the triangle search and the planar `(n, T)` interpolation of the whole
`M x pT` table once per cell, hoists everything that does not depend on the azimuth out of the inner loop, and uses binary search for the bins.
`make test`: with the exact triangulation the results agree with the original `dileptons::emmision` to 4e-14 (relative), and one cell takes about 0.015 s instead of about 0.4 s (25x).


## Data

`TriangFiles/` is 355 MB. Git tracks only the small files (`nTs.dat`, `nt.dat`, `triangles.dat`). Not tracked (see `.gitignore`):
the 29 spectral-function tables `ImDvec_T_<T>.out` (338 MB, needed to run) and 4274 files `N-N.dat` (17 MB, not read by this code).
Every table is below GitHub's 100 MB limit (largest 12 MB), so if you want them in the repository, delete the two lines from `.gitignore`.

## References

**Please cite the papers produced with this code:**

* O. Savchuk, A. Motornenko, J. Steinheimer, V. Vovchenko, M. Bleicher, M. Gorenstein and T. Galatyuk, *Enhanced dilepton emission from a phase transition in dense matter*,
  J. Phys. G **50**, 125104 (2023), [arXiv:2209.05267](https://arxiv.org/abs/2209.05267).
* T. Reichert, O. Savchuk, A. Kittiratpattana, P. Li, J. Steinheimer, M. Gorenstein and M. Bleicher, *Decoding the flow evolution in Au+Au reactions at 1.23A GeV using hadron flow correlations and dileptons*,
  Phys. Lett. B (2023), [arXiv:2302.13919](https://arxiv.org/abs/2302.13919).

**and the spectral functions** (the tables in `TriangFiles/ImDvec_T_*.out` and `rho-spectra-function.txt` are the in-medium rho, omega and phi spectral functions of Rapp, Wambach and van Hees):

* R. Rapp and J. Wambach, Eur. Phys. J. A **6**, 415 (1999), [arXiv:hep-ph/9907502](https://arxiv.org/abs/hep-ph/9907502).
* H. van Hees and R. Rapp, Nucl. Phys. A **806**, 339 (2008), [arXiv:0711.3444](https://arxiv.org/abs/0711.3444).
* R. Rapp, Adv. High Energy Phys. **2013**, 148253 (2013), [arXiv:1304.2309](https://arxiv.org/abs/1304.2309).
* Review: R. Rapp and J. Wambach, Adv. Nucl. Phys. **25**, 1 (2000), [arXiv:hep-ph/9909229](https://arxiv.org/abs/hep-ph/9909229).

BibTeX entries: [references.bib](references.bib). The tables were produced by R. Rapp's model; they are used here with his results as published in the papers above, so please also acknowledge the source of the tables when you redistribute them.

Further documentation: [docs/index.md](docs/index.md)
