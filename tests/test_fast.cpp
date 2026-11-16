// Checks that the per-cell fast path (set_cell / emission_cell) reproduces the original dileptons::emmision
// and compares the run time. Run from the repository root (needs TriangFiles/):  make test
#include "../dileptonsclass.h"
#include <chrono>
#include <cstdio>
#include <random>
#include <algorithm>

int main()
{
	dileptons d; d.init(false);                       // exact triangulation path
	dileptons dg; dg.init(true);                      // square-grid path
	std::mt19937_64 gen(1);
	std::uniform_real_distribution<double> U(0, 1);
	double maxrel = 0, tol_outside = 0; long npts = 0, nin = 0, nout = 0;
	for (int cell = 0; cell < 60; cell++)
	{
		double T = 50 + 135 * U(gen), n = (cell % 6 == 0 ? 3.0 * U(gen) : 0.35 + 5 * U(gen));   // some cells outside the table range
		double vT = 0.6 * U(gen), vz = 0.3 * (U(gen) - .5), gamma = 1 / sqrt(1 - vT * vT - vz * vz);   // physical flow, so that E* >= M
		int k = d.get_Triangle(n, T);
		(k >= 0 ? nin : nout)++;
		d.set_cell(n, T);
		if (d.cell_k != k) { printf("triangle mismatch\n"); return 1; }
		for (int p = 0; p < 4000; p++)
		{
			double M = 20 + 1480 * U(gen), pT = 20 + 2500 * U(gen), y = -.5 + U(gen), th = 2 * PI * U(gen) - PI;
			std::vector<double> q = {M, pT, y}, v = {gamma, vT, vz};
			double a = d.emmision(q, v, k, n, T, th), b = d.emission_cell(M, pT, y, gamma, vT, vz, th);
			double rel = fabs(a - b) / (fabs(a) + 1e-300);
			if (a != a || b != b) { printf("NaN: %g %g\n", a, b); return 1; }
			maxrel = std::max(maxrel, rel); npts++;
			double ct = cos(th), r1; d.emission_cell_theta(M, pT, y, gamma, vT, vz, &ct, 1, &r1);
			maxrel = std::max(maxrel, fabs(a - r1) / (fabs(a) + 1e-300));
		}
	}
	// square-grid path vs exact triangulation: whole collapsed tables and emission rates
	double tab_err = 0, em_err = 0, em_mean = 0; long ncmp = 0, nout_mismatch = 0;
	for (int cell = 0; cell < 400; cell++)
	{
		double T = 40 + 140 * U(gen), n = 0.5 + 21.5 * U(gen);
		d.set_cell(n, T); dg.set_cell(n, T);
		if ((d.cell_k < 0) != (dg.cell_k < 0)) { nout_mismatch++; continue; }
		if (d.cell_k < 0) continue;
		double mx = 0, e = 0;
		for (double x : d.cell) mx = std::max(mx, fabs(x));
		for (size_t q = 0; q < d.cell.size(); q++) e = std::max(e, fabs(d.cell[q] - dg.cell[q]) / mx);
		tab_err = std::max(tab_err, e); ncmp++;
		double vT = 0.5 * U(gen), vz = 0.2 * (U(gen) - .5), gamma = 1 / sqrt(1 - vT * vT - vz * vz);
		for (int p = 0; p < 200; p++)
		{
			double M = 20 + 1480 * U(gen), pT = 20 + 2500 * U(gen), y = -.5 + U(gen), th = 2 * PI * U(gen) - PI;
			double a = d.emission_cell(M, pT, y, gamma, vT, vz, th), b = dg.emission_cell(M, pT, y, gamma, vT, vz, th);
			double r = fabs(a - b) / (fabs(a) + 1e-300); em_err = std::max(em_err, r); em_mean += r;
		}
	}
	printf("grid vs triangulation: %ld cells inside, %ld inside/outside mismatches (points on the boundary), max table error %.2e (relative to table max)\n", ncmp, nout_mismatch, tab_err);
	printf("   emission rate: max relative difference %.2e, mean %.2e\n", em_err, em_mean / (200. * ncmp));
	printf("cells inside the table region: %ld, outside: %ld, points: %ld\n", nin, nout, npts);
	printf("max relative difference old/new: %.2e\n", maxrel);
	// timing for one cell (1 cell = 88*39*40*36 = 4.9M points)
	double n = 1.0, T = 120, gamma = 1.2, vT = .3, vz = .1; int k = d.get_Triangle(n, T);
	auto t0 = std::chrono::steady_clock::now(); double s1 = 0;
	for (int i = 1; i < 89; i += 8) for (int j = 1; j < 40; j += 4) for (int l = 0; l < 36; l += 4)
	{ std::vector<double> q = {100. + 15 * i, 50. + 100 * j, 0.1}, v = {gamma, vT, vz}; s1 += d.emmision(q, v, k, n, T, (l + .5) * PI / 18); }
	double dt_old = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	long cnt_old = 11L * 10 * 9;
	t0 = std::chrono::steady_clock::now(); double s2 = 0; d.set_cell(n, T);
	for (int rep = 0; rep < 200; rep++)
	for (int i = 1; i < 89; i += 8) for (int j = 1; j < 40; j += 4) for (int l = 0; l < 36; l += 4)
		s2 += d.emission_cell(100. + 15 * i + 1e-7 * rep, 50. + 100 * j, 0.1, gamma, vT, vz, (l + .5) * PI / 18);
	double dt_new = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	double ct[36], rr[36]; for (int l = 0; l < 36; l++) ct[l] = cos((l + .5) * PI / 18);
	t0 = std::chrono::steady_clock::now(); double s3 = 0;
	for (int rep = 0; rep < 200; rep++)
	for (int i = 1; i < 89; i += 8) for (int j = 1; j < 40; j += 4) { d.emission_cell_theta(100. + 15 * i + 1e-7 * rep, 50. + 100 * j, 0.1, gamma, vT, vz, ct, 36, rr); for (int l = 0; l < 36; l += 4) s3 += rr[l]; }
	double dt_row = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	printf("time per point, theta-row version: %.3f us (checksum %.3g)\n", 1e6 * dt_row / (200. * 11 * 10 * 36), s3 / 200);
	printf("time per point: old %.3f us, new %.3f us, speed-up %.0fx  (checksums %.3g %.3g)\n",
		1e6 * dt_old / cnt_old, 1e6 * dt_new / (200. * cnt_old), (dt_old / cnt_old) / (dt_new / (200. * cnt_old)), s1, s2 / 200);
	// one complete cell with the row version: 88 M bins x 39 pT bins x 40 rapidities x 36 angles
	t0 = std::chrono::steady_clock::now(); double s4 = 0;
	for (int i = 1; i < 89; i++) for (int j = 1; j < 40; j++) for (int k = 0; k < 40; k++)
	{ d.emission_cell_theta(20. * i, 50. + 100 * j, -.5 + (k + .5) / 40, gamma, vT, vz, ct, 36, rr); s4 += rr[0]; }
	printf("one full cell (4.9M points), row version: %.3f s  (checksum %.3g)\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), s4);
	return maxrel < 1e-10 ? 0 : 1;
}
