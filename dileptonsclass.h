#ifndef DILEPTONSCLASS_H
#define DILEPTONSCLASS_H

#define PI 3.14159265
#define alphaEM 0.0072973525693
#include <fstream>
#include "spectral_func.h"
#include <iostream>
#include <vector>
#include <math.h>
class dileptons{
	spectralFunction Rapp;
	double mT,qcm,p0,M2;
	double PI3=PI*PI*PI;
	double alphaEM2=alphaEM*alphaEM;
public:
	double emmision(std::vector<double>&, std::vector<double>&, int,double, double,double);
	double spectralFunctionR(double, double, int,double,double);
	double qCM(std::vector<double>&,std::vector<double>&,double);
	double fBE(double, double);
	double L(double);
	void init();
	// fast path: prepare a cell once, then evaluate many phase-space points
	std::vector<double> cell;
	int cell_k = -1;
	double cell_T = 1;
	void set_cell(double n, double T);   // also stores T
	// all azimuth angles at once: cosTheta[0..n-1] given, out[0..n-1] = emission rate (hoists the cosh/sinh/sqrt that do not depend on theta)
	void emission_cell_theta(double M, double pT, double y, double gamma, double vT, double vz, const double* cosTheta, int n, double* out);
	double emission_cell(double M, double pT, double y, double gamma, double vT, double vz, double theta);
	int get_Triangle(double, double);

};
#endif
