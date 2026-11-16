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
	int get_Triangle(double, double);

};
#endif
