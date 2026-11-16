
#ifndef SPECTRAL_FUNC_H
#define SPECTRAL_FUNC_H
#define PI 3.14159265
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
class spectralFunction
{
    double g = 829. * 829. / 2.7 / 4 / PI * 829 * 829;
    double e_sq = 5. / 9;
    int Mbin(double);
    int Qbin(double);
    double Mbins[89] = { 1.000e+00, 5.000e+01, 1.000e+02, 1.500e+02, 2.000e+02, 2.500e+02,
        3.000e+02, 3.500e+02, 4.000e+02, 4.250e+02, 4.500e+02, 4.750e+02,
        5.000e+02, 5.250e+02, 5.500e+02, 5.750e+02, 6.000e+02, 6.200e+02,
        6.400e+02, 6.600e+02, 6.800e+02, 7.000e+02, 7.100e+02, 7.200e+02,
        7.300e+02, 7.400e+02, 7.450e+02, 7.500e+02, 7.550e+02, 7.600e+02,
        7.620e+02, 7.640e+02, 7.660e+02, 7.680e+02, 7.700e+02, 7.720e+02,
        7.740e+02, 7.760e+02, 7.780e+02, 7.800e+02, 7.810e+02, 7.820e+02,
        7.830e+02, 7.840e+02, 7.860e+02, 7.880e+02, 7.900e+02, 7.920e+02,
        7.950e+02, 8.000e+02, 8.050e+02, 8.100e+02, 8.200e+02, 8.300e+02,
        8.400e+02, 8.500e+02, 8.600e+02, 8.700e+02, 8.800e+02, 9.000e+02,
        9.200e+02, 9.400e+02, 9.600e+02, 9.800e+02, 9.900e+02, 1.000e+03,
        1.005e+03, 1.010e+03, 1.015e+03, 1.017e+03, 1.018e+03, 1.019e+03,
        1.020e+03, 1.021e+03, 1.022e+03, 1.025e+03, 1.030e+03, 1.035e+03,
        1.040e+03, 1.050e+03, 1.060e+03, 1.075e+03, 1.100e+03, 1.150e+03,
        1.200e+03, 1.250e+03, 1.300e+03, 1.400e+03, 1.500e+03 };
    double Qbins[40] = { 1.00e+00, 5.00e+01, 1.00e+02, 1.50e+02, 2.00e+02, 2.50e+02,
        3.00e+02, 3.50e+02, 4.00e+02, 4.50e+02, 5.00e+02, 5.50e+02,
        6.00e+02, 6.50e+02, 7.00e+02, 7.50e+02, 8.00e+02, 8.50e+02,
        9.00e+02, 9.50e+02, 1.00e+03, 1.05e+03, 1.10e+03, 1.20e+03,
        1.35e+03, 1.50e+03, 1.75e+03, 2.00e+03, 2.25e+03, 2.50e+03,
        2.75e+03, 3.00e+03, 3.25e+03, 3.50e+03, 3.75e+03, 4.00e+03,
        4.25e+03, 4.50e+03, 4.75e+03, 5.00e+03 };
    int num_Mbins = 89;
    int num_Qbins = 40;
    int num_Ts = 15;
    int num_Ns = 57;
    int frequency_Ns = 5;
    double Ts[15] = { 40,  50,  60,  70,  80,  90, 100, 110, 120, 130, 140, 150, 160,
       170, 180 };
    std::vector<std::vector<double>> nT;
    std::vector<std::vector<int>> triangles;
    std::vector<std::vector<std::vector<double>>> table;
    double PlanarInterpolate(double, double, std::vector<double>&, std::vector<double>&, std::vector<double>&);
    bool isin(int*, double n, double T);
    double getSF(int, int, int, double, double);
    double LinearInterpolate(double, double, double, double, std::vector<std::vector<double>>&);
    void readEvent(std::vector<std::vector<double>>&, std::ifstream&);
    void readHeader(std::ifstream& );
public:
    void set_table(std::string, int,int);
    double spectral_function(double, double, int , double , double );
    int Triangle(double, double);
};
#endif
