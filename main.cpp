#include "dileptonsclass.h"
#include <algorithm>
#include <string>
#include <omp.h>
#include <math.h>


double lambda_pion(double n, double T)
{
	double npow=1.3;
	double a0=9.63046*1e-5;
	double a1=6.11838*1e-6;
	double a2=-7.8558*1e-7;
	double a3=1.70236*1e-8;
	double a4=-6.71064*1e-12;
	double nT0 = a0+a1*T+a2*pow(T,2)+a3*pow(T,3)+a4*pow(T,4);
	return pow(n/nT0/3,npow);
}
double v2(std::vector<double>& q, double phi0,double phi)
{
		double theta = phi + phi0;
		double px = q[1] * cos(theta);
		double py = q[1] * sin(theta);
		px = px * px;
		py = py * py;
		return (px - py) / (px + py);
}

void write_to_file(std::vector<double>& a, std::ofstream& f)
{
	f.write(reinterpret_cast<char*>(&a[0]),sizeof(double)*a.size());
}
void read_from_file(std::vector<double>& b, std::ifstream& f,int size)
{
	std::vector<double> read(size,0);
	f.read(reinterpret_cast<char*>(&read[0]),sizeof(double)*size);
	for(int i = 0;i<size;i++)
	{
		b[i]+=read[i];
	}
}

void read_cells_element(std::string element,std::vector<std::vector<std::vector<std::vector<double>>>>& cells, int size)
{
	std::ifstream file(element,std::ios::binary);
	for(int t=0;t<size;t++)
	{
		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				read_from_file(cells[t][x][y],file,size);
			}
		}
	}
}

void momenta_integration(std::vector<std::vector<std::vector<std::vector<double>>>>& yields, double Nemp, double Temp, std::vector<double> v, double phi0, double lpion, dileptons& llbar)
{
	double ymax = .5;
	int beansM = 89;
	int beansQ = 40;
	double e0 = 0.127948296;
	double n0 = 0.13868509;
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
	

	// Fast path: triangle search + planar (n,T) interpolation of the spectral function once per cell,
	// then only a bilinear (M,Q) interpolation per phase-space point (see dileptons::set_cell / emission_cell).
	llbar.set_cell(Nemp, Temp);
	double cos2[36], cosT[36], rate[36];               // v2 weight = cos(2(theta + phi0)), independent of pT
	for (int l = -18; l < 18; l++)
	{
		cos2[l + 18] = cos(2 * ((l + .5) * PI / 18 + phi0));
		cosT[l + 18] = cos((l + .5) * PI / 18);
	}
	for (int i = 1; i < beansM; i++)
	{
		double M = Mbins[i];
		for (int j = 1; j < beansQ; j++)
		{
			double pt = Qbins[j];
			for (int k = 0; k < beansQ; k++)
			{
				double y = -ymax + 2 * ymax / beansQ * (k + .5);
				llbar.emission_cell_theta(M, pt, y, v[0], v[1], v[2], cosT, 36, rate);
				double s0 = 0, s1 = 0;
				for (int l = 0; l < 36; l++)
				{
					double tR = 1. / 36 * lpion * rate[l];
					s0 += tR;
					s1 += cos2[l] * tR;
				}
				yields[0][i][j][k] += s0;
				yields[1][i][j][k] += s1;
			}
		}
	}
}
void output_PT(std::string path, int t, std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& yields, int num_threads)
{
	int beansM = 89;
	int beansQ = 40;
	double ymax = .5;
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
	double dx = .5;
	double xmax = 50.;
	double dV = dx * dx * dx;
	double d4x = dV * dx;
	double A = 1000. / (pow(197, 4));
	std::ofstream V2(path + "/rho-V2-" + std::to_string(t) + ".bin");
	for (int i = 0; i < beansM; i++)
	{
		double M = Mbins[i];


		for (int j = 0; j < beansQ; j++)
		{
			V2 << M << "	" << Qbins[j] << "	";
			double yield1 = 0;
			double yield2 = 0;
			for (int k = 0; k < beansQ; k++)
			{

				for (int l = 0; l < num_threads; l++)
				{
					yield1 += yields[l][0][i][j][k];
					yield2 += yields[l][1][i][j][k];
				}
			}
			V2 << A * yield1 * 2 * ymax / beansQ * d4x << "	" << A * yield2 * 2 * ymax / beansQ * d4x <<std::endl;
		}
		
	}
	V2.close();
}
/*
void output_M(std::string path,int t, std::vector<std::vector<std::vector<std::vector<double>>>> &yields, int num_threads)
{
	int beansM = 89;
	int beansQ = 40;
	double ymax = .5;
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
	double dx = .5;
	double xmax = 50.;
	double dV = dx * dx * dx;
	double d4x = dV * dx;
	double A = 1000. / (pow(197, 4));
	std::ofstream massSpectra(path + "/complete-massSpectra-" + std::to_string(t) + ".dat");
	for (int i = 1; i < beansM; i++)
	{
		double M = Mbins[i];
		double yield=0;
		massSpectra << M << "	";
		for (int j = 1; j < beansQ; j++)
		{
			for (int k = 0; k < beansQ; k++)
			{
				for (int l = 0; l < num_threads; l++)
				{
					yield += .5 * ( yields[l][i][j][k] + yields[l][i][j-1][k] ) * ( Qbins[j] - Qbins[j - 1] );
				}
			}
		}
		//std::cout << yield << std::endl;
		massSpectra << A * yield * 2 * ymax / beansQ * d4x  << std::endl;

	}
	massSpectra.close();
}*/
int main(int argc, char** argv)
{
	double dx = .5;
	double xmax=50.;
	double dV = dx*dx*dx;
	int nbins=xmax/dx;
	double d4x=dV*dx;
	int num_threads = atoi(argv[4]);
	omp_set_num_threads(num_threads);
	double e0 = 0.127948296;
	double n0 = 0.13868509;
	enum columns{ T,n,gamma,vx,vy,vz,pi };
	int beansM = 89;
	int beansQ = 40;
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> Tnv(7,std::vector<std::vector<std::vector<std::vector<double>>>>(nbins,std::vector<std::vector<std::vector<double>>>(nbins,std::vector<std::vector<double>>(nbins,std::vector<double>(nbins,0)))));    

	std::string path = argv[1];
	int tstart = atoi(argv[2]);
	int tend = atoi(argv[3]);
	std::cout<<"Begin reading temperatures densities etc"<<std::endl;
	read_cells_element(path+"/T.dat",Tnv[T],nbins);
	read_cells_element(path+"/n.dat",Tnv[n],nbins);
	read_cells_element(path+"/gamma.dat",Tnv[gamma],nbins);
	read_cells_element(path+"/vx.dat",Tnv[vx],nbins);
	read_cells_element(path+"/vy.dat",Tnv[vy],nbins);
	read_cells_element(path+"/vz.dat",Tnv[vz],nbins);
	read_cells_element(path+"/pion.dat",Tnv[pi],nbins);
	std::cout<<"Finished reading temperatures densities etc"<<std::endl;

	
	dileptons llbars[num_threads];
	for (int i = 0; i < num_threads; i++)
	{
		llbars[i].init();
	}
	
	
	for (int t = tstart; t <= tend; t++)
	{
		std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> yields(num_threads, std::vector<std::vector<std::vector<std::vector<double>>>>(2, std::vector<std::vector<std::vector<double>>>(beansM, std::vector<std::vector<double>>(beansQ, std::vector<double>(beansQ, 0)))));
		for (int i = 20; i < 80; i++)
		{
#pragma omp parallel
			{
				int thread_num = omp_get_thread_num();
				for (int tj = 0; tj < (80-20) / num_threads; tj++)
				{
					
					int j = 20 + thread_num + tj * num_threads;
					for (int k = 20; k < 80; k++)
					{
						if (Tnv[T][t][i][j][k] > 50 && Tnv[n][t][i][j][k] * n0 > .144)
						{
							std::vector<double> v(3,0);
							v[0] = Tnv[gamma][t][i][j][k];
							v[1] = sqrt(pow(Tnv[vx][t][i][j][k], 2) + pow(Tnv[vy][t][i][j][k], 2));
							v[2] = Tnv[vz][t][i][j][k];
							double phi0 = atan2(Tnv[vy][t][i][j][k], Tnv[vx][t][i][j][k]);
							double Temp = Tnv[T][t][i][j][k];
							double Nemp = Tnv[n][t][i][j][k] * n0 / .16;
							double lpion = lambda_pion(e0 * Tnv[pi][t][i][j][k], Temp);
							momenta_integration(yields[thread_num], Nemp, Temp, v,phi0, lpion, llbars[thread_num]);
						}
					}
				}
			}
		}
		output_PT(path, t, yields,num_threads);
	}
	return 1;
}
