#include "dileptonsclass.h"
int dileptons::get_Triangle(double n, double T)
{
	return Rapp.Triangle(n, T);
}
double dileptons::spectralFunctionR(double M,double q, int k, double n, double T)
{
	return Rapp.spectral_function(M, q, k, n, T);
}
double dileptons::qCM(std::vector<double>& q, std::vector<double>& v,double theta)
{
	
	mT=sqrt(q[0]*q[0]+q[1]*q[1]);
	p0=mT*cosh(q[2]);
	qcm = v[0]*(p0-v[1]*q[1]*cos(theta)-v[2]*mT*sinh(q[2]));
	//std::cout<<"qcm: "<<qcm<<std::endl;
	return qcm;
	
}
void dileptons::init()
{
	Rapp.set_table("./TriangFiles/",308,180);
}

double dileptons::emmision(std::vector<double>& q, std::vector<double>& v,int k, double n, double T,double theta)
{
	M2=q[0]*q[0];
	double QCM = qCM(q, v, theta);
	double e = -(2*PI*q[0]*q[1])*alphaEM2/PI3*dileptons::L(q[0])/M2*dileptons::fBE(QCM,T)*dileptons::spectralFunctionR(q[0], sqrt(QCM*QCM- q[0] * q[0]),k,n,T);
	
	return e;
}

double dileptons::fBE(double E, double T)
{

	//std::cout<<"first"<<1./(exp(E/T)-1)<<std::endl;
	return 1./(exp(E/T)-1);
}

double dileptons::L(double M)
{
	return 1.;
}


void dileptons::set_cell(double n, double T)
{
	cell_T = T;
	cell_k = Rapp.Triangle(n, T);
	Rapp.cell_table(n, T, cell_k, cell);
}

double dileptons::emission_cell(double M, double pT, double y, double gamma, double vT, double vz, double theta)
{
	// identical to emmision(q, v, k, n, T, theta) for the cell prepared by set_cell(n, T); T enters only through f_BE
	double mTv = sqrt(M * M + pT * pT);
	double QCM = gamma * (mTv * cosh(y) - vT * pT * cos(theta) - vz * mTv * sinh(y));
	return -(2 * PI * M * pT) * alphaEM2 / PI3 * L(M) / (M * M) * fBE(QCM, cell_T) *
		Rapp.spectral_function_cell(cell, cell_k, M, sqrt(QCM * QCM - M * M));
}

void dileptons::emission_cell_theta(double M, double pT, double y, double gamma, double vT, double vz, const double* cosTheta, int n, double* out)
{
	double mTv = sqrt(M * M + pT * pT);
	double A = gamma * (mTv * cosh(y) - vz * mTv * sinh(y)), B = gamma * vT * pT;
	double pref = -(2 * PI * M * pT) * alphaEM2 / PI3 * L(M) / (M * M), M2 = M * M, invT = 1. / cell_T;
	for (int l = 0; l < n; l++)
	{
		double QCM = A - B * cosTheta[l];
		out[l] = pref / (exp(QCM * invT) - 1) * Rapp.spectral_function_cell(cell, cell_k, M, sqrt(QCM * QCM - M2));
	}
}
