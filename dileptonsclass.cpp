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

