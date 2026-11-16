#include "spectral_func.h"
#include <algorithm>

int spectralFunction::Mbin(double M)
{
	for(int i=1;i<num_Mbins;i++)
	{
		if(M<Mbins[i])
		{
			return i;
		}
	}
	return num_Mbins-1;
}

int spectralFunction::Qbin(double q)
{
	for(int i=1;i<num_Qbins;i++)
	{
		if(q<Qbins[i])
		{
			return i;
		}
	}
	return num_Qbins-1;
}
double spectralFunction::PlanarInterpolate(double x, double y, std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs)
{
	double r10[3] = { xs[1] - xs[0],ys[1] - ys[0],zs[1] - zs[0] };
	double r20[3] = {xs[2] - xs[0], ys[2] - ys[0], zs[2] - zs[0]};
	double r10xr20[3] = { r10[1] * r20[2] - r10[2] * r20[1], r10[2] * r20[0] - r10[0] * r20[2], r10[0] * r20[1] - r10[1] * r20[0] };
	return -(r10xr20[0] * (x - xs[0]) + r10xr20[1] * (y - ys[0])) / r10xr20[2] + zs[0];
}
double spectralFunction::LinearInterpolate(double x, double y,double dx, double dy,std::vector<std::vector<double>>& points)
{
	
	double a1 = (points[0][1] - points[0][0])*y/dy+points[0][0];
	double a2 = (points[1][1] - points[1][0])*y/dy+points[1][0];
	//std::cout<<"dx "<< dx <<"	"<<"dy "<<dy<<std::endl;

	return (a2 - a1)*x/dx+a1;
}

void spectralFunction::readHeader(std::ifstream& file)
{
	int lines_to_skip = 197;
	for (int i = 0; i < lines_to_skip; i++)
	{
		file.ignore(1024, '\n');
	}
}

void spectralFunction::readEvent(std::vector<std::vector<double>>& cells,std::ifstream& file)
{
	int column_num = 5;
	double line[column_num];
	file.ignore(1024, '\n');
	file.ignore(1024, '\n');
	file.ignore(1024, '\n');
	for (int i = 0; i < num_Qbins; i++)
	{
		file.ignore(1024, '\n');
		//	std::cout<<particle_num<<std::endl;
		//       std::cout<<"timestep: " << i << std::endl; 
		for (int j = 0; j < num_Mbins; j++)
		{
			for (int k = 0; k < column_num; k++)
			{
				file.ignore(1, '\n');
				file >> line[k];
				//std::cout << line[k] << std::endl;
			}
			cells[j][i] = line[2];//+line[3]+line[4];
			//std::cout << cells[j][i] << std::endl;
		}
	}
	file.ignore(1024, '\n');

}

void spectralFunction::set_table(std::string path, int num1, int num2)
{
	double n;
	double T;
	std::vector<double> tnT;
	std::vector<int> triangle;
	std::ifstream fnT(path + "nTs.dat");
	std::ifstream ftriangles(path + "triangles.dat");
	double triangle_num;
	//std::cout << "begin reading triangles" << std::endl;
	for (int i = 0; i < num1; i++)
	{

		//std::cout << "begin reading triangle:  " << i << std::endl;
		triangle.clear();
		ftriangles>>triangle_num;
		triangle.push_back(triangle_num);
		ftriangles >> triangle_num;
		triangle.push_back(triangle_num);
		ftriangles >> triangle_num;
		triangle.push_back(triangle_num);
		triangles.push_back(triangle);
	}
	for (int i = 0; i < num2; i++)
	{

		//std::cout << "begin reading nts:  " << i << std::endl;
		tnT.clear();
		fnT >> n;
		fnT >> T;
		tnT.push_back(n);
		tnT.push_back(T);
		nT.push_back(tnT);
	}

	//std::cout << "begin reading spf" << std::endl;
	for (int i = 0; i < 6; i++)
	{
		std::vector<std::vector<double>> SFMQ(num_Mbins,std::vector<double>(num_Qbins,0));
		std::cout << path + "ImDvec_T_0" + std::to_string(int(Ts[i])) + ".out" << std::endl;
		std::ifstream tpfile(path +"ImDvec_T_0"+ std::to_string(int(Ts[i])) + ".out");
		spectralFunction::readHeader(tpfile);
		for (int j = 0; j < num_Ns; j++)
		{
			spectralFunction::readEvent(SFMQ, tpfile);
			if (j % frequency_Ns == 0)
			{
				table.push_back(SFMQ);
			}
		}
		
	}
	for (int i = 6; i < num_Ts; i++)
	{
		std::vector<std::vector<double>> SFMQ(num_Mbins, std::vector<double>(num_Qbins, 0));
		std::cout << path + "ImDvec_T_" + std::to_string(int(Ts[i])) + ".out" << std::endl;
		std::ifstream tpfile(path + "ImDvec_T_" + std::to_string(int(Ts[i])) + ".out");
		spectralFunction::readHeader(tpfile);
		for (int j = 0; j < num_Ns; j++)
		{
			spectralFunction::readEvent(SFMQ, tpfile);
			if (j % frequency_Ns == 0)
			{
				table.push_back(SFMQ);
			}
		}

	}
	//std::cout << "FINISHED reading spf" << std::endl;
}
double spectralFunction::getSF(int i, int j, int k, double n, double T)
{
	std::vector<double> xs;
	std::vector<double> ys;
	std::vector<double> zs;
	xs.push_back(nT[triangles[k][0]][1]);
	ys.push_back(nT[triangles[k][0]][0]);
	zs.push_back(table[triangles[k][0]][i][j]);
	xs.push_back(nT[triangles[k][1]][1]);
	ys.push_back(nT[triangles[k][1]][0]);
	zs.push_back(table[triangles[k][1]][i][j]);
	xs.push_back(nT[triangles[k][2]][1]);
	ys.push_back(nT[triangles[k][2]][0]);
	zs.push_back(table[triangles[k][2]][i][j]);
	
	//std::cout<<"planar interp 2: " <<g*spectralFunction::PlanarInterpolate(n, T, xs, ys, zs)<<std::endl;
	return spectralFunction::PlanarInterpolate(n, T, xs, ys, zs);
}
double spectralFunction::spectral_function(double M, double Q, int k, double n, double T)
{
	int i = Mbin(M);
	int j = Qbin(Q);
	std::vector<std::vector<double>> points;
	
	
	
	//std::cout << "triangle: " << k << std::endl;
	//std::cout << "M: " << Mbins[i] << std::endl;
	//std::cout << "Q: " << Qbins[j] << std::endl;
	if (k == -1) 
	{

		return -3 * M * M / 12 / PI * e_sq;
	}
	std::vector<double> x1;
	x1.push_back(getSF(i-1,j-1,k,n,T));
	x1.push_back(getSF(i-1,j,k,n,T));
	points.push_back(x1);
	x1.clear();
	x1.push_back(getSF(i,j-1,k,n,T));
	x1.push_back(getSF(i,j,k,n,T));
	points.push_back(x1);
	
	/*std::vector<int> triangle = triangles[k];
	std::cout << "triangle:" << nT[triangle[0]][0] << " " << nT[triangle[0]][1] << std::endl;
	std::cout << "triangle:" << nT[triangle[1]][0] << " " << nT[triangle[1]][1] << std::endl;
	std::cout << "triangle:" << nT[triangle[2]][0] << " " << nT[triangle[2]][1] << std::endl;*/
	return g*spectralFunction::LinearInterpolate(M-Mbins[i-1], Q - Qbins[j - 1], Mbins[i]-Mbins[i-1], Qbins[j] - Qbins[j - 1], points);
}
int spectralFunction::upper_bin(const double* edges, int num, double x)
{
	// same result as the linear scans in Mbin/Qbin: first i >= 1 with x < edges[i], else num-1
	const double* it = std::upper_bound(edges + 1, edges + num, x);
	return it == edges + num ? num - 1 : int(it - edges);
}
void spectralFunction::cell_table(double n, double T, int k, std::vector<double>& out)
{
	out.assign((size_t)num_Mbins * num_Qbins, 0.);
	if (k < 0) return;
	// barycentric weights of (n,T) in triangle k; planar interpolation through the three table values is exactly this combination
	const int* t = &triangles[k][0];
	double x0 = nT[t[0]][1], y0 = nT[t[0]][0];
	double ax = nT[t[1]][1] - x0, ay = nT[t[1]][0] - y0, bx = nT[t[2]][1] - x0, by = nT[t[2]][0] - y0;
	double det = ax * by - ay * bx;
	double dx = n - x0, dy = T - y0;
	double w1 = (dx * by - dy * bx) / det, w2 = (ax * dy - ay * dx) / det, w0 = 1 - w1 - w2;
	for (int i = 0; i < num_Mbins; i++)
		for (int j = 0; j < num_Qbins; j++)
			out[(size_t)i * num_Qbins + j] = w0 * table[t[0]][i][j] + w1 * table[t[1]][i][j] + w2 * table[t[2]][i][j];
}
double spectralFunction::spectral_function_cell(const std::vector<double>& cell, int k, double M, double Q) const
{
	if (k == -1) return -3 * M * M / 12 / PI * e_sq;
	int i = upper_bin(Mbins, num_Mbins, M), j = upper_bin(Qbins, num_Qbins, Q);
	const double* r0 = &cell[(size_t)(i - 1) * num_Qbins];
	const double* r1 = &cell[(size_t)i * num_Qbins];
	double dx = Mbins[i] - Mbins[i - 1], dy = Qbins[j] - Qbins[j - 1];
	double x = M - Mbins[i - 1], y = Q - Qbins[j - 1];
	double a1 = (r0[j] - r0[j - 1]) * y / dy + r0[j - 1];
	double a2 = (r1[j] - r1[j - 1]) * y / dy + r1[j - 1];
	return g * ((a2 - a1) * x / dx + a1);
}
int spectralFunction::Triangle(double n, double T)
{
	// try the previous triangle first (consecutive cells have similar (n,T)); the planar interpolation is continuous
	// across shared edges, so any triangle containing the point gives the same value
	if (last_triangle >= 0 && isin(&triangles[last_triangle][0], n, T)) return last_triangle;

	//std::cout << "trying to find triangle in " <<triangles.size()<< std::endl;
	for (int i = 0; i < triangles.size() ; i++)
	{

		//std::cout<<i<<std::endl;
		if (spectralFunction::isin(&triangles[i][0], n, T))
		{
			last_triangle = i;
			return i;
		}
	}
	//std::cout << "not found triangle" << std::endl;

	return -1;

}
bool spectralFunction::isin(int* triangle, double n, double T)
{
	
	double a[2] = { nT[triangle[1]][1] - nT[triangle[0]][1], nT[triangle[1]][0] - nT[triangle[0]][0] };
	double b[2] = { nT[triangle[2]][1] - nT[triangle[0]][1], nT[triangle[2]][0] - nT[triangle[0]][0] };
	double x[2] = { n - nT[triangle[0]][1], T - nT[triangle[0]][0] };
	/*std::cout << "a vector:" << a[0] << " " << a[1]<<std::endl;
	std::cout << "b vector:" << b[0] << " " << b[1] << std::endl;
	std::cout << "x vector:" << x[0] << " " << x[1] << std::endl;
	std::cout << "tri vector:" << triangle[0] << " " << triangle[1] << std::endl;*/
	double b2 = (b[0] * b[0] + b[1] * b[1]);
	double a2 = (a[0] * a[0] + a[1] * a[1]);
	double ab = (a[0] * b[0] + a[1] * b[1]);
	double ax = (a[0] * x[0] + a[1] * x[1]);
	double bx = (x[0] * b[0] + x[1] * b[1]);


	double ca = (b2*ax-ab*bx)/(a2*b2-ab*ab);
	double cb = (a2 * bx - ab * ax) / (a2 * b2 - ab * ab);
	
	if (ca >= 0 && cb >= 0 && (ca + cb) < 1)
	{
		return 1;
	}
	//std::cout<<"no"<<std::endl;
	return 0;
}
