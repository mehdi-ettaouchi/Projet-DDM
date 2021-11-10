#ifndef GRADCONJ_H
#define GRADCONJ_H

#include <vector>
#include <fstream>

class GradConj
{
private:
	std::vector<std::vector<double>> A_;
	std::vector<double> b_;
	int k_, Nx_, Ny_;

public:
	//not sure about the consts
	GradConj(std::vector<std::vector<double>> A, std::vector<double> b, int Nx, int Ny);

	static std::vector<double> product(std::vector<std::vector<double>>, std::vector<double>, int, int);

	static std::vector<double> prod_scal(std::vector<double>, double);

	static double dot_product(std::vector<double>, std::vector<double>);

	static std::vector<double> sum(std::vector<double>, std::vector<double> y, int sign);

	static double norm(std::vector<double>);

	static double normMax(std::vector<double>);

	void Solve(int, std::vector<double> &);

	void MPI_Solve(int, std::vector<double> &);
};

#endif
