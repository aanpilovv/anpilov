#pragma once
#include <vector>
#include <cstdlib>
#include "OdeSolver.h"
#include "boundary_value_problem.h"
#include "Parameters.h"
#include "plots.h"
#include <iomanip>
#include <functional>
#include <complex>
#include "OdeSolver.h"

const double pi = 3.1415926538;


//template<class T>
class layer {
	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&, std::complex<double>, double)>>
		_equations;
	std::function<double(double)>& mu;
	std::function<double(double)>& rho;
	double kappa;
	std::vector<std::complex<double>> roots;
	std::vector<std::complex<double>> residual_set;

	std::vector<std::function<std::complex<double>(double, const std::vector<std::complex<double>>&)>> get_equation(std::complex<double> alpha, size_t size) const;

public:
	layer(double kkappa);
	layer(double kkappa, const std::vector<double>& mmu, const std::vector<double>& rrho);

	double defOnTop(std::complex<double> alpha) const;
	double findRoot(double a, double b, const std::function<double(double)>& f, double epsilon = 0.1e-6) const;
	std::vector<double> find_roots(double a, double b, const std::function<double(double)>& f, size_t n, double epsilon = 0.1e-6) const;
	std::map<double, std::vector<double>> dispersionalSet(double max_kappa, double step, size_t n, const std::function<double(double, double)>& f) const;
	std::complex<double> waves(double x1, const std::vector<std::complex<double>>& roots, const std::vector<std::complex<double>>& residuals) const;
	std::map<double, std::complex<double>> waveField(double a, double b, double step, const std::vector<std::complex<double>>& residuals) const;
	std::complex<double> residual(std::complex<double> alpha) const;
	std::vector<std::vector<double>> MatrixRho(size_t columns, size_t rows) const;
	std::vector<std::vector<double>> MatrixMu(size_t columns, size_t rows);
	std::vector<double> observed(const std::vector<double>& points) const;
private:
	std::vector<std::complex<double>> getRoots();
	std::vector<std::complex<double>> residualSet();
};

