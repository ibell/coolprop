#include <vector>
#include <string>
#include <Eigen/Dense>

#ifndef SOLVERS_H
#define SOLVERS_H

class FuncWrapper1D
{
public:
	FuncWrapper1D(){};
	virtual ~FuncWrapper1D(){};
	virtual double call(double) = 0;
};

class FuncWrapperND
{
public:
	FuncWrapperND(){};
	virtual ~FuncWrapperND(){};
	virtual Eigen::VectorXd call(Eigen::VectorXd) = 0;// must be provided
	virtual Eigen::MatrixXd Jacobian(Eigen::VectorXd){return Eigen::MatrixXd();}; // optional
};

// Single-Dimensional solvers
double Brent(FuncWrapper1D *f, double a, double b, double macheps, double t, int maxiter, std::string *errstr);
double Secant(FuncWrapper1D *f, double x0, double dx, double tol, int maxiter, std::string *errstring);

// Multi-Dimensional solvers
Eigen::VectorXd NDNewtonRaphson_Jacobian(FuncWrapperND *f, Eigen::VectorXd x0, double tol, int maxiter, std::string *errstring);
#endif
