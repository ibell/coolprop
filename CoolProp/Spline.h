#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
/// A spline is a curve given by the form y = ax^3 + bx^2 + c*x + d
/// As there are 4 constants, 4 constraints are needed to create the spline.  These constraints could be the derivative or value at a point
/// Often, the value and derivative of the value are known at two points.

class SplineClass
{
protected:
	int Nconstraints;
	std::vector<std::vector<double> > A;
	std::vector<double> B;
public:
	double a,b,c,d;
	SplineClass();
	bool build(void);
	bool add_value_constraint(double x, double y);
	bool add_derivative_constraint(double x, double dydx);
	double evaluate(double x);
};

#endif