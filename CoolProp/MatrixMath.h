#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include <vector>
std::vector<double> linsolve_Gauss_Jordan(std::vector<std::vector<double> > A, std::vector<double> b);
std::vector<double> dot_product(std::vector<std::vector<double> > A, std::vector<double> b);
#endif