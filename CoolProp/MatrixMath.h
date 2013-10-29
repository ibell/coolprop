#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include <vector>
std::vector<double> linsolve_Gauss_Jordan(std::vector<std::vector<double> > const& A, std::vector<double> const& b);
//std::vector<double> dot_product(std::vector<std::vector<double> > const& A, std::vector<double> b);
std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > const& in);
std::vector<double> column(std::vector< std::vector<double> > const& in, unsigned int col);
#endif
