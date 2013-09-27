#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include <vector>
std::vector<double> linsolve_Gauss_Jordan(std::vector<std::vector<double> > A, std::vector<double> b);
std::vector<double> dot_product(std::vector<std::vector<double> > A, std::vector<double> b);
std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > in);
std::vector<double> column(std::vector< std::vector<double> > in, unsigned int col);
#endif
