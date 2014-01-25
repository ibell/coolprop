#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include <vector>
#include <string>

/// Publish the linear algebra solver
            std::vector<double>   linsolve(std::vector<std::vector<double> > const& A,             std::vector<double>   const& b);
std::vector<std::vector<double> > linsolve(std::vector<std::vector<double> > const& A, std::vector<std::vector<double> > const& B);

/// Some shortcuts and regularly needed operations
std::size_t         num_rows  (std::vector<std::vector<double> > const& in);
std::size_t         num_cols  (std::vector<std::vector<double> > const& in);
std::size_t         max_cols  (std::vector<std::vector<double> > const& in);
std::vector<double> get_row   (std::vector<std::vector<double> > const& in, size_t row);
std::vector<double> get_col   (std::vector<std::vector<double> > const& in, size_t col);
bool                is_squared(std::vector<std::vector<double> > const& in);
std::vector<std::vector<double> > make_squared(std::vector<std::vector<double> > const& in);

/// Define some basic math operations for vectors
                        double    multiply(            std::vector<double>   const& A,             std::vector<double>   const& B);
            std::vector<double>   multiply(std::vector<std::vector<double> > const& A,             std::vector<double>   const& B);
std::vector<std::vector<double> > multiply(std::vector<std::vector<double> > const& A, std::vector<std::vector<double> > const& B);

double              dot_product(std::vector<double> const& a, std::vector<double> const& b);
std::vector<double> cross_product(std::vector<double> const& a, std::vector<double> const& b);

std::vector<std::vector<double> > transpose(std::vector<std::vector<double> > const& in);
std::vector<std::vector<double> >    invert(std::vector<std::vector<double> > const& in);

std::string vec_to_string(                        double    const& a);
std::string vec_to_string(            std::vector<double>   const& a);
std::string vec_to_string(std::vector<std::vector<double> > const& A);

std::string vec_to_string(            std::vector<double>   const& a, const char *fmt);
std::string vec_to_string(std::vector<std::vector<double> > const& A, const char *fmt);
#endif
