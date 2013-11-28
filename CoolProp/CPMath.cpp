
#include <stdio.h>

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>

#include <stdlib.h>

#include "MatrixMath.h"

void print_string(std::string str) {
  std::cout << str << std::endl;
}
double get_element() {
	return (rand() % 10000)/1000.0;
}
std::string bool_to_string(bool b)
{
  return b ? std::string("true") : std::string("false");
}
//
//int main( int argc, const char* argv[] ){
//
//	print_string(std::string(" "));
//
//	std::vector<double> tmpVec;
//	int dim = 3;
//
//	std::vector<std::vector<double> > A;
//	for (int i = 0; i < dim; i++) {
//		tmpVec.clear();
//		for (int j = 0; j < dim; j++) {
//			tmpVec.push_back(get_element());
//		}
//		A.push_back(tmpVec);
//	}
//	print_string(vec_to_string(A));
//	std::vector<double> A2 = get_col(A,2);
//
//	std::vector<std::vector<double> > B;
//	for (int i = 0; i < dim; i++) {
//		tmpVec.clear();
//		for (int j = 0; j < dim+1; j++) {
//			tmpVec.push_back(get_element());
//		}
//		B.push_back(tmpVec);
//	}
//	print_string(vec_to_string(B));
//	std::vector<double> B2 = get_col(B,1);
//
//	std::vector<std::vector<double> > C;
//	for (int i = 0; i < dim-1; i++) {
//		tmpVec.clear();
//		for (int j = 0; j < dim; j++) {
//			tmpVec.push_back(get_element());
//			if (i==dim-2&&j==dim-1) {tmpVec.pop_back(); }
//		}
//		C.push_back(tmpVec);
//	}
//	print_string(vec_to_string(C));
//
//	std::vector<std::vector<double> > D;
//	tmpVec.clear();
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(2.0);
//	tmpVec.push_back(-2.0);
//	D.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(1.0);
//	D.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(1.0);
//	D.push_back(tmpVec);
//	print_string(vec_to_string(D));
//
//	std::vector<std::vector<double> > E;
//	tmpVec.clear();
//	tmpVec.push_back(8.0);
//	E.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(-1.0);
//	E.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(-3.0);
//	E.push_back(tmpVec);
//	print_string(vec_to_string(E));
//
//	std::vector<std::vector<double> > F;
//	tmpVec.clear();
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(0.0);
//	F.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(0.0);
//	F.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(0.0);
//	tmpVec.push_back(1.0);
//	F.push_back(tmpVec);
//	print_string(vec_to_string(F));
//
//	std::vector<double> G;
//	G.clear();
//	G.push_back(8.0);
//	G.push_back(-1.0);
//	G.push_back(-3.0);
//	print_string(vec_to_string(G));
//
//	std::vector<std::vector<double> > H;
//	tmpVec.clear();
//	tmpVec.push_back(1.0);
//	tmpVec.push_back(3.0);
//	tmpVec.push_back(1.0);
//	H.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(-1.0);
//	tmpVec.push_back(2.0);
//	tmpVec.push_back(0.0);
//	H.push_back(tmpVec);
//	tmpVec.clear();
//	tmpVec.push_back(2.0);
//	tmpVec.push_back(11.0);
//	tmpVec.push_back(3.0);
//	H.push_back(tmpVec);
//	print_string(vec_to_string(H));
//
//	print_string(std::string("--------------------"));
//	print_string(std::string(" Testing functions: "));
//	print_string(std::string(" "));
//	print_string(std::string(" is_squared: A,B,C"));
//	print_string(std::string(bool_to_string(is_squared(A)))+", "+bool_to_string(is_squared(B))+", "+bool_to_string(is_squared(C)));
//	print_string(std::string(" "));
//	print_string(std::string(" make_squared: A,B,C"));
//	print_string(std::string(vec_to_string(make_squared(A))));
//	print_string(std::string(vec_to_string(make_squared(B))));
//	print_string(std::string(vec_to_string(make_squared(C))));
//	print_string(std::string(" is_squared: makeSquared(A),makeSquared(B),makeSquared(C)"));
//	print_string(std::string(bool_to_string(is_squared(make_squared(A))))+", "+bool_to_string(is_squared(make_squared(B)))+", "+bool_to_string(is_squared(make_squared(C))));
//	print_string(std::string(" "));
//	print_string(std::string(" multiply: A*B"));
//	print_string(std::string(vec_to_string(multiply(A,B))));
//	print_string(std::string(" multiply: A*B[1]"));
//	print_string(std::string(vec_to_string(multiply(A,B2))));
//	print_string(std::string(" multiply: A[2]*B[1]"));
//	print_string(std::string(vec_to_string(multiply(A2,B2))));
//	print_string(std::string(" dot_product: A[2]*B[1]"));
//	print_string(std::string(vec_to_string(dot_product(A2,B2))));
//	print_string(std::string(" "));
//	print_string(std::string(" transpose: A"));
//	print_string(std::string(vec_to_string(transpose(A))));
//	print_string(std::string(" "));
//	print_string(std::string(" invert: A"));
//	print_string(std::string(vec_to_string(invert(A))));
//	print_string(std::string(" invert: D"));
//	print_string(std::string(vec_to_string(invert(D))));
//	print_string(std::string(" invert: H"));
//	print_string(std::string(vec_to_string(invert(H))));
//	print_string(std::string(" "));
//	print_string(std::string(" solve: D*x=E"));
//	print_string(vec_to_string(linsolve(D,E)));
//	print_string(std::string(" solve: D*x=F"));
//	print_string(vec_to_string(linsolve(D,F)));
//	print_string(std::string(" solve: D*x=G"));
//	print_string(vec_to_string(linsolve(D,G)));
//
//	//print_string(std::string(vec_to_string(transpose(A))));
//
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//	print_string(std::string(" "));
//
//
//}
