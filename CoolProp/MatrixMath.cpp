
#include <string>
#include <sstream>
#include <vector>
#include <numeric>
#include <math.h>
#include "CoolPropTools.h"
#include "CPExceptions.h"

#include "MatrixMath.h"

/*
Owe a debt of gratitude to http://sole.ooz.ie/en - very clear treatment of GJ
*/
void swap_rows(std::vector<std::vector<double> > *A, unsigned int row1, unsigned int row2)
{
	for (unsigned int col = 0; col < (*A)[0].size(); col++){
		std::swap((*A)[row1][col],(*A)[row2][col]);
	}
}
void subtract_row_multiple(std::vector<std::vector<double> > *A, unsigned int row, double multiple, unsigned int pivot_row)
{
	for (unsigned int col = 0; col < (*A)[0].size(); col++){
		(*A)[row][col] -= multiple*(*A)[pivot_row][col];
	}
}
void divide_row_by(std::vector<std::vector<double> > *A, unsigned int row, double value)
{
	for (unsigned int col = 0; col < (*A)[0].size(); col++){
		(*A)[row][col] /= value;
	}
}

unsigned int get_pivot_row(std::vector<std::vector<double> > *A, unsigned int col)
{
	int index = col;
	double max = 0, val;

	for (unsigned int row = col; row < (*A).size(); row++)
	{
		val = (*A)[row][col];
		if (fabs(val) > max)
		{
			max = fabs(val);
			index = row;
		}
	}
	return index;
}

std::vector<double> linsolve_Gauss_Jordan(std::vector<std::vector<double> > const& A, std::vector<double> const& b)
{
	std::vector<std::vector<double> > Ab;
	std::vector<double> x;
	unsigned int pivot_row;
	double pivot_element;

	Ab.resize(A.size(), std::vector<double>(A[0].size()+1, 0));
	x.resize(A.size(), 0);

	unsigned int Nrow = A.size();
	unsigned int Ncol = A[0].size();

	// Build the augmented matrix
	for (unsigned int row = 0; row < Nrow; row++)
	{
		for (unsigned int col = 0; col < Ncol; col++)
		{
			Ab[row][col] = A[row][col];
		}
		Ab[row][Ncol] = b[row];
	}

	for (unsigned int col = 0; col < Ncol; col++)
	{
		// Find the pivot value
		pivot_row = get_pivot_row(&Ab, col);

		if (pivot_row>=col){
			// Swap pivot row and current row
			swap_rows(&Ab, col, pivot_row);
		}
		// Get the pivot element
		pivot_element = Ab[col][col];
		// Divide the pivot row by the pivot element
		divide_row_by(&Ab,col,pivot_element);

		if (col < Nrow-1)
		{
			// All the rest of the rows, subtract the value of the [r][c] combination
			for (unsigned int row = col + 1; row < Ab.size(); row++)
			{
				subtract_row_multiple(&Ab,row,Ab[row][col],col);
			}
		}
	}
	for (int col = Ncol - 1; col > 0; col--)
	{
		for (int row = col - 1; row >=0; row--)
		{
			subtract_row_multiple(&Ab,row,Ab[row][col],col);
		}
		// Set the output value
		x[col] = Ab[col][Ncol];
	}
	// The last output value
	x[0] = Ab[0][Ncol];


	return x;
}

std::vector<double> linsolve(std::vector<std::vector<double> > const& A, std::vector<double> const& b){
	return linsolve_Gauss_Jordan(A, b);
}



/// Some shortcuts and regularly needed operations
std::size_t         num_rows  (std::vector<std::vector<double> > const& in){ return in.size(); }
std::size_t         num_cols  (std::vector<std::vector<double> > const& in){
	if (num_rows(in)>0) {
		if (is_squared(in)) {
			return in[0].size();
		} else {
			return max_cols(in);
		}
	} else {
		return 0;
	}
}
std::size_t         max_cols  (std::vector<std::vector<double> > const& in){
	std::size_t cols = 0;
	std::size_t col  = 0;
	for (std::size_t i = 0; i < in.size(); i++) {
		col = in[i].size();
		if (cols<col) {cols = col;}
    }
	return cols;
}
std::vector<double> get_row(std::vector< std::vector<double> > const& in, unsigned int row) { return in[row]; }
std::vector<double> get_col(std::vector< std::vector<double> > const& in, unsigned int col) {
	std::size_t sizeX  = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeX));
	unsigned int sizeY = in[0].size();
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeY));
	std::vector<double> out;
	for (std::size_t i = 0; i < sizeX; i++) {
		sizeY = in[i].size();
		if (sizeY-1<col) throw ValueError(format("Your matrix does not have enough entries in row %d, last index %d is less than %d. ",i,sizeY-1,col));
		out.push_back(in[i][col]);
	}
	return out;
}
bool                is_squared(std::vector<std::vector<double> > const& in){
	std::size_t cols = max_cols(in);
	if (cols!=num_rows(in)) { return false;}
	else {
		for (std::size_t i = 0; i < in.size(); i++) {
			if (cols!=in[i].size()) {return false; }
		}
	}
	return true;
}
std::vector<std::vector<double> > make_squared(std::vector<std::vector<double> > const& in){
	std::size_t cols   = max_cols(in);
	std::size_t rows   = num_rows(in);
	std::size_t maxVal = 0;
	std::vector<std::vector<double> > out;
	            std::vector<double>   tmp;

	if (cols>rows) {maxVal = cols; }
	else           {maxVal = rows; }
	out.clear();
	for (std::size_t i = 0; i < in.size(); i++) {
		tmp.clear();
		for (std::size_t j = 0; j < in[i].size(); j++) {
			tmp.push_back(in[i][j]);
		}
		while (maxVal>tmp.size()) {
			tmp.push_back(0.0);
		}
		out.push_back(tmp);
    }
	// Check rows
	tmp.clear();
	tmp.resize(maxVal,0.0);
	while (maxVal>out.size()) {
		out.push_back(tmp);
	}
	return out;
}

                        double    multiply(            std::vector<double>   const& a,             std::vector<double>   const& b){
    return dot_product(a,b);

}
//            std::vector<double>   multiply(std::vector<std::vector<double> > const& A,             std::vector<double>   const& b){
//    std::vector<std::vector<double> > tmpMat;
//                std::vector<double>   tmpVec;
//    tmpMat.clear();
//    for (size_t i = 0; i < b.size(); i++){
//		tmpVec.clear();
//		tmpVec.push_back(b[i]);
//		tmpMat.push_back(tmpVec);
//	}
//    return multiply(A,tmpMat);
//}

std::vector<std::vector<double> > multiply(std::vector<std::vector<double> > const& A, std::vector<std::vector<double> > const& B){
	if (num_cols(A) != num_rows(B)){
		throw ValueError(format("You have to provide matrices with the same columns and rows: %d is not equal to %d. ",num_cols(A),num_rows(B)));
	}
	size_t rows = num_rows(A);
	size_t cols = num_cols(B);
	double tmp;
	std::vector<std::vector<double> > outVec;
		        std::vector<double>   tmpVec;
    outVec.clear();
	for (size_t i = 0; i < rows; i++){
		tmpVec.clear();
		for (size_t j = 0; j < cols; j++){
			tmp = 0.0;
			for (size_t k = 0; k < num_cols(A); k++){
				tmp += A[i][k] * B[k][j];
			}
			tmpVec.push_back(tmp);
		}
		outVec.push_back(tmpVec);
	}
	return outVec;
}

double              dot_product(std::vector<double> const& a, std::vector<double> const& b){
	if (a.size()==b.size()){
		return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
	}
	throw ValueError(format("You have to provide vectors with the same length: %d is not equal to %d. ",a.size(),b.size()));
}

std::vector<double> cross_product(std::vector<double> const& a, std::vector<double> const& b){
	throw NotImplementedError("The cross product function has not been implemented, yet");
}

std::vector< std::vector<double> > transpose(std::vector<std::vector<double> > const& in){
	size_t sizeX = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeX));
	size_t sizeY    = in[0].size();
	size_t sizeYOld = sizeY;
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeY));
	std::vector< std::vector<double> > out(sizeY,std::vector<double>(sizeX));
	for (size_t i = 0; i < sizeX; ++i){
		sizeY = in[i].size();
		if (sizeY!=sizeYOld) throw ValueError(format("You have to provide a rectangular matrix: %d is not equal to %d. ",sizeY,sizeYOld));
		for (size_t j = 0; j < sizeY; ++j){
			out[j][i] = in[i][j];
		}
	}
	return out;
}

std::vector< std::vector<double> >    invert(std::vector<std::vector<double> > const& in){
	throw NotImplementedError("The invert function has not been implemented, yet");
}

std::string vec_to_string(            std::vector<double>   const& a) {
	if (a.size()<1) {
		return std::string("");
	} else {
		std::stringstream out;
		out << a[0];
		for (size_t j = 1; j < a.size(); j++) {
			out << ' ' << a[j];
		}
		out << '\n';
		return out.str();
	}
}

std::string vec_to_string(std::vector<std::vector<double> > const& A) {
	std::stringstream out;
	for (size_t j = 0; j < A.size(); j++) {
		out << vec_to_string(A[j]);
    }
    return out.str();
}

