
#include <string>
#include <vector>
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


std::vector< std::vector<double> > transpose(std::vector< std::vector<double> > const& in){
	unsigned int sizeX = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeX));
	unsigned int sizeY    = in[0].size();
	unsigned int sizeYOld = sizeY;
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeY));
	std::vector< std::vector<double> > out(sizeY,std::vector<double>(sizeX));
	for (unsigned int i = 0; i < sizeX; ++i){
		sizeY = in[i].size();
		if (sizeY!=sizeYOld) throw ValueError(format("You have to provide a rectangular matrix: %d is not equal to %d. ",sizeY,sizeYOld));
		for (unsigned int j = 0; j < sizeY; ++j){
			out[j][i] = in[i][j];
		}
	}
	return out;
}


std::vector<double> column(std::vector< std::vector<double> > const& in, unsigned int col){
	unsigned int sizeX = in.size();
	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeX));
	unsigned int sizeY    = in[0].size();
	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not a valid. ",sizeY));
	std::vector<double> out;
	for (unsigned int i = 0; i < sizeX; ++i){
		sizeY = in[i].size();
		if (sizeY-1<col) throw ValueError(format("Your matrix does not have enough entries in row %d, last index %d is less than %d. ",i,sizeY-1,col));
		out.push_back(in[i][col]);
	}
	return out;
}

