
#include "MatrixMath.h"
#include "math.h"
#include <iostream>

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

std::vector<double> linsolve_Gauss_Jordan(std::vector<std::vector<double> > A, std::vector<double> b)
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
		divide_row_by(&Ab,col,Ab[col][col]);

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