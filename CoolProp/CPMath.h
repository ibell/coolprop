//
//#ifndef CPMATH_H
//#define CPMATH_H
//
//#include "CPExceptions.h"
//
//#include <fstream>
//#include <iostream>
//#include <iterator>
//#include <string>
//#include <valarray>
//#include <vector>
//
///** Matrix2d class
//  * I tried to make some more general approach to all our matrix math mess.
//  * I hope to use this class as a basic wrapper for all matrix operations and
//  * some related solvers. However, I am not a C++ expert and I am sure that
//  * there is room for optimisation.
//  * I have read about using a "const" and valarray to speed things up, but
//  * I do not understand the concepts. Please help me if you can.
//  */
//
//template <typename T> class Matrix2d{
//public:
//	/// First some simple constructors
//	Matrix2d();                                             // Unknown size
//	Matrix2d(std::size_t num_rows, std::size_t num_cols);   // Known size
//	Matrix2d(            std::vector<T>   const& entries); // Known data
//	Matrix2d(std::vector<std::vector<T> > const& entries); // Known data
//
//	/// Some more shortcuts and regularly needed operations
//	std::size_t         num_rows();
//	std::size_t         num_cols();
//	std::size_t         max_cols();
//	std::vector<T>      get_row(std::size_t row);
//	std::vector<T>      get_col(std::size_t col);
//	bool                is_squared();
//
//	/// Data conversion functions
//	std::vector< std::vector<T> > to_vector_2d();
//	std::vector<T>                to_vector_1d();
//	Matrix2d<T>                   to_matrix();
//
//
//	std::vector< std::vector<T> > transpose();
//	std::vector< std::vector<T> >    invert();
//
//	static std::vector< std::vector<T> > transpose(std::vector<std::vector<T> > const& A);
//	static std::vector< std::vector<T> >    invert(std::vector<std::vector<T> > const& A);
//
//	/// Linear algebra solver (implemented is the Gauss-Jordan)
//	static std::vector<T> linsolve(std::vector<std::vector<T> > const& A, std::vector<T> const& b);
//	//static    Matrix2d<T> linsolve(Matrix2d<T> const& A, std::vector<T> const& b);
//	//static    Matrix2d<T> linsolve(Matrix2d<T> const& A, Matrix2d<T> const& b);
//	/// Fill a matrix with zeros to make it squared
//	//static    Matrix2d<T> to_squared (std::vector<std::vector<T> > const& A);
//	//static    Matrix2d<T> to_squared (Matrix2d<T> const& A);
//
//	/// Functions used during debugging
//	       std::string to_string();
//	static std::string to_string(std::vector<std::vector<T> > const& A);
//	static std::string to_string(            std::vector<T>   const& A);
//           void     print_string();
//	static void     print_string(std::vector<std::vector<T> > const& A);
//	static void     print_string(            std::vector<T>   const& A);
//
//protected:
//	bool _squared();
//	bool squared;
//	std::vector< std::vector<T> > data;
//};
//
///// First some simple constructors
//template<class T> Matrix2d<T>::Matrix2d(){
//	data = new std::vector<T>();
//	squared = true;
//}
//
//template<class T> Matrix2d<T>::Matrix2d(std::size_t rows, std::size_t columns){
//	data.clear();
//	for (size_t i = 0; i < rows; i++) {
//		data.push_back(new std::vector<T>(columns, NULL) );
//	}
//	squared = true;
//}
//
//template<class T> Matrix2d<T>::Matrix2d(            std::vector<T>   const& entries){
//	data.clear();
//	data.push_back(entries);
//	squared = _squared();
//}
//
//template<class T> Matrix2d<T>::Matrix2d(std::vector<std::vector<T> > const& entries){
//	data = entries;
//	squared = _squared();
//}
//
///// Some more shortcuts and regularly needed operations
//template<class T> std::size_t         Matrix2d<T>::num_rows(){ return data.size(); }
//template<class T> std::size_t         Matrix2d<T>::num_cols(){
//	if (this->num_rows()>0) {
//		if (this->is_squared()) {
//			return data[0].size();
//		} else {
//			return this->max_cols();
//		}
//	} else {
//		return 0;
//	}
//}
//template<class T> std::size_t         Matrix2d<T>::max_cols(){
//	std::size_t cols = 0;
//	std::size_t col  = 0;
//	for (std::size_t i = 0; i < data.size(); i++) {
//		col = data[i].size();
//		if (cols<col) {cols = col;}
//    }
//	return cols;
//}
//
//template<class T> std::vector<T>      Matrix2d<T>::get_row(std::size_t row) { return data[row]; }
//template<class T> std::vector<T>      Matrix2d<T>::get_col(std::size_t col) {
//	std::size_t sizeX  = this->num_rows();
//	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeX));
//	unsigned int sizeY = this->num_cols();
//	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeY));
//	std::vector<T> out;
//	for (std::size_t i = 0; i < data.size(); i++) {
//		sizeY = data[i].size();
//		if (sizeY-1<col) throw ValueError(format("Your matrix does not have enough entries in row %d, last index %d is less than %d. ",i,sizeY-1,col));
//		out.push_back(data[i][col]);
//	}
//	return out;
//}
//template<class T> bool               Matrix2d<T>::is_squared(){ return squared; }
//
///// Data conversion functions
//template<class T> std::vector< std::vector<T> > Matrix2d<T>::to_vector_2d(){return data; }
//template<class T> std::vector<T>                Matrix2d<T>::to_vector_1d(){return this->get_col(0); }
//template<class T> Matrix2d<T>                   Matrix2d<T>::to_matrix()   {return this; }
//template<class T> std::vector< std::vector<T> > Matrix2d<T>::transpose()   {
//	std::size_t sizeX = data.size();
//	if (sizeX<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeX));
//	std::size_t sizeY    = data[0].size();
//	std::size_t sizeYOld = sizeY;
//	if (sizeY<1) throw ValueError(format("You have to provide values, a vector length of %d is not valid. ",sizeY));
//	std::vector< std::vector<T> > out(sizeY,std::vector<T>(sizeX));
//	for (std::size_t i = 0; i < sizeX; ++i){
//		sizeY = data[i].size();
//		if (sizeY!=sizeYOld) throw ValueError(format("You have to provide a rectangular matrix: %d is not equal to %d. ",sizeY,sizeYOld));
//		for (std::size_t j = 0; j < sizeY; ++j){
//			out[j][i] = data[i][j];
//		}
//	}
//	return out;
//}
//template<class T> std::vector< std::vector<T> > Matrix2d<T>::invert()      {
//	throw NotImplementedError("The invert function has not been implemented, yet");
//	return this->to_vector_2d();
//}
//
//template<class T> static std::vector< std::vector<T> > Matrix2d<T>::transpose (std::vector<std::vector<T> > const& A){
//	Matrix2d mat = new Matrix2d(A);
//	return mat.transpose();
//}
//template<class T> static std::vector< std::vector<T> > Matrix2d<T>::invert    (std::vector<std::vector<T> > const& A){
//	Matrix2d mat = new Matrix2d(A);
//	return mat.invert();
//}
//
//
//template<class T> static std::vector<T> Matrix2d<T>::linsolve(std::vector<std::vector<T> > const& A, std::vector<T> const& b);
////template<class T> static    Matrix2d<T> Matrix2d<T>::squared (std::vector<std::vector<T> > const& A);
//
//
//
//
//template<class T> std::string Matrix2d<T>::to_string() {
//	std::stringstream out;
//    for (size_t i = 0; i < this->num_rows(); i++) {
//	    std::vector<T> row = this->get_row(i);
//	    out << row[0];
//		for (size_t j = 1; j < this->num_cols(); j++) {
//			out << ' ' << row[j];
//		}
//	out << '\n';
//    }
//    return out.str();
//}
//template<class T> static std::string Matrix2d<T>::to_string(std::vector<std::vector<T> > const& A){
//	Matrix2d mat = new Matrix2d(A);
//	return mat.to_string();
//}
//template<class T> static std::string Matrix2d<T>::to_string(            std::vector<T>   const& A){
//	Matrix2d mat = new Matrix2d(A);
//	return mat.to_string();
//}
//
//template<class T> void Matrix2d<T>::print_string() {
//  std::cout << this->to_string() << std::endl;
//}
//template<class T> static void Matrix2d<T>::print_string(std::vector<std::vector<T> > const& A){
//	Matrix2d mat = new Matrix2d(A);
//	mat.print_string();
//}
//template<class T> static void Matrix2d<T>::print_string(            std::vector<T>   const& A){
//	Matrix2d mat = new Matrix2d(A);
//	mat.print_string();
//}
//
//
//template<class T> bool Matrix2d<T>::_squared(){
//	std::size_t cols = this->max_cols();
//	for (std::size_t i = 0; i < data.size(); i++) {
//		if (cols!=data[i].size()) {return false; }
//    }
//	return true;
//}
//
//#endif
//
//
