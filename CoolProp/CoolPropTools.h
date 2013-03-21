

#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

	#define _CRT_SECURE_NO_WARNINGS

	#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
	  #define __ISWINDOWS__
	#elif __APPLE__
      #define __ISAPPLE__
    #elif __linux
      #define __ISLINUX__
    #endif

	#include <string>
	#include <vector>
	#include <cmath>

	#ifndef M_PI
	#define M_PI 3.14159265358979323846
	#endif

	#ifdef HUGE_VAL
	#define _HUGE HUGE_VAL
	#else
		// GCC Version of huge value macro
		#ifdef HUGE 
		#define _HUGE HUGE
		#endif
	#endif

    //missing string printf
    std::string format(const char* fmt, ...);
	// Missing string split - like in Python
	std::vector<std::string> strsplit(std::string s, char del);

	#define OK 1
	#define FAIL 0
	
	void MatInv_2(double A[2][2] , double B[2][2]);

	double root_sum_square(std::vector<double> x);
	double interp1d(std::vector<double> *x, std::vector<double> *y, double x0);
	double powInt(double x, int y);
	double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);
	double CubicInterp( double x0, double x1, double x2, double x3, double f0, double f1, double f2, double f3, double x);
	bool ValidNumber(double x);


#endif
