

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

	#if defined(_MSC_VER) || defined(__powerpc__)
	// Microsoft version of math.h doesn't include acosh or asinh, so we just define them here
	// Neither does the PPC version
	static double acosh(double x)
	{
 		return log(x + sqrt(x*x - 1.0) );
	}
	static double asinh(double value)
	{
		if(value>0){
			return log(value + sqrt(value * value + 1));
		}
		else{
			return -log(-value + sqrt(value * value + 1));
		}
	}
	#endif

	#include <algorithm> 
	#include <functional> 
	#include <cctype>
	#include <locale>
	#include <fstream>
	#include <cerrno>

	/// The following code for the trim functions was taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
	// trim from start
	inline std::string &strlstrip(std::string &s) {
			s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			return s;
	}
	// trim from end
	inline std::string &strrstrip(std::string &s) {
			s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
			return s;
	}
	// trim from both ends
	inline std::string &strstrip(std::string &s) {
			return strlstrip(strrstrip(s));
	}

	// Get all the contents of a file and dump into a STL string
	// Thanks to http://stackoverflow.com/questions/2602013/read-whole-ascii-file-into-c-stdstring
	std::string get_file_contents(const char *filename);

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

	std::vector<double> solve_cubic(double a, double b, double c, double d);

#endif
