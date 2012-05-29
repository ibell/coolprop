#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "FluidClass.h"

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#define __ISWINDOWS__
#endif

#ifndef COOLPROPTOOLS_H
#define COOLPROPTOOLS_H

#include <string>
#include <cstdarg>

//missing string printf
//this is safe and convenient but not exactly efficient
inline std::string format(const char* fmt, ...)
{
    int size = 512;
    char* buffer = 0;
    buffer = new char[size];
    va_list vl;
    va_start(vl,fmt);
    int nsize = vsnprintf(buffer,size,fmt,vl);
    if(size<=nsize){//fail delete buffer and try again
        delete buffer; buffer = 0;
        buffer = new char[nsize+1];//+1 for /0
        nsize = vsnprintf(buffer,size,fmt,vl);
    }
    std::string ret(buffer);
    va_end(vl);
    delete buffer;
    return ret;
}

	#define OK 1
	#define FAIL 0
	
	double powInt(double x, int y);
	double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x);
	double CubicInterp( double x0, double x1, double x2, double x3, double f0, double f1, double f2, double f3, double x);
	int ValidNumber(double x);

#endif
