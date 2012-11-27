#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <cstdio>
#include <cstdarg>
#include <stdlib.h>
#include "math.h"
#include "stdio.h"
#include <string.h>
#include "CoolProp.h"

std::string format(const char* fmt, ...)
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

std::vector<std::string> strsplit(std::string s, char del)
{
	int iL = 0, iR = 0, N;
	N = s.size();
	std::vector<std::string> v;
	// Find the first instance of the delimiter
	iR = s.find_first_of(del);
	// Delimiter not found, return the same string again
	if (iR<0){
		v.push_back(s);
		return v;
	}
	while (iR != N-1)
	{
		v.push_back(s.substr(iL,iR-iL));	
		iL = iR;
		iR = s.find_first_of(del,iR+1);
		// Move the iL to the right to avoid the delimiter
		iL += 1;
		if (iR == std::string::npos)
		{
			v.push_back(s.substr(iL,N-iL));
			return v;
		}
	}
}
    
double interp1d(std::vector<double> *x, std::vector<double> *y, double x0)
{
	unsigned int i,L,R,M;
	L=0;
	R=(*x).size()-1;
	M=(L+R)/2;
	// Use interval halving to find the indices which bracket the density of interest
	while (R-L>1)
	{
		if (x0 >= (*x)[M])
		{ L=M; M=(L+R)/2; continue;}
		if (x0 < (*x)[M])
		{ R=M; M=(L+R)/2; continue;}
	}
	i=L;
	if (i<(*x).size()-2)
    {
		// Go "forwards" with the interpolation range
		return QuadInterp((*x)[i],(*x)[i+1],(*x)[i+2],(*y)[i],(*y)[i+1],(*y)[i+2],x0);
    }
    else
    {
        // Go "backwards" with the interpolation range
		return QuadInterp((*x)[i],(*x)[i-1],(*x)[i-2],(*y)[i],(*y)[i-1],(*y)[i-2],x0);
    }
}

static int isNAN(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	return x != x;
}
static int isINFINITY(double x)
{
	// recommendation from http://www.devx.com/tips/Tip/42853
	if ((x == x) && ((x - x) != 0.0))
		return 1;//return (x < 0.0 ? -1 : 1); // This will tell you whether positive or negative infinity
	else
		return 0;
}
int ValidNumber(double x)
{
	if (!isNAN(x) && !isINFINITY(x))
		return 1;
	else
		return 0;
}
double powInt(double x, int y)
{
    // Raise a double to an integer power
    // Overload not provided in math.h
    int i;
    double product=1.0;
    double x_in;
    int y_in;
    
    if (y==0)
    {
        return 1.0;
    }
    
    if (y<0)
    {
        x_in=1/x;
        y_in=-y;
    }
    else
    {
        x_in=x;
        y_in=y;
    }

    if (y_in==1)
    {
        return x_in;
    }    
    
    product=x_in;
    for (i=1;i<y_in;i++)
    {
        product=product*x_in;
    }
    return product;
}

double QuadInterp(double x0, double x1, double x2, double f0, double f1, double f2, double x)
{
    /* Quadratic interpolation.  
    Based on method from Kreyszig, 
    Advanced Engineering Mathematics, 9th Edition 
    */
    double L0, L1, L2;
    L0=((x-x1)*(x-x2))/((x0-x1)*(x0-x2));
    L1=((x-x0)*(x-x2))/((x1-x0)*(x1-x2));
    L2=((x-x0)*(x-x1))/((x2-x0)*(x2-x1));
    return L0*f0+L1*f1+L2*f2;
}
double CubicInterp( double x0, double x1, double x2, double x3, double f0, double f1, double f2, double f3, double x)
{
	/*
	Lagrange cubic interpolation as from
	http://nd.edu/~jjwteach/441/PdfNotes/lecture6.pdf
	*/
	double L0,L1,L2,L3;
	L0=((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
	L1=((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
	L2=((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
	L3=((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));
	return L0*f0+L1*f1+L2*f2+L3*f3;
}