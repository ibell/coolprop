//You can include any C libraries that you normally use
#include "math.h"
#include "CoolProp.h"
#include "mex.h"   //--This one is required


// 
// In MATLAB, to get MEXing working, install visual studio express 2008
// then type 
//         
// >> mex -setup

// Please choose your compiler for building external interface (MEX) files: 
//  
// Would you like mex to locate installed compilers [y]/n? y
//  
// Select a compiler: 
// [1] Lcc-win32 C 2.4.1 in C:\MATLAB\R2010a\sys\lcc 
// [2] Microsoft Visual C++ 2008 Express in C:\Program Files\Microsoft Visual Studio 9.0 
//  
// [0] None 
//  
// Compiler: 2
//  
// Please verify your choices: 
//  
// Compiler: Microsoft Visual C++ 2008 Express  
// Location: C:\Program Files\Microsoft Visual Studio 9.0 
//  
// Are these correct [y]/n? y
//  
// *************************************************************************** 
//   Warning: MEX-files generated using Microsoft Visual C++ 2008 require 
//            that Microsoft Visual Studio 2008 run-time libraries be  
//            available on the computer they are run on. 
//            If you plan to redistribute your MEX-files to other MATLAB 
//            users, be sure that they have the run-time libraries. 
// *************************************************************************** 
//  
// Trying to update options file: C:\Documents and Settings\ibell\Application Data\MathWorks\MATLAB\R2010a\mexopts.bat 
// From template:              C:\MATLAB\R2010a\bin\win32\mexopts\msvc90freeopts.bat 
//  
// Done . . . 
//  
// ************************************************************************** 
//   Warning: The MATLAB C and Fortran API has changed to support MATLAB 
//            variables with more than 2^32-1 elements.  In the near future 
//            you will be required to update your code to utilize the new 
//            API. You can find more information about this at: 
//            http://www.mathworks.com/support/solutions/en/data/1-5C27B9/?solution=1-5C27B9 
//            Building with the -largeArrayDims option enables the new API. 
// ************************************************************************** 
                   
// At the prompt type
// >> mex bencode.c
// Compiles bencode.c to bencode.mexw32

// To run:
// >> hotty=[3 2; 4 5; 7 6]; 
// >> x=bencode(hotty)


//Function prototype
double mean(double *X, int u, int v);

void inv(double *S, double *Sinv)
{
    //Calculate the inverse
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int npol,Nel;
    int m,n,u,v,k,Output_len,Prop1_len,Prop2_len,Ref_len;
    int *c,status;
    double *d,Prop1,Prop2,x;
	char Output,Name1,Name2,*Ref;
    mxArray *cMat[1];
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);

	// Figure out the lengths of the input strings
	//Output_len=(mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
	//Prop1_len=(mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
	//Prop2_len=(mxGetM(prhs[4]) * mxGetN(prhs[4])) + 1;
	
	//status = mxGetChar(prhs[0], Prop1, Prop1_len);
    //status = mxGetChar(prhs[0], Prop2, Prop2_len);

	// Get the refrigerant (it is a string)
	Ref_len=(mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
    Ref = mxCalloc(Ref_len, sizeof(char));
    status = mxGetString(prhs[5], Ref, Ref_len);

	Output = mxGetScalar(prhs[0]);
	Name1 = mxGetScalar(prhs[1]);
	Name2 = mxGetScalar(prhs[3]);

    Prop1 = mxGetScalar(prhs[2]);
	Prop2 = mxGetScalar(prhs[4]);

    x = Props(Output,Name1,Prop1,Name2,Prop2,Ref);
     
    *mxGetPr(plhs[0])=x;
}
