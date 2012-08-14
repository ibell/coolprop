#include "math.h"
#include "mex.h"   //--This one is required

// 
// In MATLAB, to get MEXing working, install visual studio express 2008
// then type 
//         
// >> mex -setup
//
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


// Prototype for the Props function to be called.  Can't use the CoolProp.h header because there are a lot of
// c++ parts in the header that cannot be easily hidden when compiling 
double Props(char *Output, char Name1, double Prop1, char Name2, double Prop2, char * Ref);
double Props1(char *Output, char * Ref);
void get_errstring(char *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int npol,Nel;
    int m,n,u,v,k,Output_len,Name1_len,Name2_len,Ref_len;
    int *c,status;
    double *d,Prop1,Prop2,x,y;
	char *Output,*Name1,*Name2,*Ref,errstr[1000];
    mxArray *cMat[1];
    
    /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    
    if (nrhs == 2 && mxIsChar (prhs[0]) && mxIsChar (prhs[1]))
    {
        // Get the refrigerant (it is a string) (+1 for the NULL terminator)
        Ref_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Ref = mxCalloc(Ref_len, sizeof(char));
        status = mxGetString(prhs[0], Ref, Ref_len);
        
        // Get the output (it is a string) (+1 for the NULL terminator)
        Output_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Output = mxCalloc(Output_len, sizeof(char));
        status = mxGetString(prhs[1], Output, Output_len);
        
        *mxGetPr(plhs[0]) = Props1(Ref,Output);
    }
    else
    {
        // Get the refrigerant (it is a string) (+1 for the NULL terminator)
        Ref_len=(mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
        Ref = mxCalloc(Ref_len, sizeof(char));
        
        // Get the output (it is a string) (+1 for the NULL terminator)
        Output_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Output = mxCalloc(Output_len, sizeof(char));
        
        status = mxGetString(prhs[5], Ref, Ref_len);
        status = mxGetString(prhs[0], Output, Output_len);
        
        Name1_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Name1 = mxCalloc(Name1_len, sizeof(char));
        mxGetString(prhs[1], Name1, Name1_len);
        
        Name2_len=(mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
        Name2 = mxCalloc(Name2_len, sizeof(char));
        mxGetString(prhs[3], Name2, Name2_len);
        
        Prop1 = mxGetScalar(prhs[2]);
        Prop2 = mxGetScalar(prhs[4]);

        x = Props(Output,Name1[0],Prop1,Name2[0],Prop2,Ref);
        *mxGetPr(plhs[0])=x;
    }
}
