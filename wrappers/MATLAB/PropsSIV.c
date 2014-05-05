#include "mex.h"   /*--This one is required*/

/* Prototype for the Props function to be called.  Can't use the CoolProp.h header because there are a lot of
   c++ parts in the header that cannot be easily hidden when compiling */
double PropsSI(char *Output, char *Name1, double Prop1, char *Name2, double Prop2, char * Ref);
double Props1(char *Output, char * Ref);
long get_global_param_string(char*, char*);
long get_fluid_param_string(char *fluid, char *param, char * Output);
long get_standard_unit_system(void);
void set_standard_unit_system(long);
#include "GlobalConstants.h"
#include "float.h"

bool ValidNumber(double x)
{
    // Idea from http://www.johndcook.com/IEEE_exceptions_in_cpp.html
    return (x <= DBL_MAX && x >= -DBL_MAX);
};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int npol,Nel;
    int m,n,u,v,k,Output_len,Name1_len,Name2_len,Ref_len,Param_len;
    int *c,status;
    double *d,Prop1,Prop2,x,y;
	char *Output,*Name1,*Name2,*Ref,*Param,errstr[1000],errstr2[1000],fluidslist[10000];
    double val;
    double *rval;
    mxArray *cMat[1];
    
    double *pro1,*pro2;
    size_t mrows1,ncols1,mrows2,ncols2;
    size_t i;
    
    if (nrhs == 6 && mxIsChar (prhs[0]) && mxIsChar (prhs[1])  && mxIsChar (prhs[3])  && mxIsChar (prhs[5])  )
    {
        
        pro1 = mxGetPr(prhs[2]);mrows1 = mxGetM(prhs[2]);ncols1 = mxGetN(prhs[2]);
        pro2 = mxGetPr(prhs[4]);mrows2 = mxGetM(prhs[4]);ncols2 = mxGetN(prhs[4]);        
       
        if(!((mrows1==mrows2) && (ncols1==ncols2))){
            mexErrMsgIdAndTxt( "MATLAB:PROSV",
            "Val uses of not same size."); 
        }
        if(!(ncols1==1 && ncols2==1     )){
            mexErrMsgIdAndTxt( "MATLAB:PROSV",
                               "Not implemented for matrix input");
        }
        /* Create matrix for the return argument. */
        plhs[0] = mxCreateDoubleMatrix(mrows1,ncols1, mxREAL);
        rval = mxGetPr(plhs[0]);
        
        /* Get the refrigerant (it is a string) (+1 for the NULL terminator)*/
        Ref_len=(mxGetM(prhs[5]) * mxGetN(prhs[5])) + 1;
        Ref = mxCalloc(Ref_len, sizeof(char));
        mxGetString(prhs[5], Ref, Ref_len);
        
        /* Get the output (it is a string) (+1 for the NULL terminator)*/
        Output_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Output = mxCalloc(Output_len, sizeof(char));
        mxGetString(prhs[0], Output, Output_len);
        
        Name1_len=(mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
        Name1 = mxCalloc(Name1_len, sizeof(char));
        mxGetString(prhs[1], Name1, Name1_len);
        
        Name2_len=(mxGetM(prhs[3]) * mxGetN(prhs[3])) + 1;
        Name2 = mxCalloc(Name2_len, sizeof(char));
        mxGetString(prhs[3], Name2, Name2_len);

        /*
        Prop1 = mxGetScalar(prhs[2]);
        Prop2 = mxGetScalar(prhs[4]);
        */
       
        for(i=0;i<mrows1;++i){
            Prop1 = pro1[i];
            Prop2 = pro2[i];
            val = PropsSI(Output,Name1,Prop1,Name2,Prop2,Ref);
            if (ValidNumber(val))
            {
                rval[i]=val;
                /* *mxGetPr(plhs[0]) = val;*/
            } else{
            get_global_param_string("errstring",errstr);
            sprintf(errstr2,"CoolProp Error: %s",errstr);
            mexErrMsgTxt(errstr2);
            }

        }
        /* If it is a good value, return it */
        /* Otherwise there was an error, return the CoolProp error*/
       
    }
    else if (nrhs == 1 && mxIsChar (prhs[0]))
    {
        Ref_len=(mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
        Ref = mxCalloc(Ref_len, sizeof(char));
        status = mxGetString(prhs[0], Ref, Ref_len);
        
        if (!strcmp(Ref,"FluidsList")){
            get_global_param_string("FluidsList",fluidslist);
        }
        else if (!strcmp(Ref,"version")){
            get_global_param_string("version",fluidslist);
        }
        else if (!strcmp(Ref,"gitrevision")){
            get_global_param_string("gitrevision", fluidslist);
        }
        else
        {
            sprintf(errstr2,"single input is invalid: %s",Ref);
            mexErrMsgTxt(errstr2);
        }
        plhs[0] = mxCreateString(fluidslist);
        return;
    }
    else
    {
        mexErrMsgTxt("Props must either receive two strings or the signature Props(Output,Param1,Value1,Param2,Value2,FluidName)");
    }
}
