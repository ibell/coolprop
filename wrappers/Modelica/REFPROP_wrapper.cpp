/*
    wrapper code to use CoolProp to calculate fluid properties.  
    
    CoolProp2Modelica is based on REFPROP2Modelica written by
    Henning Francke
	francke@gfz-potsdam.de
    Helmholtz Centre Potsdam
	GFZ German Research Centre for Geosciences
	Telegrafenberg, D-14473 Potsdam
    
    This file is released under the Modelica License 2.
*/
//#define DEBUGMODE 1

#include <windows.h>
#include <stdio.h>
#include "REFPROP_wrapper.h"
#include "CoolProp.h"

double props_REFPROP(char* what, char* statevars_in, char* fluidnames, double *props, double statevar1, double statevar2, double* x, int phase, char* REFPROP_PATH, char* errormsg, int DEBUGMODE){
/*Calculates thermodynamic properties of a pure substance/mixture, returns both single value and array containing all calculated values (because the are calculated anyway)
INPUT: 
	what: character specifying return value (p,T,h,s,d,wm,q,e,w) - Explanation of variables at the end of this function
	statevars: string of any combination of two variables out of p,T,h,s,d
	fluidnames: string containing names of substances in mixtured separated by |, substance names are identical to those of *.fld-files in REFPROP program directory
	statevar1,statevar2: values of the two variables specified in statevars
 	x: array containing the mass fractions of the components of the mixture
 	REFPROP_PATH: string defining the path of the refprop.dll
OUTPUT
	return value: value of variable specified by the input variable what
	props: Array containing all calculated values (props[0] containing error number)
 	errormsg: string containing error message
    
    Keeping for now the same interface as for REFPROP2Modelica, will change the names soon
*/
	double p, T, d, val, dl,dv,q,e,h,s,cv,cp,w,wm,wmliq,wmvap,eta,tcx;
    return Props(toupper(what[0]),toupper(statevars_in[0]),statevar1,toupper(statevars_in[1]),statevar2,fluidnames);
}

//---------------------------------------------------------------------------

double satprops_REFPROP(char* what, char* statevar, char* fluidnames, double *props, double statevarval, double* x, char* REFPROP_PATH, char* errormsg, int DEBUGMODE){
/*Calculates thermodynamic saturation properties of a pure substance/mixture, returns both single value and array containing all calculated values (because the are calculated anyway)
INPUT: 
	what: character specifying return value (p,T,h,s,d,wm,q,e,w) - Explanation of variables at the end of this function
	statevar: string of 1 variable out of p,T,h,s,d
	fluidnames: string containing names of substances in mixtured separated by |, substance names are identical to those of *.fld-files in REFPROP program directory
	statevarval: values of the variable specified in statevar
 	x: array containing the mass fractions of the components of the mixture
 	REFPROP_PATH: string defining the path of the refprop.dll
OUTPUT
	return value: value of variable specified by the input variable what
	props: Array containing all calculated values
 	errormsg: string containing error message
*/
	double p, T, d, val, dl,dv,wm,wmliq,wmvap;
    // This is a bogus function call
	return Props(toupper(what[0]),toupper(statevar[0]),statevarval,toupper(statevar[1]),statevarval,fluidnames);
}