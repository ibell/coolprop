#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <crtdbg.h>
#endif

#include "CoolProp.h"
#include "CPState.h"
#include "FluidClass.h"
#include "HumidAirProp.h"

EXPORT_CODE void CONVENTION HAPropsScilab(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char * Name3, double * Prop3, double *output)
{
	*output = HAProps(Output,Name1,*Prop1,Name2,*Prop2,Name3,*Prop3);
}

EXPORT_CODE void CONVENTION PropsSIScilab(char *Output, char *Name1, double *Prop1, char *Name2, double *Prop2, char * Ref, double *output)
{
	*output = PropsSI(Output,Name1,*Prop1,Name2,*Prop2,Ref);
}

