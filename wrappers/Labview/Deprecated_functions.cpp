#include "CoolPropTools.h"
#include "CoolProp.h"
EXPORT_CODE long CONVENTION get_errstring_copy(char * Output)
{
	strcpy(Output,get_global_param_string("errstring").c_str());
	return 0;
}