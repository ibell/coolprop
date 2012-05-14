
#ifndef TYPE_TP
#define TYPE_TP 1
#endif

#ifndef TYPE_Trho
#define TYPE_Trho 2
#endif

#ifndef TYPE_TPNoLookup
#define TYPE_TPNoLookup 3
#endif

#ifndef ERRSTRLENGTH
#define ERRSTRLENGTH 255
#endif

#include <math.h>

//MINGW version of huge value macro
#ifdef HUGE_VAL
#define _HUGE HUGE_VAL
#else
	// GCC Version of huge value macro
	#ifdef HUGE 
	#define _HUGE HUGE
	#endif
#endif

