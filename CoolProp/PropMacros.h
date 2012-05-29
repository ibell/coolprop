#include <math.h>
#ifdef HUGE_VAL
#define _HUGE HUGE_VAL
#else
	// GCC Version of huge value macro
	#ifdef HUGE 
	#define _HUGE HUGE
	#endif
#endif

