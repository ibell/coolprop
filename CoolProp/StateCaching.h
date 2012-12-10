#include <list>
#include <string>
#include "CPState.h"

#ifndef STATE_CACHING_H
#define STATE_CACHING_H

class StateCache
{
private:
	std::list<CoolPropStateClass*> state_cache;
	long Nmax;
public:
	StateCache(){};
	~StateCache();

	bool add(CoolPropStateClass* CPS);
	bool use(long iInput1, double Prop1, long iInput2, double Prop2, double *T0, double *rho0);
};

#endif
