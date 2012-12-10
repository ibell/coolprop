#include <vector>
#include "Solvers.h"
#include "math.h"
#include <Eigen/Dense>
#include <iostream>
#include "CPExceptions.h"
#include "CoolPropTools.h"
#include "StateCaching.h"

StateCache::~StateCache()
{
	while (!state_cache.empty())
	{
		delete state_cache.back();
		state_cache.pop_back();
	}
}

bool StateCache::add(CoolPropStateClass *CPS)
{
	if (state_cache.size()== Nmax)
	{
		state_cache.pop_front();
	}
	state_cache.push_back(CPS); 
	return true;
}

std::map<double, CoolPropStateClass*> distance_map;

bool StateCache::use(long iInput1, double Value1, long iInput2, double Value2, double *T, double *rho)
{
	double dist, dist_min = _HUGE;
	CoolPropStateClass *CPS;

	if (match_pair(iInput1,iInput2,iP,iH))
	{
		// Sort so they are in the order P, H
		sort_pair(&iInput1,&Value1,&iInput2,&Value2,iP,iH);
		// Create a map from distance to element in cache - map is sorted by key, so elements will be sorted by distance
		// Have to use a multimap because you might have equal distances for multiple states and don't want to throw out that information!
		for (std::list<CoolPropStateClass*>::iterator it = state_cache.begin(); it != state_cache.end(); it++)
		{
			dist = sqrt(pow((*it)->p()-Value1,2)+pow((*it)->h()-Value2,2));
			if (dist < dist_min)
			{
				CPS = (*it);
			}
		}
		*T = CPS->T();
		*rho = CPS->rho();
		return true;
	}
	else
	{
		return false;
	}
}