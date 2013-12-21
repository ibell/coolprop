/*
 * FluidCache.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef FLUIDCACHE_H_
#define FLUIDCACHE_H_

#include "AbstractCache.h"

namespace CoolProp {

class FluidCache: public CoolProp::AbstractCache {
public:
	FluidCache();
	virtual ~FluidCache();
};

} /* namespace CoolProp */
#endif /* FLUIDCACHE_H_ */
