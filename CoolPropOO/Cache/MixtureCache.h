/*
 * MixtureCache.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef MIXTURECACHE_H_
#define MIXTURECACHE_H_

#include "AbstractCache.h"

namespace CoolProp {

class MixtureCache: public CoolProp::AbstractCache {
public:
	MixtureCache();
	virtual ~MixtureCache();
};

} /* namespace CoolProp */
#endif /* MIXTURECACHE_H_ */
