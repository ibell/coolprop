/*
 * CoolPropMixture.h
 *
 *  Created on: 23 Dec 2013
 *      Author: jowr
 */

#ifndef COOLPROPMIXTURE_H_
#define COOLPROPMIXTURE_H_

#include "AbstractMixture.h"

namespace CoolProp {

class CoolPropMixture: public CoolProp::AbstractMixture {
public:
	CoolPropMixture();
	virtual ~CoolPropMixture();
};

} /* namespace CoolProp */
#endif /* COOLPROPMIXTURE_H_ */
