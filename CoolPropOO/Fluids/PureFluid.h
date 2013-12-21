/*
 * PureFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef PUREFLUID_H_
#define PUREFLUID_H_

#include "CoolPropFluid.h"

namespace CoolProp {

class PureFluid: public CoolProp::CoolPropFluid {
public:
	PureFluid();
	virtual ~PureFluid();
};

} /* namespace CoolProp */
#endif /* PUREFLUID_H_ */
