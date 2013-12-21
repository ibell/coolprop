/*
 * CoolPropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef COOLPROPFLUID_H_
#define COOLPROPFLUID_H_

#include "AbstractFluid.h"

namespace CoolProp {

class CoolPropFluid: public CoolProp::AbstractFluid {
public:
	CoolPropFluid();
	virtual ~CoolPropFluid();
};

} /* namespace CoolProp */
#endif /* COOLPROPFLUID_H_ */
