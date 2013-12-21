/*
 * RefpropFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPFLUID_H_
#define REFPROPFLUID_H_

#include "AbstractFluid.h"

namespace CoolProp {

class RefpropFluid: public CoolProp::AbstractFluid {
public:
	RefpropFluid();
	virtual ~RefpropFluid();
};

} /* namespace CoolProp */
#endif /* REFPROPFLUID_H_ */
