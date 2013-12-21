/*
 * PseudoPureFluid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef PSEUDOPUREFLUID_H_
#define PSEUDOPUREFLUID_H_

#include "CoolPropFluid.h"

namespace CoolProp {

class PseudoPureFluid: public CoolProp::CoolPropFluid {
public:
	PseudoPureFluid();
	virtual ~PseudoPureFluid();
};

} /* namespace CoolProp */
#endif /* PSEUDOPUREFLUID_H_ */
