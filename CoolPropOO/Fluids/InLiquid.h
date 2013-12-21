/*
 * InLiquid.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INLIQUID_H_
#define INLIQUID_H_

#include "Incompressible.h"

namespace CoolProp {

class InLiquid: public CoolProp::Incompressible {
public:
	InLiquid();
	virtual ~InLiquid();
};

} /* namespace CoolProp */
#endif /* INLIQUID_H_ */
