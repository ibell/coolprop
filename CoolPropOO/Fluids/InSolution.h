/*
 * InSolution.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef INSOLUTION_H_
#define INSOLUTION_H_

#include "Incompressible.h"

namespace CoolProp {

class InSolution: public CoolProp::Incompressible {
public:
	InSolution();
	virtual ~InSolution();
};

} /* namespace CoolProp */
#endif /* INSOLUTION_H_ */
