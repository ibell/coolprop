/*
 * RefpropMixture.h
 *
 *  Created on: 23 Dec 2013
 *      Author: jowr
 */

#ifndef REFPROPMIXTURE_H_
#define REFPROPMIXTURE_H_

#include "AbstractMixture.h"

namespace CoolProp {

class RefpropMixture: public CoolProp::AbstractMixture {
public:
	RefpropMixture();
	virtual ~RefpropMixture();
};

} /* namespace CoolProp */
#endif /* REFPROPMIXTURE_H_ */
