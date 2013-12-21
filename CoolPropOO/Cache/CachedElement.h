/*
 * CachedElement.h
 *
 *  Created on: 21 Dec 2013
 *      Author: jowr
 */

#ifndef CACHEDELEMENT_H_
#define CACHEDELEMENT_H_

#include "../DataStructures.h"

namespace CoolProp {

/*!
A class that contains the magic to cache a value.

Includes an "=" assignment operator and casting to boolean
so you can do something like::

double CoolPropStateClassSI::d3phir_dTau3(double tau, double delta){
	if (cache.d3phir_dTau3)	{
		return cache.d3phir_dTau3;
	} else {
		cache.d3phir_dTau3 = pFluid->d3phir_dTau3(tau,delta);
		return cache.d3phir_dTau3;
	}
};
*/

class CachedElement {
private:
	bool is_cached;
	double value;
public:
	/// Default constructor
	CachedElement() {
		this->clear();
	};
	/// Assignment operator - sets the value and sets the flag
	void operator=(const double& value) {
		this->value = value;
		this->is_cached = true;
	};
	/// Cast to boolean, for checking if cached
	operator bool() const {return is_cached;};
	/// Cast to double, for returning value
	operator double() const {return value;};
	/// Clear the flag and the value
	void clear() {
		is_cached = false;
		this->value = _HUGE;};
};

} /* namespace CoolProp */
#endif /* CACHEDELEMENT_H_ */
