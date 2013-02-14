#ifndef COOLPROPSOLVER_H_
#define COOLPROPSOLVER_H_

#include "basesolver.h"

//! CoolProp solver class
/*!
  This class defines a solver that calls out to the open-source CoolProp property database.

  libraryName = "CoolProp";

  Ian Bell (ian.h.bell@gmail.com)
  2012-2013
  University of Liege, Liege, Belgium
*/
class CoolPropSolver : public BaseSolver{
protected:
	class CoolPropStateClass *state;
	bool enable_TTSE, calc_transport;
	int debug_level;
public:
	CoolPropSolver(const string &mediumName, const string &libraryName, const string &substanceName);
	~CoolPropSolver(){};
	virtual void setFluidConstants();

	virtual void setSat_p(double &p, ExternalSaturationProperties *const properties);
	virtual void setSat_T(double &T, ExternalSaturationProperties *const properties);

	virtual void setState_ph(double &p, double &h, int &phase, ExternalThermodynamicState *const properties);
	virtual void setState_pT(double &p, double &T, ExternalThermodynamicState *const properties);
	virtual void setState_dT(double &d, double &T, int &phase, ExternalThermodynamicState *const properties);
	virtual void setState_ps(double &p, double &s, int &phase, ExternalThermodynamicState *const properties);
};

#endif // COOLPROPSOLVER_H_
