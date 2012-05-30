#include "FluidClass.h"
#include <complex>
#include "PengRobinson.h"

std::vector<double> PRGuess_rho(Fluid * pFluid, double T, double p)
{
	// See http://ascend4.org/PengRobinson_EOS_in_FPROPS or 
	std::complex<double> cQ,cdelta,cC,cZ1,cZ2,cZ3,cplussqrt3,cminussqrt3;
	std::vector<double> solns;
	double alpha,kappa,omega,Tc,pc,a,b,c,d,A,B,DELTA,R,delta,Z,Q,Z1,Z2,Z3,sqrt3,rho;
	omega = pFluid->params.accentricfactor;
	Tc = pFluid->reduce.T;
	pc = pFluid->reduce.p;
	R = pFluid->params.R_u/pFluid->params.molemass;

	kappa = 0.37464 + 1.54225 * omega - 0.26992 * omega*omega;
	alpha = pow(1+kappa*(1-sqrt(T/Tc)),2);
	a = 0.45724 * R*R * Tc*Tc/pc; //[R*T]
	b = 0.077796*R*Tc/pc; //[m^3/kg]

	A = a*alpha*p/pow(R*T,2);
	B = b*p/(R*T);
	
	// Peng-Robinson EOS is of the form
	// Z^3 + (B-1)*Z^2 +(A-3B^2-2B)*Z-AB+B^2+B^3=0
	// coeffs are
	a = 1;
	b = B-1;
	c = (A - 3*B*B - 2*B);
	d = -A*B + B*B + B*B*B;
	DELTA = 18*a*b*c*d-4*pow(b,3)*d+pow(b,2)*pow(c,2)-4*a*pow(c,3)-27*pow(a,2)*pow(d,2);
	delta = 2*pow(b,3)-9*a*b*c+27*pow(a,2)*d;

	Q = sqrt(-27*pow(a,2)*DELTA);
	cQ = sqrt(-27*pow(a,2)*std::complex<double>(DELTA,0));
	if (Q==0 && (pow(b,2)-3*a*c)!=0)
	{
		// Three roots; one multiple, and a pair
		Z1 = (b*c-9*a*d)/(2*(3*a*c-pow(b,2)));
		Z2 = Z1;
		Z3 = (9*pow(a,2)*d - 4*a*b*c + pow(b,3))/(a*(3*a*c-pow(b,2)));
		// If both roots are positive, use the one closer to one
		if (Z1>0 && Z3>0)
		{
			if (fabs(1-Z1)<fabs(1-Z3))
			{
				Z = Z1;
			}
			else
			{
				Z = Z3;
			}
		}
		else if (Z1>0 && Z3<0)
		{
			Z=Z1;
		}
		else
		{
			Z=Z3;
		}
		rho = p/(R*T*Z);
		solns.push_back(rho);
	}
	else if (Q==0 && (pow(b,2)-3*a*c)==0)
	{
		// All three roots are equal and real
		Z = -b/(3*a);
		rho = p/(R*T*Z);
		solns.push_back(rho);
	}
	else if (pow(b,2)-3*a*c!=0)
	{ 
		// Some unpleasant mixture of complex roots
		cdelta = std::complex<double>(delta,0);
		cC = pow(0.5*(cdelta-cQ),1.0/3.0);
		cZ1 = -1/(3.0*a)*(b+cC+(pow(b,2)-3*a*c)/cC);
		sqrt3 = sqrt((double)3);
		cplussqrt3 = std::complex<double>(1,sqrt3);
		cminussqrt3 = std::complex<double>(1,-sqrt3);

		cZ2 = -1.0/(3.0*a)*(b-cC*cplussqrt3/2.0-cminussqrt3*(pow(b,2)-3.0*a*c)/(2.0*cC));
		cZ3 = -1.0/(3.0*a)*(b-cC*cminussqrt3/2.0-cplussqrt3*(pow(b,2)-3.0*a*c)/(2.0*cC));

		// Find the densities that are real and positive
		if (fabs(cZ1.imag())<1e-14 && cZ1.real()>0)
			solns.push_back(p/(cZ1.real()*R*T));
		if (fabs(cZ2.imag())<1e-14 && cZ2.real()>0)
			solns.push_back(p/(cZ2.real()*R*T));
		if (fabs(cZ3.imag())<1e-14 && cZ3.real()>0)
			solns.push_back(p/(cZ3.real()*R*T));
	}
	else if (pow(b,2)-3*a*c==0)
	{
		std::cout << "Somehow b^2-3*a*c in PRGuess_Delta is equal to zero, using ideal gas" << std::endl;
		rho=p/(R*T);
	}
	
	return solns;
}