#include "CoolProp.h"
#include "CPState.h"
#include "TTSE.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "time.h"

TTSESinglePhaseTable::TTSESinglePhaseTable(Fluid *pFluid, int Nrow, int Ncol)
{
	this->Nrow = Nrow;
	this->Ncol = Ncol;
	this->pFluid = pFluid;

	h.resize(Nrow);
	p.resize(Ncol);

	// Seed the generator
	srand((unsigned int)time(NULL));

	s.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	dsdh.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	dsdp.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2sdh2.resize(Nrow, std::vector<double>(Ncol, _HUGE));	
	d2sdp2.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2sdhdp.resize(Nrow, std::vector<double>(Ncol, _HUGE));

	T.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	dTdh.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	dTdp.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2Tdh2.resize(Nrow, std::vector<double>(Ncol, _HUGE));	
	d2Tdp2.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2Tdhdp.resize(Nrow, std::vector<double>(Ncol, _HUGE));

	rho.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	drhodh.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	drhodp.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2rhodh2.resize(Nrow, std::vector<double>(Ncol, _HUGE));	
	d2rhodp2.resize(Nrow, std::vector<double>(Ncol, _HUGE));
	d2rhodhdp.resize(Nrow, std::vector<double>(Ncol, _HUGE));

	// make the CoolPropStateClass instance that will be used throughout
	pCPS = new CoolPropStateClass(pFluid);

	//// Enthalpy at saturated liquid at triple point temperature
	//double hmin = Props("H",'T',pFluid->params.Ttriple+1e-8,'Q',0,pFluid->get_name());
	//double hmax = Props("H",'T',pFluid->params.Ttriple+1e-8,'Q',1,pFluid->get_name())*2;
	//double pmin = pFluid->params.ptriple;
	//double pmax = 2*pFluid->reduce.p;
}

double TTSESinglePhaseTable::build(double hmin, double hmax, double pmin, double pmax, bool logp)
{
	this->hmin = hmin;
	this->hmax = hmax;
	this->pmin = pmin;
	this->pmax = pmax;

	CoolPropStateClass CPS = *pCPS;

	dh = (hmax - hmin)/(Nrow - 1);
	dp = (pmax - pmin)/(Ncol - 1);

	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i<Nrow; i++)
	{
		double hval = hmin + i*dh;
		h[i] = hval;
		for (unsigned int j = 0; j<Ncol; j++)
		{
			double pval = pmin + j*dp;
			p[j] = pval;

			CPS.update(iH, hval, iP, pval);
			double T = CPS.T();
			double rho = CPS.rho();
			double cp = CPS.cp();

			// Matrices for entropy as a function of pressure and enthalpy
			s[i][j] = CPS.s();
			dsdh[i][j] = 1/T;
			dsdp[i][j] = -1/(T*rho);
			d2sdh2[i][j] = -1/(T*T)/CPS.dhdT_constp();
			d2sdhdp[i][j] = -1/(T*T)/CPS.dpdT_consth();
			d2sdp2[i][j] = 1/(T*T*rho)/CPS.dhdT_constp()+1/(T*rho*rho)/CPS.dpdrho_consth();

			// These are common terms needed for a range of terms for T(h,p) as well as rho(h,p)
			double A = CPS.dpdT_constrho()*CPS.dhdrho_constT()-CPS.dpdrho_constT()*CPS.dhdT_constrho();
			double dAdT_constrho = CPS.d2pdT2_constrho()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrhodT()-CPS.d2pdrhodT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdT2_constrho();
			double dAdrho_constT = CPS.d2pdrhodT()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrho2_constT()-CPS.d2pdrho2_constT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdrhodT();

			//Matrices for temperature as a function of pressure and enthalpy
			double ddT_dTdp_h_constrho = 1/A*CPS.d2hdrhodT()-1/(A*A)*dAdT_constrho*CPS.dhdrho_constT();
			double ddrho_dTdp_h_constT = 1/A*CPS.d2hdrho2_constT()-1/(A*A)*dAdrho_constT*CPS.dhdT_constrho();
			double ddT_dTdh_p_constrho = -1/(cp*cp)*(CPS.d2hdT2_constrho()-CPS.dhdp_constT()*CPS.d2pdT2_constrho()+CPS.d2hdrhodT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrhodT());
			double ddrho_dTdh_p_constT = -1/(cp*cp)*(CPS.d2hdrhodT()-CPS.dhdp_constT()*CPS.d2pdrhodT()+CPS.d2hdrho2_constT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrho2_constT());
			this->T[i][j] = T;
			dTdh[i][j] = 1/cp;
			dTdp[i][j] = 1/A*CPS.dhdrho_constT();
			d2Tdh2[i][j]  = ddT_dTdh_p_constrho/CPS.dhdT_constp()+ddrho_dTdh_p_constT/CPS.dhdrho_constp();
			d2Tdhdp[i][j] = ddT_dTdp_h_constrho/CPS.dhdT_constp()+ddrho_dTdp_h_constT/CPS.dhdrho_constp();
			d2Tdp2[i][j]  = ddT_dTdp_h_constrho/CPS.dpdT_consth()+ddrho_dTdp_h_constT/CPS.dpdrho_consth();

			//// Matrices for density as a function of pressure and enthalpy
			double ddT_drhodp_h_constrho = -1/A*CPS.d2hdT2_constrho()+1/(A*A)*dAdT_constrho*CPS.dhdT_constrho();
			double ddrho_drhodp_h_constT = -1/A*CPS.d2hdrhodT()+1/(A*A)*dAdrho_constT*CPS.dhdrho_constT();
			double ddT_drhodh_p_constrho = 1/A*CPS.d2pdT2_constrho()-1/(A*A)*dAdT_constrho*CPS.dpdT_constrho();
			double ddrho_drhodh_p_constT = 1/A*CPS.d2pdrhodT()-1/(A*A)*dAdrho_constT*CPS.dpdrho_constT();
			this->rho[i][j] = rho;
			drhodh[i][j] = 1/A*CPS.dpdT_constrho();
			drhodp[i][j] = -1/A*CPS.dhdT_constrho();
			d2rhodh2[i][j]  = ddT_drhodh_p_constrho/CPS.dhdT_constp()+ddrho_drhodh_p_constT/CPS.dhdrho_constp();
			d2rhodhdp[i][j] = ddT_drhodp_h_constrho/CPS.dhdT_constp()+ddrho_drhodp_h_constT/CPS.dhdrho_constp();
			d2rhodp2[i][j]  = ddT_drhodp_h_constrho/CPS.dpdT_consth()+ddrho_drhodp_h_constT/CPS.dpdrho_consth();
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	return elap;
}
double TTSESinglePhaseTable::check_randomly(long iParam, unsigned int N, std::vector<double> *h, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE)
{	
	double val=0;
	h->resize(N);
	p->resize(N);
	EOS->resize(N);
	TTSE->resize(N);
	
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;
		
		CoolPropStateClass CPS = *pCPS;
		CPS.update(iH,h1,iP,p1);
		double sEOS = CPS.s();
		double cpEOS = CPS.cp();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*h)[i] = h1;
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,h1,p1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOS)[i] = sEOS; break;
		case iT:
			(*EOS)[i] = TEOS; break;
		case iC:
			(*EOS)[i] = cpEOS; break;
		case iD:
			(*EOS)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g TTSE (h,p,EOS,TTSE)\n",h1,p1,(*EOS)[i],(*TTSE)[i]);
	}
	return val;
}

double TTSESinglePhaseTable::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;

		// Get the value from TTSE
		evaluate(iParam,h1,p1);
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}


double TTSESinglePhaseTable::evaluate(long iParam, double h, double p)
{
	int i = (int)((h-hmin)/(hmax-hmin)*(Nrow-1));
	int j = (int)((p-pmin)/(pmax-pmin)*(Ncol-1));
	double deltah = h-this->h[i];
	double deltap = p-this->p[j];
	
	switch (iParam)
	{
	case iS:
		return s[i][j]+deltah*dsdh[i][j]+deltap*dsdp[i][j]+0.5*deltah*deltah*d2sdh2[i][j]+0.5*deltap*deltap*d2sdp2[i][j]+deltap*deltah*d2sdhdp[i][j]; break;
	case iT:
		return T[i][j]+deltah*dTdh[i][j]+deltap*dTdp[i][j]+0.5*deltah*deltah*d2Tdh2[i][j]+0.5*deltap*deltap*d2Tdp2[i][j]+deltap*deltah*d2Tdhdp[i][j]; break;
	case iD:
		return rho[i][j]+deltah*drhodh[i][j]+deltap*drhodp[i][j]+0.5*deltah*deltah*d2rhodh2[i][j]+0.5*deltap*deltap*d2rhodp2[i][j]+deltap*deltah*d2rhodhdp[i][j]; break;
	default:
		throw ValueError();
	}
	return 0;
}

TTSETwoPhaseTable::TTSETwoPhaseTable(Fluid *pFluid, int N, double Q)
{

	this->N = N;
	this->pFluid = pFluid;
	this->Q = Q;

	// Seed the generator
	srand((unsigned int)time(NULL));

	// Resize all the arrays
	h.resize(N);
	p.resize(N);
	T.resize(N);
	dTdp.resize(N);
	d2Tdp2.resize(N);
	rho.resize(N);
	drhodp.resize(N);
	d2rhodp2.resize(N);
	s.resize(N);
	dsdp.resize(N);
	d2sdp2.resize(N);
	h.resize(N);
	dhdp.resize(N);
	d2hdp2.resize(N);

	// make the CoolPropStateClass instance that will be used throughout
	pCPS = new CoolPropStateClass(pFluid);
}

double TTSETwoPhaseTable::build(double pmin, double pmax)
{
	CoolPropStateClass CPS = *pCPS;

	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);
	this->logpmax = log(pmax);

	double dlogp = (logpmax-logpmin)/(N-1);
	clock_t t1,t2;
	t1 = clock();
	// Linear distribution of pressures
	for (unsigned int i = 0; i < N; i++)
	{
		p[i] = exp(logpmin + i*dlogp);
		CPS.update(iP,p[i],iQ,Q);
		T[i] = CPS.T();
		dTdp[i] = CPS.dTdp_along_sat();
		d2Tdp2[i] = CPS.d2Tdp2_along_sat();
		h[i] = CPS.h();
		dhdp[i] = CPS.dhdp_along_sat_vapor();
		d2hdp2[i] = CPS.d2hdp2_along_sat_vapor();
		s[i] = CPS.s();
		dsdp[i] = CPS.dsdp_along_sat_vapor();
		d2sdp2[i] = CPS.d2sdp2_along_sat_vapor();
		rho[i] = CPS.rho();
		drhodp[i] = CPS.drhodp_along_sat_vapor();
		d2rhodp2[i] = CPS.d2rhodp2_along_sat_vapor();
	}
	t2 = clock();
	return double(t2-t1)/CLOCKS_PER_SEC;
}

double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

double TTSETwoPhaseTable::evaluate(long iParam, double p)
{
	CoolPropStateClass CPS = *pCPS;

	int i = (int)round(((log(p)-logpmin)/(logpmax-logpmin)*(N-1)));
	double deltap = log(p/this->p[i]);
	double log_PI_PIi = log(p/this->p[i]);
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iS:
		return s[i]+log_PI_PIi*pi*dsdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2sdp2[i]*pi*pi;
	case iT:
		return T[i]+log_PI_PIi*pi*dTdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2Tdp2[i]*pi*pi;
	case iH:
		return h[i]+log_PI_PIi*pi*dhdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2hdp2[i]*pi*pi;
	case iD:
		// log(p) v. log(rho) gives close to a line for most of the curve
		return exp(log(rho[i])+log_PI_PIi*(1.0+0.5*log_PI_PIi*(1-pi/rho[i]*drhodp[i]))*pi/rho[i]*drhodp[i]+0.5*log_PI_PIi*log_PI_PIi*d2rhodp2[i]*pi*pi/rho[i]);
	default:
		throw ValueError();
	}
}

double TTSETwoPhaseTable::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;

		// Get the value from TTSE
		evaluate(iParam,p1);
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}


double TTSETwoPhaseTable::check_randomly(long iParam, double Q, unsigned int N, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE)
{	
	double val=0;
	p->resize(N);
	EOS->resize(N);
	TTSE->resize(N);
	
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		
		CoolPropStateClass CPS = *pCPS;
		CPS.update(iP,p1,iQ,Q);
		double hEOS = CPS.h();
		double sEOS = CPS.s();
		double cpEOS = CPS.cp();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,p1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOS)[i] = sEOS; break;
		case iT:
			(*EOS)[i] = TEOS; break;
		case iH:
			(*EOS)[i] = hEOS; break;
		case iD:
			(*EOS)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g TTSE (p,EOS,TTSE, delta [mK])\n",p1,(*EOS)[i],(*TTSE)[i],((*EOS)[i]-(*TTSE)[i]));
	}
	return val;
}

