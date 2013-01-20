#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "CoolProp.h"
#include "CPState.h"
#include "TTSE.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "time.h"

double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

TTSESinglePhaseTableClass::TTSESinglePhaseTableClass(Fluid *pFluid)
{
	this->pFluid = pFluid;

	// Seed the generator for random number generation
	srand((unsigned int)time(NULL));
}
void TTSESinglePhaseTableClass::set_size(unsigned int Nrow, unsigned int Ncol)
{
	this->Nrow = Nrow;
	this->Ncol = Ncol;

	h.resize(Nrow);
	p.resize(Ncol);

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
}

//double TTSESinglePhaseTableClass::build(double hmin, double hmax, double pmin, double pmax, bool logp)
//{
//	this->hmin = hmin;
//	this->hmax = hmax;
//	this->pmin = pmin;
//	this->pmax = pmax;
//
//	CoolPropStateClass CPS = CoolPropStateClass(pFluid);
//
//	dh = (hmax - hmin)/(Nrow - 1);
//	dp = (pmax - pmin)/(Ncol - 1);
//
//	clock_t t1,t2;
//	t1 = clock();
//	for (unsigned int i = 0; i<Nrow; i++)
//	{
//		double hval = hmin + i*dh;
//		h[i] = hval;
//		for (unsigned int j = 0; j<Ncol; j++)
//		{
//			double pval = pmin + j*dp;
//			p[j] = pval;
//
//			CPS.update(iH, hval, iP, pval);
//			double T = CPS.T();
//			double rho = CPS.rho();
//			double cp = CPS.cp();
//
//			// Matrices for entropy as a function of pressure and enthalpy
//			s[i][j] = CPS.s();
//			dsdh[i][j] = 1/T;
//			dsdp[i][j] = -1/(T*rho);
//			d2sdh2[i][j] = -1/(T*T)/CPS.dhdT_constp();
//			d2sdhdp[i][j] = -1/(T*T)/CPS.dpdT_consth();
//			d2sdp2[i][j] = 1/(T*T*rho)/CPS.dhdT_constp()+1/(T*rho*rho)/CPS.dpdrho_consth();
//
//			// These are common terms needed for a range of terms for T(h,p) as well as rho(h,p)
//			double A = CPS.dpdT_constrho()*CPS.dhdrho_constT()-CPS.dpdrho_constT()*CPS.dhdT_constrho();
//			double dAdT_constrho = CPS.d2pdT2_constrho()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrhodT()-CPS.d2pdrhodT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdT2_constrho();
//			double dAdrho_constT = CPS.d2pdrhodT()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrho2_constT()-CPS.d2pdrho2_constT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdrhodT();
//
//			//Matrices for temperature as a function of pressure and enthalpy
//			double ddT_dTdp_h_constrho = 1/A*CPS.d2hdrhodT()-1/(A*A)*dAdT_constrho*CPS.dhdrho_constT();
//			double ddrho_dTdp_h_constT = 1/A*CPS.d2hdrho2_constT()-1/(A*A)*dAdrho_constT*CPS.dhdT_constrho();
//			double ddT_dTdh_p_constrho = -1/(cp*cp)*(CPS.d2hdT2_constrho()-CPS.dhdp_constT()*CPS.d2pdT2_constrho()+CPS.d2hdrhodT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrhodT());
//			double ddrho_dTdh_p_constT = -1/(cp*cp)*(CPS.d2hdrhodT()-CPS.dhdp_constT()*CPS.d2pdrhodT()+CPS.d2hdrho2_constT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrho2_constT());
//			this->T[i][j] = T;
//			dTdh[i][j] = 1/cp;
//			dTdp[i][j] = 1/A*CPS.dhdrho_constT();
//			d2Tdh2[i][j]  = ddT_dTdh_p_constrho/CPS.dhdT_constp()+ddrho_dTdh_p_constT/CPS.dhdrho_constp();
//			d2Tdhdp[i][j] = ddT_dTdp_h_constrho/CPS.dhdT_constp()+ddrho_dTdp_h_constT/CPS.dhdrho_constp();
//			d2Tdp2[i][j]  = ddT_dTdp_h_constrho/CPS.dpdT_consth()+ddrho_dTdp_h_constT/CPS.dpdrho_consth();
//
//			//// Matrices for density as a function of pressure and enthalpy
//			double ddT_drhodp_h_constrho = -1/A*CPS.d2hdT2_constrho()+1/(A*A)*dAdT_constrho*CPS.dhdT_constrho();
//			double ddrho_drhodp_h_constT = -1/A*CPS.d2hdrhodT()+1/(A*A)*dAdrho_constT*CPS.dhdrho_constT();
//			double ddT_drhodh_p_constrho = 1/A*CPS.d2pdT2_constrho()-1/(A*A)*dAdT_constrho*CPS.dpdT_constrho();
//			double ddrho_drhodh_p_constT = 1/A*CPS.d2pdrhodT()-1/(A*A)*dAdrho_constT*CPS.dpdrho_constT();
//			this->rho[i][j] = rho;
//			drhodh[i][j] = 1/A*CPS.dpdT_constrho();
//			drhodp[i][j] = -1/A*CPS.dhdT_constrho();
//			d2rhodh2[i][j]  = ddT_drhodh_p_constrho/CPS.dhdT_constp()+ddrho_drhodh_p_constT/CPS.dhdrho_constp();
//			d2rhodhdp[i][j] = ddT_drhodp_h_constrho/CPS.dhdT_constp()+ddrho_drhodp_h_constT/CPS.dhdrho_constp();
//			d2rhodp2[i][j]  = ddT_drhodp_h_constrho/CPS.dpdT_consth()+ddrho_drhodp_h_constT/CPS.dpdrho_consth();
//		}
//	}
//	t2 = clock();
//	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
//	return elap;
//}

void TTSESinglePhaseTableClass::nearest_neighbor(int i, int j, double *T0, double *rho0)
{
	// Left
	if (i>0 && ValidNumber(rho[i-1][j]) && ValidNumber(T[i-1][j])){
		*T0 = T[i-1][j];
		*rho0 = rho[i-1][j];
		return;
	}
	// Right
	else if (i<Nrow-1 && ValidNumber(rho[i+1][j]) && ValidNumber(T[i+1][j])){
		*T0 = T[i+1][j];
		*rho0 = rho[i+1][j];
		return;
	}
	// Down
	else if (j>0 && ValidNumber(rho[i][j-1]) && ValidNumber(T[i][j-1])){
		*T0 = T[i][j-1];
		*rho0 = rho[i][j-1];
		return;
	}
	// Up
	else if (j<Ncol-1 && ValidNumber(rho[i][j+1]) && ValidNumber(T[i][j+1])){
		*T0 = T[i][j+1];
		*rho0 = rho[i][j+1];
		return;
	}
	else
		*T0 = -1;
		*rho0 = -1;
		return;
	}


	
double TTSESinglePhaseTableClass::build(double hmin, double hmax, double pmin, double pmax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV)
{
	bool SinglePhase = false;

	this->hmin = hmin;
	this->hmax = hmax;
	this->pmin = pmin;
	this->pmax = pmax;

	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	

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
			
			// Check whether the point is single phase
			// If pressure between ptriple point and pcrit, might be two-phase or sigle phase, otherwise definitely single phase
			if (pval <= pFluid->reduce.p && pval >= pFluid->params.ptriple)
			{
				if (SatL == NULL || SatV == NULL){
					// Not using TTSE method, use saturation (slow...)
					CPS.update(iP,pval,iQ,0.5);
					SinglePhase = (hval < CPS.hL()  || hval > CPS.hV());
				}
				else{
					// Using the TTSE method, nice and fast
					SinglePhase = (hval < SatL->evaluate(iH,pval)  || hval > SatV->evaluate(iH,pval));
				}
			}
			else
			{
				SinglePhase = true;
			}
			
			// If enthalpy is outside the saturation region, continue and do the calculation as a function of p,h
			if (SinglePhase)
			{
				double T0=-1,rho0=-1,T,rho,rhoL,rhoV,TsatL,TsatV;

				// Find a good point around this point that is single-phase if any of its neighbors have been calculated
				nearest_neighbor(i,j,&T0,&rho0);

				// If good T,rho was returned, use it as a guess value to calculate T,rho from p,h more quickly
				if (T0 > 0 && rho0 > 0)
				{
					// Get T,rho from p,h using our guess value
					CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
					CPS.flag_SinglePhase = true;
					CPS.update(iT, T, iD, rho);
				}
				else
				{
					// Probably here you are close to the saturation boundary since TTSE says it is single phase, but no neighbors to lean on
					CPS.flag_SinglePhase = true;

					double hsatLTTSE = SatL->evaluate(iH,pval);
					double hsatVTTSE = SatV->evaluate(iH,pval);

					if (fabs(hval-hsatLTTSE)<10)
					{
						// Close to the saturated liquid boundary
						T0 = SatL->evaluate(iT,pval);
						rho0 = SatL->evaluate(iD,pval);
						CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
						CPS.flag_SinglePhase = true;
						CPS.update(iT, T, iD, rho);
					}
					else if (fabs(hval-hsatVTTSE)<10)
					{
						// Close to the saturated vapor boundary
						T0 = SatV->evaluate(iT,pval);
						rho0 = SatV->evaluate(iD,pval);
						CPS.pFluid->temperature_ph(pval,hval,&T,&rho,&rhoL,&rhoV,&TsatL,&TsatV,T0,rho0);
						CPS.flag_SinglePhase = true;
						CPS.update(iT, T, iD, rho);
					}
					else
					{
						CPS.update(iP, pval, iH, hval);
					}
					//std::cout << format("%d %d %g %g\n",i,j,pval,hval);
				}
				
				//// TEMP
				//this->T[i][j] = T;
				//this->rho[i][j] = rho;
				//continue;
				
				T = CPS.T();
				rho = CPS.rho();
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
			else
			{
				s[i][j] = _HUGE;
				dsdh[i][j] = _HUGE;
				dsdp[i][j] = _HUGE;
				d2sdh2[i][j] = _HUGE;
				d2sdhdp[i][j] = _HUGE;
				d2sdp2[i][j] = _HUGE;

				this->T[i][j] = _HUGE;
				dTdh[i][j] = _HUGE;
				dTdp[i][j] = _HUGE;
				d2Tdh2[i][j]  = _HUGE;
				d2Tdhdp[i][j] = _HUGE;
				d2Tdp2[i][j]  = _HUGE;

				this->rho[i][j] = _HUGE;
				drhodh[i][j] = _HUGE;
				drhodp[i][j] = _HUGE;
				d2rhodh2[i][j]  = _HUGE;
				d2rhodhdp[i][j] = _HUGE;
				d2rhodp2[i][j]  = _HUGE;
			}
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout << elap << " to build single phase table" << std::endl;
	return elap;
}

void TTSESinglePhaseTableClass::write_dotdrawing_tofile(char fName[])
{
	FILE *fp;
	fp = fopen(fName,"w");
	for (int j = Ncol-1; j>=0; j--)
	{
		for (unsigned int i = 0; i<Nrow; i++)
		{
			if (ValidNumber(rho[i][j]))
			{
				fprintf(fp,".");
			}
			else
			{
				fprintf(fp,"X");
			}
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

double TTSESinglePhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *h, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE)
{	
	double val=0;
	h->resize(N);
	p->resize(N);
	EOS->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;
		
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

double TTSESinglePhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
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

double TTSESinglePhaseTableClass::evaluate(long iParam, double h, double p)
{
	int i = (int)round(((h-hmin)/(hmax-hmin)*(Nrow-1)));
	int j = (int)round(((p-pmin)/(pmax-pmin)*(Ncol-1)));
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

TTSETwoPhaseTableClass::TTSETwoPhaseTableClass(Fluid *pFluid, double Q)
{
	this->pFluid = pFluid;
	this->Q = Q;
}
void TTSETwoPhaseTableClass::set_size(unsigned int N)
{
	this->N = N;
	
	// Seed the generator
	srand((unsigned int)time(NULL));

	// Resize all the arrays
	h.resize(N);
	p.resize(N);
	logp.resize(N);
	T.resize(N);
	dTdp.resize(N);
	d2Tdp2.resize(N);
	rho.resize(N);
	logrho.resize(N);
	drhodp.resize(N);
	d2rhodp2.resize(N);
	s.resize(N);
	dsdp.resize(N);
	d2sdp2.resize(N);
	h.resize(N);
	dhdp.resize(N);
	d2hdp2.resize(N);
}

double TTSETwoPhaseTableClass::build(double pmin, double pmax, TTSETwoPhaseTableClass *other)
{
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);
	this->logpmax = log(pmax);

	double dlogp = (logpmax-logpmin)/(N-1);
	clock_t t1,t2;
	t1 = clock();
	// Logarithmic distribution of pressures
	for (unsigned int i = 0; i < N-1; i++)
	{
		// Calculate the pressure
		p[i] = exp(logpmin + i*dlogp);
		logp[i] = log(p[i]);
		// Update the class
		CPS.update(iP,p[i],iQ,Q);

		// rho for the other saturation curve
		double rho_other = (this->Q>0.5) ? CPS.rhoL() : CPS.rhoV();
		double Q_other = (this->Q>0.5) ? 0 : 1;

		// Set the variables
		T[i] = CPS.T();
		dTdp[i] = CPS.dTdp_along_sat();
		d2Tdp2[i] = CPS.d2Tdp2_along_sat();
		h[i] = CPS.h();
		dhdp[i] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
		d2hdp2[i] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
		s[i] = CPS.s();
		dsdp[i] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
		d2sdp2[i] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
		rho[i] = CPS.rho();
		logrho[i] = log(CPS.rho());
		drhodp[i] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
		d2rhodp2[i] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

		// If other is provided
		if (other != NULL)
		{
			other->pmin = pmin;
			other->pmax = pmax;
			other->logpmin = log(pmin);
			other->logpmax = log(pmax);

			// Set the variables
			other->T[i] = CPS.T();
			other->dTdp[i] = CPS.dTdp_along_sat();
			other->d2Tdp2[i] = CPS.d2Tdp2_along_sat();
			other->h[i] = (other->Q>0.5) ? CPS.hV() : CPS.hL();
			other->dhdp[i] = (other->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
			other->d2hdp2[i] = (other->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
			other->s[i] = (other->Q>0.5) ? CPS.sV() : CPS.sL();
			other->dsdp[i] = (other->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
			other->d2sdp2[i] = (other->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
			other->rho[i] = (other->Q>0.5) ? CPS.rhoV() : CPS.rhoL();
			other->logrho[i] = log(other->rho[i]);
			other->drhodp[i] = (other->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
			other->d2rhodp2[i] = (other->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();
		}
	}
	//// At the last point
	//CPS.flag_SinglePhase = true; // Don't have it check the state or do a saturation call
	//p[N-1] = CPS.pFluid->reduce.p;
	//logp[N-1] = log(p[N-1]);
	//CPS.update(iT,CPS.pFluid->reduce.T,iD,CPS.pFluid->reduce.rho);
	//T[N-1] = CPS.T();
	//dTdp[N-1] = CPS.dTdp_along_sat();
	//d2Tdp2[N-1] = CPS.d2Tdp2_along_sat();
	//h[N-1] = CPS.h();
	//dhdp[N-1] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
	//d2hdp2[N-1] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
	//s[N-1] = CPS.s();
	//dsdp[N-1] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
	//d2sdp2[N-1] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
	//rho[N-1] = CPS.rho();
	//logrho[N-1] = log(CPS.rho());
	//drhodp[N-1] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
	//d2rhodp2[N-1] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

	t2 = clock();
	std::cout << double(t2-t1)/CLOCKS_PER_SEC << " to build both two phase tables" << std::endl;
	return double(t2-t1)/CLOCKS_PER_SEC;
}

double TTSETwoPhaseTableClass::evaluate(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit below the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	// If the pressure is just barely above the critical pressure
	// or between the critical pressure and the next lowest point,
	// just just the next lowest point to avoid some of the 
	// derivatives that are infinite at the critical point
	if (i == N || i == N-1 ) i = N-2;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range, yielded index of %d",p,i));
	}
		
	double log_PI_PIi = logp-this->logp[i];
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iS:
		return s[i]+log_PI_PIi*pi*dsdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2sdp2[i]*pi*pi;
	case iT:
		return T[i]+log_PI_PIi*pi*dTdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2Tdp2[i]*pi*pi;
	case iH:
		//return h[i]+(pi-this->p[i])*dhdp[i]+0.5*(pi-this->p[i])*(pi-this->p[i])*d2hdp2[i];
		return h[i]+log_PI_PIi*pi*dhdp[i]*(1.0+0.5*log_PI_PIi)+0.5*log_PI_PIi*log_PI_PIi*d2hdp2[i]*pi*pi;
	case iD:
		{
		// log(p) v. log(rho) gives close to a line for most of the curve
		double dRdPI = pi/rho[i]*drhodp[i];
		return exp(logrho[i]+log_PI_PIi*(1.0+0.5*log_PI_PIi*(1-dRdPI))*dRdPI+0.5*log_PI_PIi*log_PI_PIi*d2rhodp2[i]*pi*pi/rho[i]);
		}
	default:
		throw ValueError();
	}
}
double TTSETwoPhaseTableClass::evaluate_sat_derivative(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit out of the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	if (i == N) i = N-1;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range",p));
	}
		
	double log_PI_PIi = logp-this->logp[i];
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iT:
		{
			/// First order expansion of dTdp around point of interest
			return dTdp[i]+(p-pi)*d2Tdp2[i];
		}
	case iH:
		{	
			/// First order expansion of dhdp around point of interest
			return dhdp[i]+(p-pi)*d2hdp2[i];
		}
	case iS:
		{	
			/// First order expansion of dsdp around point of interest
			return dsdp[i]+(p-pi)*d2sdp2[i];
		}
	case iD:
		{
			/// First order expansion of drhodp around point of interest
			return drhodp[i]+(p-pi)*d2rhodp2[i];
		}
	default:
		throw ValueError();
	}
}

double TTSETwoPhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
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


double TTSETwoPhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *p, std::vector<double> *EOS, std::vector<double> *TTSE)
{	
	double val=0;
	p->resize(N);
	EOS->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClass CPS = CoolPropStateClass(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		
		CPS.update(iP,p1,iQ,this->Q);
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
		
		std::cout << format("%g %g %g %g TTSE (p,EOS,TTSE, diff)\n",p1,(*EOS)[i],(*TTSE)[i],((*EOS)[i]-(*TTSE)[i]));
	}
	return val;
}

