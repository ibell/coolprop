#ifndef _STRUCTS_MACROS_H
#define _STRUCTS_MACROS_H

	#define Ntheta_MAX 50000
	#define T_p_xL	0
	#define T_m_xL	1
	#define true 1
	#define false 0
	#define NTHETA_ADISC 100



	// ***************************************************************
	//                      Structure Prototypes 
	// ***************************************************************

	struct phiVals{
		double phi_fi0;
		double phi_fis;
		double phi_fie;
		double phi_fo0;
		double phi_fos;
		double phi_foe;
		double phi_oi0;
		double phi_ois;
		double phi_oie;
		double phi_oo0;
		double phi_oos;
		double phi_ooe;
	};
	
	struct discVals{
		double x0;
		double y0;
		double R;

		double xa_arc1;
		double ya_arc1;
		double ra_arc1;
		double t1_arc1;
		double t2_arc1;

		double m_line;
		double b_line;
		double t1_line;
		double t2_line;

		double xa_arc2;
		double ya_arc2;
		double ra_arc2;
		double t1_arc2;
		double t2_arc2;

		char Type[100];

		double thetaAdisc[NTHETA_ADISC];
		double Adisc[NTHETA_ADISC];
	};

	struct suctVals{
		double A_sa_suction;
	};

	struct wallVals{
		double x0;
		double y0;
		double r;
	};
	
	struct flowModelVals{
		int flank;
		int radial;
		int s_sa;
		int suction;
		int d_dd;
		int discharge;
	};
	
	struct geoVals{
		struct phiVals phi;
		struct discVals disc;
		struct suctVals suct;
		struct wallVals wall;
		struct flowModelVals flowModels;
		double rb;
		double ro;
		double t;
		double hs;
		double delta_flank;
		double delta_radial;
	};

	struct flowVecVals{
		int *CV1;
		int *CV2;
		int *CVup;
		double *A;
		int *flowModel;
		double *mdot;
		double *mdot_L;
		double *mdot_g;
		double *mdot_c;
		double *p_up;
		double *p_down;
		double *h_up;
		double *h_down;
		double *T_up;
		double *T_down;
		double *xL;
		double *Ed;
		double *Re;
		double *Ma;
		int N;
	};

	struct scrollInputVals 
	{
		double T_in;
		double T_out;
		double T_amb;
		double p_in;
		double p_out;
		double xL_in;
		char Ref[200];	
		char Liq[200];
		double omega;
	};
	struct ExperVals
	{
		double mdot;
		double T_d;
		double P_shaft;
		double eta_c;
		double eta_v;
	};
	struct FlowVals
	{
		double w_ent;
		double sigma; //Area ratio used for nozzle model
		double Z_D_bends;
		double L_inlet;
		double D_inlet; // Inlet and outlet ports are assumed to be same diameter and length
		double L_flank;
		double delta_flank;
		double delta_radial;
		double phi_flank;
		int flowModel_flank;
		int flowModel_radial;
		int flowModel_s_sa;
		int flowModel_d_dd;
		int flowModel_suction;
		int flowModel_discharge;
	};
	struct MLVals
	{
		double m;
		double b;
		double c;
		double eta_m;
		char Type[200];
		double UA_amb;
		double etac_guess;
	};
	struct flagVals
	{
		int useDDD;
		int LeftDischarge; // At theta just to the left of the discharge angle
		int lastRotation;
	};
	struct PowerEffVals
	{
		double P_gas;
		double P_shaft;
		double P_ML;
		double eta_m;
		double eta_c;
		double eta_v;
	};
	struct HTVals
	{
		double T_scroll;
		double T_amb;
		double *hc;
		double *A_wall_i;
		double *A_wall_o;
		double *Tm_plate;
		double *Tm_wall_i;
		double *Tm_wall_o;
		double Q_scroll_gas;
		double Q_scroll_inlet;
		double Q_scroll_outlet;
		double Q_scroll_plenum;
		double Q_scroll_amb;
	};

	struct massFlowVals
	{
		struct FlowVals Inputs;
		double mdot_tot;
	};

	struct DebugVals
	{
		double wrap_error_rel;
		double mdot_error_abs;
		int Ntheta;
		double ElapsedTime;
	};

	struct LossesVals
	{
		double suction;
		double discharge;
		double leakage_flank;
		double leakage_radial;
		double mechanical;
		double Wdot_adiabatic;
	};

	struct scrollVals
	{
		struct geoVals geo;
		struct flagVals flags;
		struct HTVals HT;
		struct MLVals ML;
		struct massFlowVals massFlow;
		struct flowVecVals *flowVec;
		struct scrollInputVals Inputs;
		struct ExperVals Exper;
		struct PowerEffVals PowerEff;
		struct DebugVals Debug;
		struct LossesVals Losses;
		double *V;
		double *dV;
		double *T;
		double *p;
		double *xL;
		double *m;
		double *Q;
		double *rho;
		double *error;
		double *theta;
		char Ref[200];						/*	The working fluid	  */
		char Liq[200];						/*	The flooding liquid	  */
		double omega;
		int N;
		int Ntheta;
	};


#endif



#ifdef HUGE_VAL
#define _HUGE HUGE_VAL
#else
	// GCC Version of huge value macro
	#ifdef HUGE 
	#define _HUGE HUGE
	#endif
#endif


