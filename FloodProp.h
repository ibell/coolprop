/* File FloodProp.h */
 #ifndef FLOODPROP_H
 #define FLOODPROP_H
 
	//Define macros for different types of input parameters for EOS
	#define TP 1
	#define Trho 2

	//  ****************************
	//  ***** Gas Variables ********
	//  ****************************

	double cp_A[6], cp_B[6], cp_C[6], cp_D[6], cp_E[6];
	double kg_A[7], kg_B[7], kg_C[7];
	double w[6], Pc[6], Tc[6], MM_g[6];
	double mug_A[7], mug_B[7], mug_C[7];

	//  ****************************
	//  ***** Liq Variables ********
	//  ****************************
	#define NL 8

	double rhol_A[NL], rhol_B[NL], rhol_n[NL], rhol_Tc[NL];
	double kl_A[NL], kl_B[NL], kl_C[NL], kl_D[NL];
	double cl_A[NL], cl_B[NL], cl_C[NL], cl_D[NL];
	double mul_A[NL], mul_B[NL], mul_C[NL], mul_D[NL];
	double MM_l[NL];
	
	static int I_N2=0, I_He=1, I_Ne=2, I_Ar=3, I_Kr=4, I_Xe=5, I_CO2=6;
	static int I_Methanol=0, I_Ethanol=1, I_Propanol=2, I_Butanol=3, I_Water=4, I_NH3=5, I_Zerol=6,I_POE=7;

	double hm2(double T, double P, double xL);
	double cK_e(double v_l, double v_g, double x, double w, double flag);
	double cv_e(double v_l, double v_g, double K_e, double x, double w, double flag);
	double R(char *Gas);
	double rho_l(char *Liq, double T);
	double rho_m(char *Gas, char *Liq, double T, double P, double xL);
	double u_l(char *Liq, double T);
	double h_g(char *Gas,double T, double P);
	double u_m(char *Gas, char *Liq, double T, double P, double xL);
	double h_m(char *Gas, char *Liq, double T, double P, double xL);
	double c_v(char *Gas, double T, double P);
	double c_l(char *Liq, double T);
	double dudT_m(char *Gas, char *Liq, double T, double P, double xL);
	double dudP_m(char *Gas, char *Liq, double T, double P, double xL);
	double dudxL_m(char *Gas, char *Liq, double T, double P, double xL);
	double drhodP_m(char *Gas, char *Liq, double T, double P, double xL);
	double drhodT_m(char *Gas, char *Liq, double T, double P, double xL);
	double drhodxL_m(char *Gas, char *Liq, double T, double P, double xL);
	double cp_mix(char *Gas, char *Liq, double T, double P, double xL);
	double mu_mix(char *Gas, char *Liq, double T, double p, double xL);
	double k_mix(char *Gas, char *Liq, double T, double P, double xL);
	double Pr_mix(char *Gas, char *Liq, double T,double P,double xL);
	double s_m(char *Gas, char *Liq, double T, double P, double xL);
	double e_m(char *Gas, char *Liq, double T, double P, double xL);
	double dvdT_m(char *Gas, char *Liq, double T, double P, double xL);
	double dvdP_m(char *Gas, char *Liq, double T, double P, double xL);
	int getIndex(char *Fluid);
	double mu_l(char *Liq, double T);
	double mu_g(char *Gas, double T, double p);
	double k_l(char *Liq, double T);
	double k_g(char *Gas, double T, double p);
	double s_l(char *Liq, double T);
	double s_g(char *Gas, double T, double P);
	double VoidFrac(char *Gas, char *Liq, double T, double P, double xL);
	double c_p(char *Gas, double T, double P);
	double kstar_m(char *Gas, char *Liq, double T, double P, double xL);
	double rho_g(char *Gas, double T, double P);
	double u_g(char *Gas, double T, double P);
	double T_hp(char *Gas, char *Liq, double h, double p, double xL, double T_guess);
	double T_Up(char *Gas, char *Liq, double U, double p, double xL, double V, double T_guess);
	double h_sp(char *Gas, char *Liq, double s, double p, double xL, double T_guess);
	double T_sp(char *Gas, char *Liq, double s, double p, double xL, double T_guess);
	double p_Trho(char *Gas, char *Liq, double rho, double T, double xL, double p_guess);
	/*double c_g(char *Gas, double T, double p)*/;
	int isNAN_FP(double x);
	int isINFINITY_FP(double x);

	void setGas();
	void setLiq();

	double dpdT_const_v(char *Gas, char *Liq, double T, double p1, double xL);

 #endif 

