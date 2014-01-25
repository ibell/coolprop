# -*- coding: utf-8 -*-
"""
Created on Sun May 19 22:23:19 2013

@author: Belli
"""
from __future__ import division
from math import exp

n=[0,    
		 0.125335479355233e-1, #[1]
		 0.78957634722828e1, #[2]
		-0.87803203303561e1, #[3]
		 0.31802509345418, #[4]
		-0.26145533859358, #[5]
		-0.78199751687981e-2, #[6]
		 0.88089493102134e-2, #[7]
		-0.66856572307965, #[8]
		 0.20433810950965, #[9]
		-0.66212605039687e-4, #[10]
		-0.19232721156002, #[11]
		-0.25709043003438, #[12]
		 0.16074868486251, #[13]
		-0.40092828925807e-1, #[14]
		 0.39343422603254e-6, #[15]
		-0.75941377088144e-5, #[16]
		 0.56250979351888e-3, #[17]
		-0.15608652257135e-4, #[18]
		 0.11537996422951e-8, #[19]
		 0.36582165144204e-6, #[20]
		-0.13251180074668e-11, #[21]
		-0.62639586912454e-9, #[22]
		-0.10793600908932, #[23]
		 0.17611491008752e-1, #[24]
		 0.22132295167546, #[25]
		-0.40247669763528, #[26]
		 0.58083399985759, #[27]
		 0.49969146990806e-2, #[28]
		-0.31358700712549e-1, #[29]
		-0.74315929710341, #[30]
		 0.47807329915480, #[31]
		 0.20527940895948e-1, #[32]
		-0.13636435110343, #[33]
		 0.14180634400617e-1, #[34]
		 0.83326504880713e-2, #[35]
		-0.29052336009585e-1, #[36]
		 0.38615085574206e-1, #[37]
		-0.20393486513704e-1, #[38]
		-0.16554050063734e-2, #[39]
		 0.19955571979541e-2, #[40]
		 0.15870308324157e-3, #[41]
		-0.16388568342530e-4, #[42]
		 0.43613615723811e-1, #[43]
		 0.34994005463765e-1, #[44]
		-0.76788197844621e-1, #[45]
		 0.22446277332006e-1, #[46]
		-0.62689710414685e-4, #[47]
		-0.55711118565645e-9, #[48]
		-0.19905718354408, #[49]
		 0.31777497330738, #[50]
		-0.11841182425981, #[51]
		-0.31306260323435e2, #[52]
		 0.31546140237781e2, #[53]
		-0.25213154341695e4, #[54]
		-0.14874640856724, #[55]
		 0.31806110878444, #[56]
	]

d=[0,
		1, #[1]
		1, #[2]
		1, #[3]
		2, #[4]
		2, #[5]
		3, #[6]
		4, #[7]
		1, #[8]
		1, #[9]
		1, #[10]
		2, #[11]
		2, #[12]
		3, #[13]
		4, #[14]
		4, #[15]
		5, #[16]
		7, #[17]
		9, #[18]
		10, #[19]
		11, #[20]
		13, #[21]
		15, #[22]
		1, #[23]
		2, #[24]
		2, #[25]
		2, #[26]
		3, #[27]
		4, #[28]
		4, #[29]
		4, #[30]
		5, #[31]
		6, #[32]
		6, #[33]
		7, #[34]
		9, #[35]
		9, #[36]
		9, #[37]
		9, #[38]
		9, #[39]
		10, #[40]
		10, #[41]
		12, #[42]
		3, #[43]
		4, #[44]
		4, #[45]
		5, #[46]
		14, #[47]
		3, #[48]
		6, #[49]
		6, #[50]
		6, #[51]
		3, #[52]
		3, #[53]
		3, #[54]
	]

t= [0.00,
		-0.5, #[1]
		0.875, #[2]
		1, #[3]
		0.5, #[4]
		0.75, #[5]
		0.375, #[6]
		1, #[7]
		4, #[8]
		6, #[9]
		12, #[10]
		1, #[11]
		5, #[12]
		4, #[13]
		2, #[14]
		13, #[15]
		9, #[16]
		3, #[17]
		4, #[18]
		11, #[19]
		4, #[20]
		13, #[21]
		1, #[22]
		7, #[23]
		1, #[24]
		9, #[25]
		10, #[26]
		10, #[27]
		3, #[28]
		7, #[29]
		10, #[30]
		10, #[31]
		6, #[32]
		10, #[33]
		10, #[34]
		1, #[35]
		2, #[36]
		3, #[37]
		4, #[38]
		8, #[39]
		6, #[40]
		9, #[41]
		8, #[42]
		16, #[43]
		22, #[44]
		23, #[45]
		23, #[46]
		10, #[47]
		50, #[48]
		44, #[49]
		46, #[50]
		50, #[51]
		0, #[52]
		1, #[53]
		4, #[54]
	]

c=[0,
		0, #[1]
		0, #[2]
		0, #[3]
		0, #[4]
		0, #[5]
		0, #[6]
		0, #[7]
		1, #[8]
		1, #[9]
		1, #[10]
		1, #[11]
		1, #[12]
		1, #[13]
		1, #[14]
		1, #[15]
		1, #[16]
		1, #[17]
		1, #[18]
		1, #[19]
		1, #[20]
		1, #[21]
		1, #[22]
		2, #[23]
		2, #[24]
		2, #[25]
		2, #[26]
		2, #[27]
		2, #[28]
		2, #[29]
		2, #[30]
		2, #[31]
		2, #[32]
		2, #[33]
		2, #[34]
		2, #[35]
		2, #[36]
		2, #[37]
		2, #[38]
		2, #[39]
		2, #[40]
		2, #[41]
		2, #[42]
		3, #[43]
		3, #[44]
		3, #[45]
		3, #[46]
		4, #[47]
		6, #[48]
		6, #[49]
		6, #[50]
		6, #[51]
		0, #[52]
		0, #[53]
		0, #[54]
	]

alpha=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 20,		20,		20]

beta=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,150,		150,
		250,
		0.3,
		0.3,
	]

GAMMA=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		1.21, #[52]
		1.21, #[53]
		1.25, #[54]
	]

epsilon=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		1,
		1,
		1,
	]

a=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		3.5,
		3.5,]

b=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.85,
		0.95,
	]

A=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.32, #[55]
		0.32, #[56]
	]

B=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		0.2, #[55]
		0.2, #[56]
	]

C=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
		28., #[55]
		32., #[56]
	]

D=[
		0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,700.0,800.0]
  
  
  
tau = 1.01
delta = 0.95
dtau = 1e-10
ddelta = 1e-10
for i in [55]:
    
    def f(tau, delta):
        
        theta = (1.0-tau)+A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]))
        DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
        PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
        dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
        dDELTA_dDelta=(delta-1.0)*(A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
        dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
        dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
        dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
        dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
        dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
        dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
        dDELTAbi2_dDelta_dTau=-A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
        dPSI3_dDelta2_dTau = (2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*dPSI_dTau;
        dDELTA_dTau = -2*theta;
        dDELTA2_dDelta_dTau = 2.0*A[i]/(beta[i])*pow(pow(delta-1,2),1.0/(2.0*beta[i])-0.5);
        dDELTA3_dDelta2_dTau = 2.0*A[i]*(beta[i]-1)/(beta[i]*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1.0);
        
        dDELTAbim1_dTau = (b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
        dDELTAbim2_dTau = (b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
        Line1 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
        Line2 = (b[i]-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
        dDELTAbi3_dDelta2_dTau = b[i]*(Line1+Line2);
        
        ddelta2 = n[i]*(pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);        
        
        dPSI3_dDelta3=2.0*C[i]*PSI*(-4*C[i]*C[i]*pow(delta-1.0,3)+6*C[i]*(delta-1));
        dtheta_dDelta = A[i]/(2*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1)*2*(delta-1);
        PI = 4*B[i]*a[i]*(a[i]-1)*pow(pow(delta-1,2),a[i]-2)+2*pow(A[i]/beta[i],2)*pow(pow(delta-1,2),1/beta[i]-2)+4*A[i]*theta/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2);
        dPI_dDelta = -8*B[i]*a[i]*(a[i]-1)*(a[i]-2)*pow(pow(delta-1,2),a[i]-5.0/2.0)-8*pow(A[i]/beta[i],2)*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/beta[i]-5.0/2.0)-(8*A[i]*theta)/beta[i]*(1/(2*beta[i])-1)*(1/(2*beta[i])-2)*pow(pow(delta-1,2),1/(2*beta[i])-5.0/2.0)+4*A[i]/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2)*dtheta_dDelta;
        
        dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
        dDELTAbi3_dDelta3 = b[i]*(pow(DELTA,b[i]-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+(b[i]-1)*(pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta));
  
        Line1 = pow(DELTA,b[i])*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
        Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
        Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        ddelta2_dtau = n[i]*(Line1+Line2+Line3);
        
        ddelta3 = n[i]*(pow(DELTA,b[i])*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
        
        return locals()
    
    
    base = f(tau, delta)
    plus_tau = f(tau+dtau, delta)
    minus_tau = f(tau-dtau, delta)
    plus_delta = f(tau, delta+ddelta)
    minus_delta = f(tau, delta-ddelta)
    
#    print base['dDELTA_dDelta'],(plus_delta['DELTA']-minus_delta['DELTA'])/(2*ddelta)    
#    print base['dDELTA_dTau'],(plus_tau['DELTA']-minus_tau['DELTA'])/(2*dtau)   
#    print base['dDELTA2_dDelta_dTau'],(plus_tau['dDELTA_dDelta']-minus_tau['dDELTA_dDelta'])/(2*dtau)
#    print base['dDELTA3_dDelta2_dTau'],(plus_delta['dDELTA2_dDelta_dTau']-minus_delta['dDELTA2_dDelta_dTau'])/(2*ddelta)
#    print base['dDELTAbi3_dDelta2_dTau'],(plus_delta['dDELTAbi2_dDelta_dTau']-minus_delta['dDELTAbi2_dDelta_dTau'])/(2*ddelta)
#    print base['ddelta2_dtau'],(plus_tau['ddelta2']-minus_tau['ddelta2'])/(2*dtau)
#    print base['dPSI3_dDelta3'],(plus_delta['dPSI2_dDelta2']-minus_delta['dPSI2_dDelta2'])/(2*ddelta)
#    print base['dtheta_dDelta'],(plus_delta['theta']-minus_delta['theta'])/(2*ddelta)
    print base['dPI_dDelta'],(plus_delta['PI']-minus_delta['PI'])/(2*ddelta)
    print base['dDELTA3_dDelta3'],(plus_delta['dDELTA2_dDelta2']-minus_delta['dDELTA2_dDelta2'])/(2*ddelta)    
    print base['dDELTAbi3_dDelta3'],(plus_delta['dDELTAbi2_dDelta2']-minus_delta['dDELTAbi2_dDelta2'])/(2*ddelta)
    print base['ddelta3'],(plus_delta['ddelta2']-minus_delta['ddelta2'])/(2*ddelta)