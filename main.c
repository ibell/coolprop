//#include "R134a.h"

#include "string.h"
#include <stdio.h>
#include "PropErrorCodes.h"
#include "PropMacros.h"
#include "R717.h"
#include "R410A.h"

int main()
{
	int c;
	double rho,T,v,k;
	
	/*rho=rho_R134a(300,400,TYPE_TP);
	printf("%d\n",errCode_R134a());
	rho=rho_R134a(10000,400,TYPE_TP);
	printf("%d\n",errCode_R134a());
	printf("rho: %g\n",rho);*/

	T=240;
	rho=rhosatL_R717(T,TYPE_TPNoLookup);
	printf("rho: %g\n",rho);

	v=visc_R717(T,rho,TYPE_Trho);
	printf("visc: %g\n",v);

	k=k_R717(T,rho,TYPE_Trho);
	printf("k: %g\n",k);

	printf("%d\n",errCode_R717());



	T=240;
	rho=rhosatV_R410A(T,TYPE_TPNoLookup);
	printf("rho: %g\n",rho);

	v=visc_R410A(T,rho,TYPE_Trho);
	printf("visc: %g\n",v);

	k=k_R410A(T,rho,TYPE_Trho);
	printf("k: %g\n",k);
	
	printf("Any key to end\n");
	c=getchar();

}