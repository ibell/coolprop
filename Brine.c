
#if defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <crtdbg.h>
#else
#include <stdlib.h>
#endif

#include "math.h"
#include "stdio.h"
#include <string.h>

static double powInt(double x, int y);

// inputs in C, %,
// outputs in kg/m^3, J/kg-K, mW/m-K, Pa-s
int Brine(char * Mix, double T, double C, /*in --- out */double *Tfreeze, double *Tmax, double *rho, double *cp, double *k, double *visc)
{

	static double PG[18][5]={
		{-23.6,1042,3679,0.3806,2.274},
		{-0.000297,-0.4907,1.571,0.0005765,-0.05342},
		{0.0000001871,-0.002819,0.01331,-0.0000003477,0.0005372},
		{0.000000159,-0.0000005895,0.0000001975,-0.000000006041,-0.000004955},
		{-1.088,0.8081,-19.33,-0.003815,0.045},
		{0.00004693,-0.009652,0.1118,-0.00001423,-0.0005488},
		{-0.0000009089,0.00007168,-0.001108,-0.00000001203,0.000001845},
		{-0.0000000009845,0.0000002404,0.000004924,-0.0000000005854,0.0000001192},
		{-0.02815,-0.007156,-0.04879,0.00000842,-0.00007808},
		{0.000001613,0.0001088,-0.0002338,0.0000001081,0.000001453},
		{0.00000005066,-0.000003328,0.00002753,0.000000001959,-0.0000002816},
		{-0.000000001251,0.0000001153,-0.0000003148,0.0000000001271,0.000000008562},
		{-0.0005285,0.000119,0.004749,-0.00000111,0.000006565},
		{-0.0000001786,-0.000006226,-0.00002621,-0.000000001612,-0.0000004032},
		{0.000000002051,-0.00000003026,0.000001286,0.0000000003005,-0.000000001212},
		{0.000009794,-0.0000117,-0.0002871,0.000000005503,0.0000006441},
		{-0.000000005668,-0.0000002915,-0.0000000905,0.0000000001437,-0.0000000143},
		{0.000000427,-0.0000006033,-0.00001068,0.00000000129,0.00000001092}
	};

	static double EG[18][5]={
		{-21.96,1056,3501,0.4211,1.453},
		{0.000009186,-0.3987,3.954,0.0007995,-0.03747},
		{-0.000001136,-0.003068,0.00006065,-0.00000005509,0.0002842},
		{0.00000001798,0.00001233,-0.000005979,-0.0000000146,-0.0000008025},
		{-1.015,1.505,-24.19,-0.003694,0.0292},
		{-0.000003391,-0.008953,0.1031,-0.00001751,-0.0001131},
		{0.00000007548,0.00006378,0.00004312,0.00000006656,0.000001729},
		{0.0000000008493,-0.0000001152,0.000005168,0.000000002017,-0.00000005073},
		{-0.01414,-0.001634,0.004613,0.00002095,0.0001264},
		{0.0000002355,0.0001541,-0.00006595,0.0000002078,0.000000006785},
		{0.000000002159,-0.000001874,0.0000162,-0.000000002394,-0.00000001685},
		{-0.0000000001311,-0.000000009809,-0.000000325,-0.00000000006772,-0.000000001082},
		{-0.00003727,-0.0002317,0.006028,0.0000003663,0.000004386},
		{0.0000000006091,0.000002549,0.00005642,-0.000000005272,-0.0000002191},
		{-0.0000000001465,-0.00000005523,-0.0000007777,-0.0000000001126,-0.00000000009117},
		{-0.0000001398,-0.00000851,-0.00007977,-0.000000006389,-0.00000009223},
		{-0.0000000002269,-0.00000003848,0.000000519,-0.0000000001112,-0.000000004294},
		{-0.00000002512,-0.0000001128,-0.00000338,-0.000000000182,-0.000000003655},
	};

	static double EA[18][5]={
		{-28.42,954.4,3925,0.3545,2.214},
		{0.000009753,-0.6416,3.876,0.0004421,-0.0571},
		{-0.00001236,-0.002495,0.00023,-0.0000002942,0.0004679},
		{0.0000006378,0.00001729,0.00001322,-0.00000001115,-0.000001374},
		{-0.8563,-1.729,-27.95,-0.004334,0.0008025},
		{0.00005274,-0.01824,0.1773,-0.00002021,0.0002618},
		{0.000001843,0.0003116,0.00004769,-0.000000004865,-0.000008472},
		{-0.0000001428,-0.0000006425,0.000003008,0.0000000002972,0.0000001478},
		{0.00405,-0.02193,-0.0962,0.00003021,-0.000733},
		{-0.000003058,0.0005847,-0.003908,0.0000004239,0.000007056},
		{-0.0000001531,-0.000002517,0.00001951,0.000000001007,0.0000002473},
		{0.000000005543,-0.00000002875,0.00000003366,-0.000000000007325,-0.00000001329},
		{-0.0001179,0.0006217,0.00758,0.0000006904,0.0000004285},
		{-0.00000009416,0.000004208,0.00002283,-0.000000003203,0.0000003239},
		{0.000000004676,-0.000000346,-0.0000009149,-0.00000000001439,-0.00000001234},
		{-0.000001992,0.000002288,-0.0001213,-0.00000001512,0.00000004313},
		{0.000000005409,-0.0000004141,0.000002545,-0.0000000003486,0.000000008582},
		{0.0000002951,-0.0000006412,0.0000002235,-0.000000001012,0.000000007654}
	};
	static double Glycerol[18][5]={
		{-21.38,1127,3211,0.4145,2.191},
		{-0.0002444,-0.4447,6.014,0.0007262,-0.04289},
		{0.000003632,-0.0014,-0.001326,0.0000005127,0.0004082},
		{0.00000005741,0.0000193,0.0001286,-0.000000008461,-0.000004723},
		{-0.8642,2.759,-21.48,-0.00303,0.05512},
		{0.00002169,-0.008192,0.1291,-0.00001776,-0.0004473},
		{-0.0000003779,0.000149,0.00002542,0.00000006843,0.000004772},
		{0.000000001157,-0.00000005755,-0.000002559,0.000000001317,0.00000004525},
		{-0.0154,-0.01341,0.09064,0.00001202,0.0001183},
		{0.000001725,0.00001003,-0.004195,0.0000002484,-0.000002661},
		{-0.0000000007435,-0.000004063,0.000004541,0.0000000009114,0.00000009049},
		{-0.0000000004038,0.00000004613,-0.0000005808,-0.00000000001729,0.000000007788},
		{-0.0004301,0.0001226,-0.003417,-0.0000005136,0.000004659},
		{-0.00000008476,0.000003235,0.0001241,-0.000000002727,-0.0000006555},
		{0.0000000003455,-0.0000002545,-0.0000004252,-0.00000000004282,0.000000004048},
		{-0.00001378,0.00005531,-0.0001272,-0.00000002006,0.000001682},
		{-0.000000003722,-0.00000001115,0.000007831,-0.0000000003081,-0.00000002797},
		{-0.0000001875,0.000001224,0.000002625,0.00000000006019,0.00000004123}
	};

	static double K2CO3[18][5]={
		{-15.68,1280,3008,0.5455,1.222},
		{0.000206,-0.454,1.565,0.001629,-0.02823},
		{-0.00000564,-0.002752,-0.002776,0.0000001773,0.0002125},
		{-0.00000006561,0.00005248,0.00004242,0.000000007906,-0.0000005715},
		{-1.157,11.18,-32.91,-0.0011,0.04577},
		{-0.00001398,-0.0004651,0.007871,-0.000002182,-0.0002322},
		{-0.0000005322,0.00007013,0.0003626,-0.0000002358,0.000002438},
		{0.00000002323,-0.000001974,0.00001842,-0.0000000009042,-0.0000001208},
		{-0.04603,0.009163,0.6281,-0.00002304,0.0007762},
		{-0.000001274,0.0004144,0.001919,-0.0000002244,-0.00001477},
		{0.0000001011,0.000001759,0.00001417,0.00000001281,0.0000002174},
		{-0.000000001636,0.0000000737,-0.000002517,-0.0000000001298,-0.000000005302},
		{-0.0009071,-0.0005183,0.01198,-0.0000006505,0.000006944},
		{0.00000003652,-0.00005161,0.000008158,0.00000001117,0.0000004662},
		{0.000000001217,0.0000003456,-0.000006241,0.0000000008923,0.000000002624},
		{0.00004312,0.0001101,-0.001027,0.00000005789,0.0000003892},
		{0.000000001571,-0.000001652,-0.000009209,0.0000000009352,0.00000001933},
		{0.000001577,0.000003241,-0.00004063,0.000000002245,0.00000001149}
	};

	static double Methanol[18][5]={
		{-29.03,959.2,3787,0.4017,1.548},
		{-0.0006352,-0.4154,2.599,0.0006268,-0.04187},
		{0.00000698,-0.002414,-0.001128,0.0000002004,0.0002657},
		{0.0000006072,0.000008179,0.00007089,0.00000004663,-0.000001831},
		{-1.237,-1.307,-23.25,-0.003897,0.000751},
		{0.00001642,-0.01618,0.1823,-0.00002469,0.0001718},
		{0.0000005963,0.0001165,-0.0001854,0.00000004831,-0.0000005091},
		{-0.00000003149,-0.000001234,0.00001831,-0.000000000669,-0.00000009014},
		{-0.01318,-0.0171,0.07125,0.00001292,-0.0004126},
		{0.000004175,0.0001804,0.001137,0.0000008417,-0.000003052},
		{-0.00000008982,-0.0000001331,0.00008856,-0.000000008036,-0.00000004946},
		{-0.000000000918,0.00000004142,-0.000001894,-0.0000000002734,0.000000004482},
		{-0.00008161,-0.0001412,0.0001995,-0.0000003086,0.000004679},
		{-0.00000007024,0.000004157,0.0001754,-0.00000002512,0.0000001484},
		{0.0000000000975,0.00000007795,-0.000005306,0.0000000002278,-0.000000000953},
		{-0.000001414,0.00001192,-0.0003196,0.000000003034,-0.0000002455},
		{-0.000000004866,0.00000001199,-0.000009812,-0.0000000005562,0.00000001915},
		{-0.00000005328,0.00000009182,-0.000004758,0.000000000446,-0.000000002372}
	};

	static double MgCl2[18][5]={
		{-14.3,1120,3381,0.5412,1.04},
		{0.0003996,-0.2119,2.296,0.001749,-0.03494},
		{-0.00001789,-0.003251,-0.002678,0.0000009893,0.0002209},
		{0,0,0,0,0},
		{-1.88,9.281,-51.22,-0.002517,0.06105},
		{0.000009842,-0.01125,0.07385,-0.00000001043,0.00008266},
		{-0.000001506,0.0001878,-0.00002043,0.00000009471,-0.000004232},
		{0,0,0,0,0},
		{-0.07412,0.007164,0.371,-0.00004859,0.000313},
		{-0.00000783,-0.0006966,-0.01552,-0.0000004932,0.00001185},
		{0.0000003937,0.00001006,0.0001504,-0.00000003196,-0.00000005275},
		{0,0,0,0,0},
		{0.00001866,-0.005212,-0.07234,-0.000004771,-0.00009447},
		{-0.0000004194,0.00009488,0.001424,-0.00000006611,-0.000003069},
		{0.00000003098,0.0000006481,-0.000009528,-0.0000000005892,0.00000001665},
		{0.000179,0.0003832,0.004349,0.0000005374,0.00001709},
		{0.000000004594,0.00001119,0.00008803,0.00000000004494,-0.0000001421},
		{0.000007885,0.00004374,0.0006,0.00000005168,0.000001484}
	};

	static double NaCl[18][5]={
		{-8.45,1094,3619,0.5692,0.4951},
		{-0.0001238,-0.3624,1.893,0.001677,-0.02743},
		{0.000003543,-0.001457,-0.0002804,-0.000002661,0.0002397},
		{0,0,0,0,0},
		{-0.8559,6.958,-33.84,-0.0008528,0.02277},
		{0.000004486,-0.01599,0.06473,-0.00001519,-0.000009952},
		{-0.0000002517,0.00005737,-0.001467,0.0000003244,0.000004419},
		{0,0,0,0,0},
		{-0.0211,0.02889,0.7992,-0.000009082,0.0004907},
		{0.000000736,0.0003118,-0.01458,-0.00000004241,-0.000009974},
		{0,0,0,0,0},
		{0,0,0,0,0},
		{-0.0005851,0.006133,-0.01959,-0.0000003147,-0.000002524},
		{0,0,0,0,0},
		{0,0,0,0,0},
		{0,0,0,0,0},
		{0,0,0,0,0},
		{0,0,0,0,0}
	};

	static double NH3[18][5]={
		{-30.57,939.9,4240,0.4418,0.952},
		{-0.0001414,-0.3043,-1.679,0.001753,-0.03449},
		{0.000002449,-0.002693,0.01681,0.000000571,0.0002861},
		{0.0000000374,0.000008072,0.000002898,-0.000000001917,-0.000002255},
		{-2.715,-2.886,1.168,-0.004573,0.009597},
		{-0.00002229,-0.02814,0.03685,0.000001102,-0.000005456},
		{-0.0000004379,0.0002623,0.001072,-0.0000001749,0.000003612},
		{0.00000006155,-0.0000004941,-0.00004204,-0.000000002836,-0.00000007025},
		{-0.1079,0.005131,0.02129,0.0001131,0.0001751},
		{0.00001571,-0.0001273,0.003609,-0.0000007165,0.00001688},
		{0.0000001631,-0.00002052,-0.00007538,-0.000000003195,-0.0000001284},
		{-0.00000001285,0.0000003845,-0.00000173,0.0000000003352,0.000000004655},
		{-0.004386,0.001575,-0.0008558,0.000006136,-0.0000128},
		{-0.00000007756,0.00005527,0.000137,0.0000001029,-0.0000001217},
		{-0.00000003555,-0.0000002289,-0.00000965,0.000000001549,-0.00000005541},
		{0.0002768,0.0005629,0.0004583,-0.000001119,0.0000002954},
		{-0.0000001975,0.000006525,-0.000008621,0.00000001056,-0.000000217},
		{0.00003206,0.00002825,0.00004901,-0.0000001582,0.0000003166}
	};

	static double CaCl2[18][5]={
		{-23.88,1216,2940,0.5423,1.293},
		{0.0002256,-0.3703,2.991,0.001455,-0.02799},
		{-0.0000031,-0.003431,-0.01335,0.0000003436,0.0001994},
		{-0.00000008525,0.0000405,0.00001808,-0.00000002934,-0.000001548},
		{-2.158,10.25,-37.73,-0.001068,0.0517},
		{-0.00006304,-0.01716,-0.01591,-0.00001181,0.0001623},
		{0.00000117,0.0003797,0.001381,-0.0000001678,0.000005934},
		{-0.0000000009079,-0.000002133,0.000004197,0.000000004407,-0.0000001066},
		{-0.103,-0.01559,0.9467,-0.00003676,0.00097},
		{-0.000004898,-0.0003008,-0.00866,-0.0000006941,-0.00002344},
		{-0.00000005161,0.00001497,0.0002843,0.00000001032,0.0000003511},
		{0.000000002914,-0.0000001062,0.0000008145,0.0000000002433,0.000000009953},
		{-0.007262,0.001609,0.006369,-0.000004215,0.00001344},
		{0.000001012,0.00003185,-0.0001969,0.00000006028,-0.000003604},
		{-0.000000006249,-0.000003868,0.000008495,0.000000001223,0.00000003177},
		{-0.0004101,0.001165,-0.001997,0.0000004362,0.0000006513},
		{0.00000006819,-0.000002937,-0.000003355,0.000000007269,-0.0000001297},
		{-0.000008832,0.00006015,-0.0001423,0.00000004104,-0.0000000124}
	};
	
	static int a[18][2]={
		{0,0},
		{0,1},
		{0,2},
		{0,3},
		{1,0},
		{1,1},
		{1,2},
		{1,3},
		{2,0},
		{2,1},
		{2,2},
		{2,3},
		{3,0},
		{3,1},
		{3,2},
		{4,0},
		{4,1},
		{5,0}
	};

	double (*A)[18][5];
	double xm,ym,x,y;
	double f_Tfreeze=0.0, f_rho=0.0,f_cp=0.0,f_k=0.0,f_visc=0.0;
	/*double Tfreeze, rho,cp,k,visc;*/
	double Cmin,Cmax;
	int i;

	if (!strcmp(Mix,"EG"))
		{ A=&EG; xm=38.1615; ym=6.3333; Cmin=0.0; Cmax=56.0; *Tmax=40;}
	if (!strcmp(Mix,"PG"))
		{ A=&PG; xm=42.7686; ym=5.3571; Cmin=15.0; Cmax=57.0; *Tmax=40;}
	if (!strcmp(Mix,"EA"))
		{ A=&EA; xm=38.9250; ym=-4.9038; Cmin=11.0; Cmax=60.0; *Tmax=20;}
	if (!strcmp(Mix,"NH3"))
		{ A=&NH3; xm=17.9664; ym=-7.1429; Cmin=7.8; Cmax=23.6; *Tmax=20;}
	if (!strcmp(Mix,"K2CO3"))
		{ A=&K2CO3; xm=27.6708; ym=4.1667; Cmin=0.0; Cmax=39.0; *Tmax=30;}
	if (!strcmp(Mix,"Methanol"))
		{ A=&Methanol; xm=32.9283; ym=-6.25; Cmin=7.8; Cmax=47.4; *Tmax=20;}
	if (!strcmp(Mix,"Glycerol"))
		{ A=&Glycerol; xm=47.8467; ym=7.000; Cmin=19.5; Cmax=63; *Tmax=40;}
	if (!strcmp(Mix,"CaCl2"))
		{ A=&CaCl2; xm=22.9180; ym=0.2459; Cmin=9.0; Cmax=29.4; *Tmax=30;}
	if (!strcmp(Mix,"MgCl2"))
		{ A=&MgCl2; xm=13.7000; ym=5.8750; Cmin=0.0; Cmax=20.5; *Tmax=30;}
	if (!strcmp(Mix,"NaCl"))
		{ A=&NaCl; xm=12.3539; ym=9.2581; Cmin=0.0; Cmax=23.0; *Tmax=30;}
	
	if (C<Cmin || C>Cmax)
	{
		printf("Conc.(%g %%) out of range [%g,%g]",C,Cmin,Cmax);
		*Tfreeze=-1;
		*Tmax=-1;
		*rho=-1;
		*cp=-1;
		*k=-1;
		*visc=-1;
		return -1;
	}

	for(i=0;i<18;i++)
	{
		x=C; y=T;
		f_Tfreeze 	+= (*A)[i][0] * powInt((x-xm),a[i][0]) * powInt((y-ym),a[i][1]);
		f_rho		+= (*A)[i][1] * powInt((x-xm),a[i][0]) * powInt((y-ym),a[i][1]);
		f_cp		+= (*A)[i][2] * powInt((x-xm),a[i][0]) * powInt((y-ym),a[i][1]);
		f_k			+= (*A)[i][3] * powInt((x-xm),a[i][0]) * powInt((y-ym),a[i][1]);
		f_visc		+= (*A)[i][4] * powInt((x-xm),a[i][0]) * powInt((y-ym),a[i][1]);
	}

	if (T < f_Tfreeze || T > *Tmax)
	{
		printf("Brine Temp.(%g C) out of range [%g,%g] for fluid %s\n",T,f_Tfreeze,*Tmax,Mix);
		*Tfreeze=f_Tfreeze;
		*rho=-2;
		*cp=-2;
		*k=-2;
		*visc=-2;
		return -1;
	}
	else
	{
		*Tfreeze=f_Tfreeze;
		*rho=f_rho;
		*cp=f_cp;
		*k=f_k;
		*visc=exp(f_visc);
		return 0;
	}
}

static double powInt(double x, int y)
{
    int i;
    double product=1.0;
    double x_in;
    int y_in;
    
    if (y==0)
    {
        return 1.0;
    }
    
    if (y<0)
    {
        x_in=1/x;
        y_in=-y;
    }
	else
	{
		x_in=x;
		y_in=y;
	}

    if (y_in==1)
    {
        return x_in;
    }    
    
    product=x_in;
    for (i=1;i<y_in;i++)
    {
        product=product*x_in;
    }
    
    return product;
}
