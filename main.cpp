#include "CoolProp.h"
#include <iostream>
#include <stdlib.h>

int main()
{
	double T = Props("T","H",246.532409342343,"P",1896.576573868160,"R410A");
	std::cout << T << std::endl;
}