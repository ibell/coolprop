
#ifndef REFPROP_H
#define REFPROP_H

	// Only add REFPROP if build on Windows or Linux platform
	#if defined(__ISWINDOWS__)||defined(__ISLINUX__)
	double REFPROP(char Output,        char Name1,        double Prop1, char Name2,        double Prop2, char * Ref);
	double REFPROP(std::string Output, std::string Name1, double Prop1, std::string Name2, double Prop2, std::string Ref);
	#endif

#endif
