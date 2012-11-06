#include <exception>
#include <iostream>

#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

class NotImplementedError: public std::exception
{
private:
	std::string err;
public:
	NotImplementedError() throw(){};
	NotImplementedError(std::string errstring) throw(){err=errstring;};
    ~NotImplementedError() throw(){};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class SolutionError: public std::exception
{
private:
	std::string err;
public:
	SolutionError() throw(){};
	SolutionError(std::string errstring) throw(){err=errstring;};
    ~SolutionError() throw(){};
	virtual const char* what() const throw(){ return err.c_str(); }
};


class ValueError: public std::exception
{
private:
	std::string err;
public:
	ValueError() throw(){};
	ValueError(std::string errstring) throw(){err=errstring;};
    ~ValueError() throw(){};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class AttributeError: public std::exception
{
private:
	std::string err;
public:
	AttributeError() throw(){};
	AttributeError(std::string errstring) throw(){err=errstring;};
    ~AttributeError() throw(){};
	virtual const char* what() const throw(){ return err.c_str(); }
};

#endif