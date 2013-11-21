

#ifndef CPEXCEPTIONS_H
#define CPEXCEPTIONS_H

#include <exception>
#include <iostream>

class CoolPropBaseError: public std::exception
{
protected:
	std::string err; // Can be accessed by subclasses since it is protected
public:
	CoolPropBaseError() throw() {};
	~CoolPropBaseError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class NotImplementedError: public CoolPropBaseError
{
public:
	NotImplementedError() throw() {};
	NotImplementedError(std::string errstring) throw(){err=errstring;};
    ~NotImplementedError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class SolutionError: public CoolPropBaseError
{
public:
	SolutionError()throw() {};
	SolutionError(std::string errstring) throw(){err=errstring;};
    ~SolutionError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class ValueError: public CoolPropBaseError
{
public:
	ValueError() throw() {};
	ValueError(std::string errstring) throw(){err=errstring;};
    ~ValueError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

class AttributeError: public CoolPropBaseError
{
public:
	AttributeError() throw() {};
	AttributeError(std::string errstring) throw() {err=errstring;};
    ~AttributeError() throw() {};
	virtual const char* what() const throw(){ return err.c_str(); }
};

#endif
