#ifndef MPLPLOT_H
#define MPLPLOT_H

#include <Python.h>
#include <string>
#include <map>
#include <vector>
#include "MatrixMath.h"


/*! A C++ equivalent of a dictionary to hold values for the plotting functions
*/
class Dictionary
{
public:
	Dictionary(void){};
	std::map<std::string, std::string> string_values;
	std::map<std::string, double> double_values;

	/// Add a floating point value to the dictionary
	void add(std::string key, double value){
		double_values.insert(std::pair<std::string, double>(key, value));
	};

	/// Add a string to the dictionary
	void add(std::string key, std::string value){
		string_values.insert(std::pair<std::string,std::string>(key,value));
	};

	/// print in the form for the call to plot function in matplotlib.pyplot.plot function
	/// k1 = v1, k2 = v2, k3 = v3 ..... with vi properly quote escaped
	std::string print_call()
	{
		std::string s;

		int i = 0;
		std::string el;
		for (std::map<std::string, double>::iterator it = double_values.begin(); it != double_values.end(); it ++)
		{
			el = ", " + (*it).first + " = " + format("%0.12g",(*it).second);
			s += el;
			i++;
		}
		el = "";
		for (std::map<std::string, std::string>::iterator it = string_values.begin(); it != string_values.end(); it ++)
		{
			el = ", " + (*it).first + " = '" + (*it).second + "'" ;
			s += el;
			i++;
		}
		
		// Determine what to return (empty if there is nothing here)
		if (i > 0)
		{
			return s;
		}
		else
		{
			return std::string("");
		}
	}

};

/*! A class to contain the information for a call to plot() function
*/
class PlotCall
{
private:
	bool dict_set;
public:
	Dictionary dict;
	std::vector<double> x,y;

	/// Save internal variables for this call to plot function, including the keyword argument function as a Dictionary instance
	PlotCall(std::vector<double> x, std::vector<double> y, Dictionary *dict = NULL)
	{
		if (dict != NULL)
		{
			this->dict = *dict;
			dict_set = true;
		}
		else
		{
			dict_set = false;
		}

		this->x = x;
		this->y = y;
	};
	
	/// return the string for this call, including the preceding import matplotlib.pyplot as plt and the following plt.show()
	std::string tostring()
	{
		std::string s;
		s += "x = " + vec_to_string(x,"%0.16g") + "\n";
		s += "y = " + vec_to_string(y,"%0.16g") + "\n";
		if (&dict != NULL)
		{
			s += "plt.plot(x, y" + dict.print_call()+")\n";	
		}
		return s;
	}
};

class PyPlotter
{	
private:
	std::vector<PlotCall> PlotCalls;
public:
	void plot(std::vector<double> x, std::vector<double> y, Dictionary *dict = NULL)
	{
		PlotCalls.push_back(PlotCall(x,y,dict));
	};
	void additional_code(std::string){};

	std::string print_calls()
	{
		std::string s;
		for (std::vector<PlotCall>::iterator it = PlotCalls.begin(); it!=PlotCalls.end(); it++)
		{
			s += (*it).tostring();
		}
		return s;
	};
	void show()
	{
		std::string s = "import matplotlib\nmatplotlib.use('TkAgg')\nimport matplotlib.pyplot as plt\n";
		s += this->print_calls();
		s += "plt.savefig('AA.png')\n";
		s += "plt.show()\n";
		if (PyRun_SimpleString(s.c_str()) != 0)
		{
			FILE *fp;
			fp = fopen("errored_plot.py","w");
			fprintf(fp,"%s",s.c_str());
			fclose(fp);
			std::cout << format("plot failed.  Written log to errored_plot.py\n");
		};
	};
};

#endif
