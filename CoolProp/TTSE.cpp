#if defined(_MSC_VER)
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "CoolPropTools.h"

#if defined(__ISWINDOWS__)
#include <windows.h> // for the CreateDirectory function
#else
	#if !defined(__powerpc__)
	#include <pwd.h>
	#endif
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif



#include "CoolProp.h"
#include "CPState.h"
#include "TTSE.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>	
#include <errno.h>
#include <cerrno>
#include <float.h>

// An ugly hack to disable the timing function on PPC since the sysClkRate function not found
#if defined(__powerpc__)
#define CLOCKS_PER_SEC 1000
#endif

// The revision of the TTSE tables, only use tables with the same revision.  Increment this macro if any non-forward compatible changes are made
#define TTSEREV 6

std::string get_home_dir(void)
{
	
	// See http://stackoverflow.com/questions/2552416/how-can-i-find-the-users-home-dir-in-a-cross-platform-manner-using-c
	#if defined(__ISLINUX__)
        char *home = NULL;
		home = getenv("HOME");
		return std::string(home);
	#elif defined(__ISAPPLE__)
        char *home = NULL;
		home = getenv("HOME");
		if (home==NULL) {
		  struct passwd* pwd = getpwuid(getuid());
		  if (pwd) {
			home = pwd->pw_dir;
		  }
		}
		if (home==NULL) {
		  throw NotImplementedError("Could not detect home directory.");
		} 
//		throw NotImplementedError("This function is not defined for your platform.");
		return std::string(home);
	#elif defined(__ISWINDOWS__)
		char * pUSERPROFILE = getenv("USERPROFILE");
		if (pUSERPROFILE != NULL) {
			return std::string(pUSERPROFILE);
		} else {
			char * pHOMEDRIVE = getenv("HOMEDRIVE");
			char * pHOMEPATH = getenv("HOMEPATH");
			if (pHOMEDRIVE != NULL && pHOMEPATH != NULL) {
				return std::string(pHOMEDRIVE) + std::string(pHOMEPATH);
			} else {
				return std::string("");
			}
		}
	#else
		throw NotImplementedError("This function is not defined for your platform.");
	#endif
}

/// The inverse of the A matrix for the bicubic interpolation (http://en.wikipedia.org/wiki/Bicubic_interpolation)
const static double Ainv[16][16] = {
	{ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
	{ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
	{-3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
	{ 2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
	{ 0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0},
	{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0},
	{ 0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0},
	{ 0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0},
	{-3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0},
	{ 0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0, -2,  0, -1,  0},
	{ 9, -9, -9,  9,  6,  3, -6, -3,  6, -6,  3, -3,  4,  2,  2,  1},
	{-6,  6,  6, -6, -3, -3,  3,  3, -4,  4, -2,  2, -2, -2, -1, -1},
	{ 2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0},
	{ 0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  1,  0,  1,  0},
	{-6,  6,  6, -6, -4, -2,  4,  2, -3,  3, -3,  3, -2, -1, -2, -1},
	{ 4, -4, -4,  4,  2,  2, -2, -2,  2, -2,  2, -2,  1,  1,  1,  1}
	};

double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

void matrix_vector_product(std::vector<double> *x, std::vector<double> *y)
{
	double sum;
	for (unsigned int i = 0; i < 16; i++)
	{
		sum = 0;
		for (unsigned int j = 0; j < 16; j++)
		{
			sum += Ainv[i][j]*(*x)[j];
		}
		(*y)[i] = sum;
	}
}

TTSESinglePhaseTableClass::TTSESinglePhaseTableClass(){
	this->enable_writing_tables_to_files = true;
	this->enable_transport = true;
	SatL = NULL; 
	SatV = NULL;

	// Default to use the "normal" TTSE evaluation
	mode = TTSE_MODE_TTSE;
}

TTSESinglePhaseTableClass::TTSESinglePhaseTableClass(Fluid *pFluid) {
	this->pFluid = pFluid;
	// The default data location for the LUT
	#if defined(__ISWINDOWS__)
		this->root_path = get_home_dir()+std::string("\\.CoolProp-TTSEData\\")+pFluid->get_name();
	#else
		this->root_path = get_home_dir()+std::string("/.CoolProp-TTSEData/")+pFluid->get_name();
	#endif

	// Seed the generator for random number generation
	srand((unsigned int)time(NULL));
	this->enable_writing_tables_to_files = true;
	this->enable_transport = true;
	SatL = NULL;
	SatV = NULL;

	// Instantiate the storage vectors for bicubic interpolation; instantiated once since instantation is quite slow
	alpha_bicubic = std::vector<double>(16);
	z_bicubic = std::vector<double>(16);

	// Default to use the "normal" TTSE evaluation
	mode = TTSE_MODE_TTSE;
}

int TTSESinglePhaseTableClass::get_mode() {
	return this->mode;
}

int TTSESinglePhaseTableClass::set_mode(int mode) {
	if (mode == TTSE_MODE_TTSE || mode == TTSE_MODE_BICUBIC) {
		this->mode = mode; 
		return true;
	} else {
		return false;
	}
}

void TTSESinglePhaseTableClass::set_size_ph(unsigned int Np, unsigned int Nh) {
	this->Nh = Nh;
	this->Np = Np;

	h.resize(Nh);
	p.resize(Np);

	s.resize(Nh, std::vector<double>(Np, _HUGE));
	dsdh.resize(Nh, std::vector<double>(Np, _HUGE));
	dsdp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2sdh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2sdp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2sdhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	T.resize(Nh, std::vector<double>(Np, _HUGE));
	dTdh.resize(Nh, std::vector<double>(Np, _HUGE));
	dTdp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2Tdh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2Tdp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2Tdhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	rho.resize(Nh, std::vector<double>(Np, _HUGE));
	drhodh.resize(Nh, std::vector<double>(Np, _HUGE));
	drhodp.resize(Nh, std::vector<double>(Np, _HUGE));
	d2rhodh2.resize(Nh, std::vector<double>(Np, _HUGE));	
	d2rhodp2.resize(Nh, std::vector<double>(Np, _HUGE));
	d2rhodhdp.resize(Nh, std::vector<double>(Np, _HUGE));

	IL.resize(Np);
	IV.resize(Np);
	TL.resize(Np);
	TV.resize(Np);
	SL.resize(Np);
	SV.resize(Np);
	DL.resize(Np);
	DV.resize(Np);

	// Instantiate the cell matrices for the bicubic interpolation
	bicubic_cells = BiCubicCellsContainerClass();

	for(unsigned int i = 0; i < Nh; i++){
		bicubic_cells.cells.push_back(std::vector<BiCubicCellClass>(Np));
	}
}

void TTSESinglePhaseTableClass::set_size_Trho(unsigned int NT, unsigned int Nrho) {
	this->NT = NT;
	this->Nrho = Nrho;

	T_Trho.resize(NT);
	rho_Trho.resize(Nrho);

	s_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dsdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dsdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2sdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2sdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2sdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	p_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dpdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dpdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2pdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2pdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2pdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	h_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dhdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dhdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2hdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2hdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2hdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	k_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dkdT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dkdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2kdT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2kdrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2kdTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));

	mu_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dmudT_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	dmudrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2mudT2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));	
	d2mudrho2_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
	d2mudTdrho_Trho.resize(NT, std::vector<double>(Nrho, _HUGE));
}

void TTSESinglePhaseTableClass::nearest_good_neighbor(int *i, int *j) {
	// Left
	if (*i>0 && ValidNumber(rho[*i-1][*j]) && ValidNumber(T[*i-1][*j])){
		*i -= 1;
		return;
	}
	// Right
	else if (*i<(int)Nh-1 && ValidNumber(rho[*i+1][*j]) && ValidNumber(T[*i+1][*j])){
		*i += 1;
		return;
	}
	// Down
	else if (*j>0 && ValidNumber(rho[*i][*j-1]) && ValidNumber(T[*i][*j-1])){
		*j -= 1;
		return;
	}
	// Up
	else if (*j<(int)Np-1 && ValidNumber(rho[*i][*j+1]) && ValidNumber(T[*i][*j+1])){
		*j += 1;
		return;
	}
	// Two Left
	if (*i>1 && ValidNumber(rho[*i-2][*j]) && ValidNumber(T[*i-2][*j])){
		*i -= 2;
		return;
	}
	// Two Right
	else if (*i<(int)Nh-2 && ValidNumber(rho[*i+2][*j]) && ValidNumber(T[*i+2][*j])){
		*i += 2;
		return;
	}
	// Two Down
	else if (*j>1 && ValidNumber(rho[*i][*j-2]) && ValidNumber(T[*i][*j-2])){
		*j -= 2;
		return;
	}
	// Two Up
	else if (*j<(int)Np-2 && ValidNumber(rho[*i][*j+2]) && ValidNumber(T[*i][*j+2])){
		*j += 2;
		return;
	}

	else
	{
		throw ValueError(format("No neighbors found for %d,%d",i,j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_good_neighbor_Trho_interpolate(int *i, int *j)
{
	if (bicubic_cells.cells[*i][*j].valid_Trho){
		return;
	}
	else if (bicubic_cells.cells[*i+1][*j].valid_Trho)
	{
		*i += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j].valid_Trho)
	{
		*i -= 1; return;
	}
	else if (bicubic_cells.cells[*i][*j-1].valid_Trho)
	{
		*j -= 1; return;
	}
	else if (bicubic_cells.cells[*i][*j+1].valid_Trho)
	{
		*j += 1; return;
	}
	
	else if (bicubic_cells.cells[*i+1][*j+1].valid_Trho)
	{
		*i += 1; *j += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j-1].valid_Trho)
	{
		*i -= 1; *j -= 1; return;
	}
	else if (bicubic_cells.cells[*i+1][*j-1].valid_Trho)
	{
		*j -= 1; *i += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j+1].valid_Trho)
	{
		*j += 1; *i -= 1; return;
	}
	else
	{
		throw ValueError(format("No neighbors found for %d,%d",*i,*j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_good_neighbor_ph_interpolate(int *i, int *j)
{
	if (bicubic_cells.cells[*i][*j].valid_hp){
		return;
	}
	else if (bicubic_cells.cells[*i+1][*j].valid_hp)
	{
		*i += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j].valid_hp)
	{
		*i -= 1; return;
	}
	else if (bicubic_cells.cells[*i][*j-1].valid_hp)
	{
		*j -= 1; return;
	}
	else if (bicubic_cells.cells[*i][*j+1].valid_hp)
	{
		*j += 1; return;
	}
	
	else if (bicubic_cells.cells[*i+1][*j+1].valid_hp)
	{
		*i += 1; *j += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j-1].valid_hp)
	{
		*i -= 1; *j -= 1; return;
	}
	else if (bicubic_cells.cells[*i+1][*j-1].valid_hp)
	{
		*j -= 1; *i += 1; return;
	}
	else if (bicubic_cells.cells[*i-1][*j+1].valid_hp)
	{
		*j += 1; *i -= 1; return;
	}
	else
	{
		throw ValueError(format("No neighbors found for %d,%d",*i,*j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_good_neighbor_Trho(int *i, int *j)
{
	// Left
	if (*i>0 && ValidNumber(h_Trho[*i-1][*j]) && ValidNumber(p_Trho[*i-1][*j])){
		*i -= 1;
		return;
	}
	// Right
	else if (*i<(int)Nh-1 && ValidNumber(h_Trho[*i+1][*j]) && ValidNumber(p_Trho[*i+1][*j])){
		*i += 1;
		return;
	}
	// Down
	else if (*j>0 && ValidNumber(h_Trho[*i][*j-1]) && ValidNumber(p_Trho[*i][*j-1])){
		*j -= 1;
		return;
	}
	// Up
	else if (*j<(int)Np-1 && ValidNumber(h_Trho[*i][*j+1]) && ValidNumber(p_Trho[*i][*j+1])){
		*j += 1;
		return;
	}
	else
	{
		throw ValueError(format("No neighbors found for %d,%d",i,j));
		return;
	}
}

void TTSESinglePhaseTableClass::nearest_neighbor_ph(int i, int j, double *T0, double *rho0)
{
	// Left
	if (i>0 && ValidNumber(rho[i-1][j]) && ValidNumber(T[i-1][j])){
		*T0 = T[i-1][j];
		*rho0 = rho[i-1][j];
		return;
	}
	// Right
	else if (i<(int)Nh-1 && ValidNumber(rho[i+1][j]) && ValidNumber(T[i+1][j])){
		*T0 = T[i+1][j];
		*rho0 = rho[i+1][j];
		return;
	}
	// Down
	else if (j>0 && ValidNumber(rho[i][j-1]) && ValidNumber(T[i][j-1])){
		*T0 = T[i][j-1];
		*rho0 = rho[i][j-1];
		return;
	}
	// Up
	else if (j<(int)Np-1 && ValidNumber(rho[i][j+1]) && ValidNumber(T[i][j+1])){
		*T0 = T[i][j+1];
		*rho0 = rho[i][j+1];
		return;
	}
	else
	{
		*T0 = -1;
		*rho0 = -1;
		return;
	}
}
std::string join(std::vector<std::string> strings, char delim)
{
	std::string output = strings[0];
	for (unsigned int i = 1; i < strings.size(); i++)
	{
		output += format("%c%s",delim,strings[i].c_str());
	}
	return output;
}

void TTSESinglePhaseTableClass::matrix_to_file(std::string fName, std::vector< std::vector<double> > *A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"wb");
	if (pFile != NULL)
	{
		for (unsigned int i = 0; i < Nh; i++)
		{
			fwrite(&(((*A)[i])[0]), sizeof(double),Np,pFile);
		}
		fclose(pFile);
	}
}

void TTSESinglePhaseTableClass::matrix_from_file(std::string fName, std::vector<std::vector<double> > *A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"rb");
	if (pFile != NULL)
	{
		for (unsigned int i = 0; i < Nh; i++)
		{
			/*std::vector<double> row(Nh,0);
			fread(&(row[0]), sizeof(double), Np, pFile);
			double tgregt = 1;*/
			fread(&(((*A)[i])[0]), sizeof(double),Np,pFile);
		}
		fclose(pFile);
	}
}

void TTSESinglePhaseTableClass::vector_to_file(std::string fName, std::vector<double>*A)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"wb");
	if (pFile != NULL)
	{
		fwrite((const char *)&(*A).front(), sizeof(double),(*A).size(),pFile);
		fclose(pFile);
	}
}
void TTSESinglePhaseTableClass::vector_from_file(std::string fName, int N, std::vector<double> *vec)
{
	FILE *pFile;
	pFile = fopen(fName.c_str(),"rb");
	if (pFile != NULL)
	{
		fread(&((*vec)[0]), sizeof(double), N, pFile);
		fclose(pFile);
	}
}
bool pathExists(std::string path)
{
	#if defined(__ISWINDOWS__) // Defined for 32-bit and 64-bit windows
		struct _stat buf;
		// Get data associated with path using the windows libraries, 
		// and if you can (result == 0), the path exists
		if ( _stat( path.c_str(), &buf) == 0)
			return true;
		else
			return false;
	#elif defined(__ISLINUX__) || defined(__ISAPPLE__)
		struct stat st;
//		if(stat(path.c_str(),&st) == 0) {
//			if(st.st_mode & S_IFDIR != 0) return true;
//		    if(st.st_mode & S_IFREG != 0) return true;
//		    return false;
//		} else {
//			return false;
//		}
		if(lstat(path.c_str(),&st) == 0) {
			if(S_ISDIR(st.st_mode)) return true;
			if(S_ISREG(st.st_mode)) return true;
			return false;
		} else {
			return false;
		}

	#else
		throw NotImplementedError("This function is not defined for your platform.");
	#endif
}
//bool fileExists(const char *fileName)
//{
//	std::ifstream infile(fileName);
//    return infile.good();
//}
void make_dirs(std::string file_path)
{
	std::vector<std::string> pathsplit = strsplit(file_path,'/');
	std::string path = pathsplit[0];
	if (pathsplit.size()>0)
	{
		for (unsigned int i = 0; i < pathsplit.size(); i++)
		{
			if (!pathExists(path))
			{
			#ifdef _WIN32
				#if defined(_UNICODE)
					CreateDirectoryA((LPCSTR)path.c_str(),NULL);
				#else
					CreateDirectory((LPCSTR)path.c_str(),NULL);
				#endif
			#else 
				#if defined(__powerpc__)
				#else
					mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				#endif
			#endif
			}
			if (i < pathsplit.size()-1)
				path += format("/%s",pathsplit[i+1].c_str());
		}
	}
	else
	{
		throw ValueError(format("Could not make path [%s]",file_path.c_str()));
	}
}

bool TTSESinglePhaseTableClass::read_all_from_file(std::string root_path)
{
	std::string Fluid;
	double hmin,hmax,pmin,pmax;
	int Np, Nh, TTSERev;

	// Replace any '\' with '/' in the path
	for (unsigned int i = 0; i<root_path.length(); i++)
	{
		if (root_path[i] == '\\') root_path[i] = '/';
	}
	// If it ends with a '/', remove it for now
	if (root_path[root_path.size()-1] == '/')
		root_path = std::string(root_path,0,root_path.size()-1);

	if (!pathExists(root_path)) return false;

	// Append a '/'
	root_path += format("/");

	if (!pathExists(root_path+format("Info_ph.txt"))) return false;

	// Parse the info file to find the dimensions and boundaries and check that they are the same
	std::string info = get_file_contents((root_path+format("Info_ph.txt")).c_str());

	// Split into lines
#if defined(__ISWINDOWS__)
	std::vector<std::string> lines = strsplit(info,'\r');
#else
	std::vector<std::string> lines = strsplit(info,'\n');
#endif

	for (unsigned int i = 0; i< lines.size(); i++) {

		// Split at ':'
		std::vector<std::string> line = strsplit(lines[i],':');
		if (line.size() == 1) // No : found
			continue;

		if (line[0].find("Fluid")!=  std::string::npos){
			Fluid = line[1];}
		else if (line[0].find("hmin")!=  std::string::npos) {	
			hmin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("hmax")!=  std::string::npos) {	
			hmax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("pmin")!=  std::string::npos) {	
			pmin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("pmax")!=  std::string::npos) {	
			pmax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("Np")!=  std::string::npos) {	
			Np = (int)strtol(line[1].c_str(),NULL,0);}
		else if (line[0].find("Nh")!=  std::string::npos) {	
			Nh = (int)strtol(line[1].c_str(),NULL,0);}
		else if (line[0].find("Tmin")!=  std::string::npos) {
			Tmin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("Tmax")!=  std::string::npos) {	
			Tmax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("rhomin")!=  std::string::npos) {	
			rhomin = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("rhomax")!=  std::string::npos) {	
			rhomax = strtod(line[1].c_str(),NULL);}
		else if (line[0].find("NT")!=  std::string::npos) {	
			NT = (int)strtol(line[1].c_str(),NULL,0);}
		else if (line[0].find("Nrho")!=  std::string::npos) {	
			Nrho = (int)strtol(line[1].c_str(),NULL,0);}
		else if (line[0].find("TTSERev")!=  std::string::npos) {
			TTSERev = (int)strtol(line[1].c_str(),NULL,0);}
	}

	// Didn't work since at least one of the parameters was different
	// so we need to build the tables again
	if (!(!Fluid.compare(pFluid->get_name())
		  && Nh == (int)this->Nh
		  && Np == (int)this->Np
		  && fabs(pmin - this->pmin)<10*DBL_EPSILON
		  && fabs(pmax - this->pmax)<10*DBL_EPSILON
		  && fabs(hmin - this->hmin)<10*DBL_EPSILON
		  && fabs(hmax - this->hmax)<10*DBL_EPSILON
		  && NT == this->NT
		  && Nrho == this->Nrho
		  && fabs(Tmin - this->Tmin)<10*DBL_EPSILON
		  && fabs(Tmax - this->Tmax)<10*DBL_EPSILON
		  && fabs(rhomin - this->rhomin)<10*DBL_EPSILON
		  && fabs(rhomax - this->rhomax)<10*DBL_EPSILON
		  && TTSERev == TTSEREV
		)) return false;

	// Read all the data from the binary files
	vector_from_file(root_path + std::string("p_ph.ttse"),Np,&p);
	vector_from_file(root_path + std::string("h_ph.ttse"),Nh,&h);
	matrix_from_file(root_path + std::string("T_ph.ttse"),&T);
	matrix_from_file(root_path + std::string("dTdh_ph.ttse"),&dTdh);
	matrix_from_file(root_path + std::string("dTdp_ph.ttse"),&dTdp);
	matrix_from_file(root_path + std::string("d2Tdh2_ph.ttse"),&d2Tdh2);
	matrix_from_file(root_path + std::string("d2Tdp2_ph.ttse"),&d2Tdp2);
	matrix_from_file(root_path + std::string("d2Tdhdp_ph.ttse"),&d2Tdhdp);
	matrix_from_file(root_path + std::string("s_ph.ttse"),&s);
	matrix_from_file(root_path + std::string("dsdh_ph.ttse"),&dsdh);
	matrix_from_file(root_path + std::string("dsdp_ph.ttse"),&dsdp);
	matrix_from_file(root_path + std::string("d2sdh2_ph.ttse"),&d2sdh2);
	matrix_from_file(root_path + std::string("d2sdp2_ph.ttse"),&d2sdp2);
	matrix_from_file(root_path + std::string("d2sdhdp_ph.ttse"),&d2sdhdp);
	matrix_from_file(root_path + std::string("rho_ph.ttse"),&rho);
	matrix_from_file(root_path + std::string("drhodh_ph.ttse"),&drhodh);
	matrix_from_file(root_path + std::string("drhodp_ph.ttse"),&drhodp);
	matrix_from_file(root_path + std::string("d2rhodh2_ph.ttse"),&d2rhodh2);
	matrix_from_file(root_path + std::string("d2rhodp2_ph.ttse"),&d2rhodp2);
	matrix_from_file(root_path + std::string("d2rhoTdhdp_ph.ttse"),&d2rhodhdp);

	vector_from_file(root_path + std::string("T_Trho.ttse"),NT,&T_Trho);
	vector_from_file(root_path + std::string("rho_Trho.ttse"),Nrho,&rho_Trho);
	matrix_from_file(root_path + std::string("p_Trho.ttse"),&p_Trho);
	matrix_from_file(root_path + std::string("dpdT_Trho.ttse"),&dpdT_Trho);
	matrix_from_file(root_path + std::string("dpdrho_Trho.ttse"),&dpdrho_Trho);
	matrix_from_file(root_path + std::string("d2pdT2_Trho.ttse"),&d2pdT2_Trho);
	matrix_from_file(root_path + std::string("d2pdrho2_Trho.ttse"),&d2pdrho2_Trho);
	matrix_from_file(root_path + std::string("d2pdTdrho_Trho.ttse"),&d2pdTdrho_Trho);
	matrix_from_file(root_path + std::string("s_Trho.ttse"),&s_Trho);
	matrix_from_file(root_path + std::string("dsdT_Trho.ttse"),&dsdT_Trho);
	matrix_from_file(root_path + std::string("dsdrho_Trho.ttse"),&dsdrho_Trho);
	matrix_from_file(root_path + std::string("d2sdT2_Trho.ttse"),&d2sdT2_Trho);
	matrix_from_file(root_path + std::string("d2sdrho2_Trho.ttse"),&d2sdrho2_Trho);
	matrix_from_file(root_path + std::string("d2sdTdrho_Trho.ttse"),&d2sdTdrho_Trho);
	matrix_from_file(root_path + std::string("h_Trho.ttse"),&h_Trho);
	matrix_from_file(root_path + std::string("dhdT_Trho.ttse"),&dhdT_Trho);
	matrix_from_file(root_path + std::string("dhdrho_Trho.ttse"),&dhdrho_Trho);
	matrix_from_file(root_path + std::string("d2hdT2_Trho.ttse"),&d2hdT2_Trho);
	matrix_from_file(root_path + std::string("d2hdrho2_Trho.ttse"),&d2hdrho2_Trho);
	matrix_from_file(root_path + std::string("d2hdTdrho_Trho.ttse"),&d2hdTdrho_Trho);
	matrix_from_file(root_path + std::string("mu_Trho.ttse"),&mu_Trho);
	matrix_from_file(root_path + std::string("dmudT_Trho.ttse"),&dmudT_Trho);
	matrix_from_file(root_path + std::string("dmudrho_Trho.ttse"),&dmudrho_Trho);
	matrix_from_file(root_path + std::string("d2mudT2_Trho.ttse"),&d2mudT2_Trho);
	matrix_from_file(root_path + std::string("d2mudrho2_Trho.ttse"),&d2mudrho2_Trho);
	matrix_from_file(root_path + std::string("d2mudTdrho_Trho.ttse"),&d2mudTdrho_Trho);
	matrix_from_file(root_path + std::string("k_Trho.ttse"),&k_Trho);
	matrix_from_file(root_path + std::string("dkdT_Trho.ttse"),&dkdT_Trho);
	matrix_from_file(root_path + std::string("dkdrho_Trho.ttse"),&dkdrho_Trho);
	matrix_from_file(root_path + std::string("d2kdT2_Trho.ttse"),&d2kdT2_Trho);
	matrix_from_file(root_path + std::string("d2kdrho2_Trho.ttse"),&d2kdrho2_Trho);
	matrix_from_file(root_path + std::string("d2kdTdrho_Trho.ttse"),&d2kdTdrho_Trho);

	this->pratio = pow(pmax/pmin,1/((double)Np-1));
	this->logpratio = log(pratio); // For speed since log() is a slow function
	this->logpmin = log(pmin);
	this->rhoratio = pow(rhomax/rhomin,1/((double)Nrho-1));
	this->logrhoratio = log(rhoratio); // For speed since log() is a slow function
	this->logrhomin = log(rhomin);
	this->jpcrit_floor = (int)floor((log(pFluid->reduce.p.Pa)-logpmin)/logpratio);
	this->jpcrit_ceil = (int)ceil((log(pFluid->reduce.p.Pa)-logpmin)/logpratio);
	
	update_saturation_boundary_indices();

	update_cell_validity();

	return true;
}
void TTSESinglePhaseTableClass::write_all_to_file(std::string root_path)
{
	// Replace any '\' with '/' in the path
	for (unsigned int i = 0; i<root_path.length(); i++)
	{
		if (root_path[i] == '\\') root_path[i] = '/';
	}
	// If it ends with a '/', remove it for now
	if (root_path[root_path.size()-1] == '/')
		root_path = std::string(root_path,0,root_path.size()-1);

	if (!pathExists(root_path))
		make_dirs(root_path);

	// Append a '/'
	root_path += format("/");

	std::string header = std::string("Data for the TTSE method\nDO NOT CHANGE ANY OF THESE PARAMETERS FOR ANY REASON!\n\n");
		
	header += format("TTSERev:%d\nFluid:%s\npmin:%23.19g\npmax:%23.19g\nNp:%25d\nhmin:%23.19g\nhmax:%23.19g\nNh:%25d\nTmin:%23.19g\nTmax:%23.19g\nNT:%25d\nrhomin:%23.19g\nrhomax:%23.19g\nNrho:%25d\n",TTSEREV,pFluid->get_name().c_str(),pmin,pmax,Np,hmin,hmax,Nh,Tmin,Tmax,NT,rhomin,rhomax,Nrho);
	
	clock_t t1,t2;
	t1 = clock();

	// Write the header information to a text file
	FILE *fp;
	fp = fopen((root_path+std::string("Info_ph.txt")).c_str(),"w");
	fprintf(fp,"%s",header.c_str());
	fclose(fp);

	// Write each of these files in binary mode
	vector_to_file(root_path + std::string("p_ph.ttse"),&p);
	vector_to_file(root_path + std::string("h_ph.ttse"),&h);
	matrix_to_file(root_path + std::string("T_ph.ttse"),&T);
	matrix_to_file(root_path + std::string("dTdh_ph.ttse"),&dTdh);
	matrix_to_file(root_path + std::string("dTdp_ph.ttse"),&dTdp);
	matrix_to_file(root_path + std::string("d2Tdh2_ph.ttse"),&d2Tdh2);
	matrix_to_file(root_path + std::string("d2Tdp2_ph.ttse"),&d2Tdp2);
	matrix_to_file(root_path + std::string("d2Tdhdp_ph.ttse"),&d2Tdhdp);
	matrix_to_file(root_path + std::string("s_ph.ttse"),&s);
	matrix_to_file(root_path + std::string("dsdh_ph.ttse"),&dsdh);
	matrix_to_file(root_path + std::string("dsdp_ph.ttse"),&dsdp);
	matrix_to_file(root_path + std::string("d2sdh2_ph.ttse"),&d2sdh2);
	matrix_to_file(root_path + std::string("d2sdp2_ph.ttse"),&d2sdp2);
	matrix_to_file(root_path + std::string("d2sdhdp_ph.ttse"),&d2sdhdp);
	matrix_to_file(root_path + std::string("rho_ph.ttse"),&rho);
	matrix_to_file(root_path + std::string("drhodh_ph.ttse"),&drhodh);
	matrix_to_file(root_path + std::string("drhodp_ph.ttse"),&drhodp);
	matrix_to_file(root_path + std::string("d2rhodh2_ph.ttse"),&d2rhodh2);
	matrix_to_file(root_path + std::string("d2rhodp2_ph.ttse"),&d2rhodp2);
	matrix_to_file(root_path + std::string("d2rhoTdhdp_ph.ttse"),&d2rhodhdp);
	
	vector_to_file(root_path + std::string("T_Trho.ttse"),&T_Trho);
	vector_to_file(root_path + std::string("rho_Trho.ttse"),&rho_Trho);
	matrix_to_file(root_path + std::string("p_Trho.ttse"),&p_Trho);
	matrix_to_file(root_path + std::string("dpdT_Trho.ttse"),&dpdT_Trho);
	matrix_to_file(root_path + std::string("dpdrho_Trho.ttse"),&dpdrho_Trho);
	matrix_to_file(root_path + std::string("d2pdT2_Trho.ttse"),&d2pdT2_Trho);
	matrix_to_file(root_path + std::string("d2pdrho2_Trho.ttse"),&d2pdrho2_Trho);
	matrix_to_file(root_path + std::string("d2pdTdrho_Trho.ttse"),&d2pdTdrho_Trho);
	matrix_to_file(root_path + std::string("s_Trho.ttse"),&s_Trho);
	matrix_to_file(root_path + std::string("dsdT_Trho.ttse"),&dsdT_Trho);
	matrix_to_file(root_path + std::string("dsdrho_Trho.ttse"),&dsdrho_Trho);
	matrix_to_file(root_path + std::string("d2sdT2_Trho.ttse"),&d2sdT2_Trho);
	matrix_to_file(root_path + std::string("d2sdrho2_Trho.ttse"),&d2sdrho2_Trho);
	matrix_to_file(root_path + std::string("d2sdTdrho_Trho.ttse"),&d2sdTdrho_Trho);
	matrix_to_file(root_path + std::string("h_Trho.ttse"),&h_Trho);
	matrix_to_file(root_path + std::string("dhdT_Trho.ttse"),&dhdT_Trho);
	matrix_to_file(root_path + std::string("dhdrho_Trho.ttse"),&dhdrho_Trho);
	matrix_to_file(root_path + std::string("d2hdT2_Trho.ttse"),&d2hdT2_Trho);
	matrix_to_file(root_path + std::string("d2hdrho2_Trho.ttse"),&d2hdrho2_Trho);
	matrix_to_file(root_path + std::string("d2hdTdrho_Trho.ttse"),&d2hdTdrho_Trho);
	matrix_to_file(root_path + std::string("k_Trho.ttse"),&k_Trho);
	matrix_to_file(root_path + std::string("dkdT_Trho.ttse"),&dkdT_Trho);
	matrix_to_file(root_path + std::string("dkdrho_Trho.ttse"),&dkdrho_Trho);
	matrix_to_file(root_path + std::string("d2kdT2_Trho.ttse"),&d2kdT2_Trho);
	matrix_to_file(root_path + std::string("d2kdrho2_Trho.ttse"),&d2kdrho2_Trho);
	matrix_to_file(root_path + std::string("d2kdTdrho_Trho.ttse"),&d2kdTdrho_Trho);
	matrix_to_file(root_path + std::string("mu_Trho.ttse"),&mu_Trho);
	matrix_to_file(root_path + std::string("dmudT_Trho.ttse"),&dmudT_Trho);
	matrix_to_file(root_path + std::string("dmudrho_Trho.ttse"),&dmudrho_Trho);
	matrix_to_file(root_path + std::string("d2mudT2_Trho.ttse"),&d2mudT2_Trho);
	matrix_to_file(root_path + std::string("d2mudrho2_Trho.ttse"),&d2mudrho2_Trho);
	matrix_to_file(root_path + std::string("d2mudTdrho_Trho.ttse"),&d2mudTdrho_Trho);

	t2 = clock();
	std::cout << "write time: " << (double)(t2-t1)/CLOCKS_PER_SEC << std::endl;
}

double TTSESinglePhaseTableClass::build_ph(double hmin, double hmax, double pmin, double pmax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV)
{
	bool SinglePhase = false;

	this->hmin = hmin;
	this->hmax = hmax;
	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);

	CoolPropStateClassSI CPS(pFluid);

	long iFluid = get_Fluid_index(CPS.pFluid->get_name());

	double dh = (hmax - hmin)/(Nh - 1);
	pratio = pow(pmax/pmin,1/((double)Np-1));
	logpratio = log(pratio);

	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i<Nh; i++)
	{
		double hval = hmin + i*dh;
		h[i] = hval;
		for (unsigned int j = 0; j<Np; j++)
		{
			double pval = pmin*pow(pratio,(int)j);
			p[j] = pval;
			
			// Check whether the point is single phase
			// If pressure between ptriple point and pcrit, might be two-phase or single phase, otherwise definitely single phase
			if (pval <= pFluid->reduce.p.Pa && pval >= pFluid->params.ptriple)
			{
				if (SatL == NULL || SatV == NULL){
					// Not using TTSE method, use saturation (slow...)
					CPS.update(iP,pval,iQ,0.5);
					SinglePhase = (hval < CPS.hL() || hval > CPS.hV());
				}
				else{
					// Using the TTSE method, nice and fast
					SinglePhase = (hval < SatL->evaluate(iH,pval)  || hval > SatV->evaluate(iH,pval));
				}
			}
			else
			{
				SinglePhase = true;
			}
			
			// If enthalpy is outside the saturation region, continue and do the calculation as a function of p,h
			if (SinglePhase)
			{
				try
				{
					double T0=-1,rho0=-1,T,rho,rhoL,rhoV,TsatL,TsatV;

					// Find a good point around this point that is single-phase if any of its neighbors have been calculated
					nearest_neighbor_ph(i,j,&T0,&rho0);

					// If good T,rho was returned, use it as a guess value to calculate T,rho from p,h more quickly
					if (T0 > 0 && rho0 > 0)
					{
						// Get T,rho from p,h using our guess value
						CPS.pFluid->temperature_ph(pval,hval,T,rho,rhoL,rhoV,TsatL,TsatV,T0,rho0);
						CPS.flag_SinglePhase = true;
						CPS.update(iT, T, iD, rho);
					}
					else
					{
						// Probably here you are close to the saturation boundary since TTSE says it is single phase, but no neighbors to lean on
						CPS.flag_SinglePhase = true;

						double hsatLTTSE,hsatVTTSE;
						if (SatL == NULL || SatV == NULL){
							// Not using TTSE method, use saturation (slow...)
							CPS.update(iP,pval,iQ,0.5);
							hsatLTTSE = CPS.hL();
							hsatVTTSE = CPS.hV();
						}	
						else
						{
							hsatLTTSE = SatL->evaluate(iH,pval);
							hsatVTTSE = SatV->evaluate(iH,pval);
						}

						if (fabs(hval-hsatLTTSE)<10)
						{
							// Close to the saturated liquid boundary
							T0 = SatL->evaluate(iT,pval);
							rho0 = SatL->evaluate(iD,pval);
							CPS.pFluid->temperature_ph(pval,hval,T,rho,rhoL,rhoV,TsatL,TsatV,T0,rho0);
							CPS.flag_SinglePhase = true;
							CPS.update(iT, T, iD, rho);
						}
						else if (fabs(hval-hsatVTTSE)<10)
						{
							// Close to the saturated vapor boundary
							T0 = SatV->evaluate(iT,pval);
							rho0 = SatV->evaluate(iD,pval);
							CPS.pFluid->temperature_ph(pval,hval,T,rho,rhoL,rhoV,TsatL,TsatV,T0,rho0);
							CPS.flag_SinglePhase = true;
							CPS.update(iT, T, iD, rho);
						}
						else
						{
							CPS.update(iP, pval, iH, hval);
						}
						//std::cout << format("%d %d %g %g\n",i,j,pval,hval);
					}
					
					T = CPS.T();
					rho = CPS.rho();
					double cp = CPS.cp();

					double A = CPS.dpdT_constrho()*CPS.dhdrho_constT()-CPS.dpdrho_constT()*CPS.dhdT_constrho();

					this->T[i][j] = T;
					dTdh[i][j] = 1/cp;
					dTdp[i][j] = 1/A*CPS.dhdrho_constT();
					this->rho[i][j] = rho;
					drhodh[i][j] = 1/A*CPS.dpdT_constrho();
					drhodp[i][j] = -1/A*CPS.dhdT_constrho();
					s[i][j] = CPS.s();
					dsdh[i][j] = 1/T;
					dsdp[i][j] = -1/(T*rho);
					
					// Matrices for second derivatives of entropy as a function of pressure and enthalpy
					d2sdh2[i][j] = -1/(T*T)*dTdh[i][j];
					d2sdhdp[i][j] = -1/(T*T)*dTdp[i][j];
					d2sdp2[i][j] = 1/(T*T*rho)*dTdp[i][j]+1/(T*rho*rho)*drhodp[i][j];

					// These are common terms needed for a range of terms for T(h,p) as well as rho(h,p)
					double dAdT_constrho = CPS.d2pdT2_constrho()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrhodT()-CPS.d2pdrhodT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdT2_constrho();
					double dAdrho_constT = CPS.d2pdrhodT()*CPS.dhdrho_constT()+CPS.dpdT_constrho()*CPS.d2hdrho2_constT()-CPS.d2pdrho2_constT()*CPS.dhdT_constrho()-CPS.dpdrho_constT()*CPS.d2hdrhodT();

					//Matrices for temperature as a function of pressure and enthalpy
					double ddT_dTdp_h_constrho = 1/A*CPS.d2hdrhodT()-1/(A*A)*dAdT_constrho*CPS.dhdrho_constT(); //[check]
					double ddrho_dTdp_h_constT = 1/A*CPS.d2hdrho2_constT()-1/(A*A)*dAdrho_constT*CPS.dhdrho_constT(); //[check
					double ddT_dTdh_p_constrho = -1/(cp*cp)*(CPS.d2hdT2_constrho()-CPS.dhdp_constT()*CPS.d2pdT2_constrho()+CPS.d2hdrhodT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrhodT());
					double ddrho_dTdh_p_constT = -1/(cp*cp)*(CPS.d2hdrhodT()-CPS.dhdp_constT()*CPS.d2pdrhodT()+CPS.d2hdrho2_constT()*CPS.drhodT_constp()-CPS.dhdp_constT()*CPS.drhodT_constp()*CPS.d2pdrho2_constT());
					
					d2Tdh2[i][j]  = ddT_dTdh_p_constrho/CPS.dhdT_constp()+ddrho_dTdh_p_constT/CPS.dhdrho_constp();
					d2Tdhdp[i][j] = ddT_dTdp_h_constrho/CPS.dhdT_constp()+ddrho_dTdp_h_constT/CPS.dhdrho_constp();
					d2Tdp2[i][j]  = ddT_dTdp_h_constrho/CPS.dpdT_consth()+ddrho_dTdp_h_constT/CPS.dpdrho_consth();

					//// Matrices for density as a function of pressure and enthalpy
					double ddT_drhodp_h_constrho = -1/A*CPS.d2hdT2_constrho()+1/(A*A)*dAdT_constrho*CPS.dhdT_constrho();
					double ddrho_drhodp_h_constT = -1/A*CPS.d2hdrhodT()+1/(A*A)*dAdrho_constT*CPS.dhdT_constrho();
					double ddT_drhodh_p_constrho = 1/A*CPS.d2pdT2_constrho()-1/(A*A)*dAdT_constrho*CPS.dpdT_constrho();
					double ddrho_drhodh_p_constT = 1/A*CPS.d2pdrhodT()-1/(A*A)*dAdrho_constT*CPS.dpdT_constrho();
					
					d2rhodh2[i][j]  = ddT_drhodh_p_constrho/CPS.dhdT_constp()+ddrho_drhodh_p_constT/CPS.dhdrho_constp();
					d2rhodhdp[i][j] = ddT_drhodp_h_constrho/CPS.dhdT_constp()+ddrho_drhodp_h_constT/CPS.dhdrho_constp();
					d2rhodp2[i][j]  = ddT_drhodp_h_constrho/CPS.dpdT_consth()+ddrho_drhodp_h_constT/CPS.dpdrho_consth();
				}
				catch(std::exception &)
				{
					s[i][j] = _HUGE;
					dsdh[i][j] = _HUGE;
					dsdp[i][j] = _HUGE;
					d2sdh2[i][j] = _HUGE;
					d2sdhdp[i][j] = _HUGE;
					d2sdp2[i][j] = _HUGE;

					this->T[i][j] = _HUGE;
					dTdh[i][j] = _HUGE;
					dTdp[i][j] = _HUGE;
					d2Tdh2[i][j]  = _HUGE;
					d2Tdhdp[i][j] = _HUGE;
					d2Tdp2[i][j]  = _HUGE;

					this->rho[i][j] = _HUGE;
					drhodh[i][j] = _HUGE;
					drhodp[i][j] = _HUGE;
					d2rhodh2[i][j]  = _HUGE;
					d2rhodhdp[i][j] = _HUGE;
					d2rhodp2[i][j]  = _HUGE;
				}
			}
			else
			{
				s[i][j] = _HUGE;
				dsdh[i][j] = _HUGE;
				dsdp[i][j] = _HUGE;
				d2sdh2[i][j] = _HUGE;
				d2sdhdp[i][j] = _HUGE;
				d2sdp2[i][j] = _HUGE;

				this->T[i][j] = _HUGE;
				dTdh[i][j] = _HUGE;
				dTdp[i][j] = _HUGE;
				d2Tdh2[i][j]  = _HUGE;
				d2Tdhdp[i][j] = _HUGE;
				d2Tdp2[i][j]  = _HUGE;

				this->rho[i][j] = _HUGE;
				drhodh[i][j] = _HUGE;
				drhodp[i][j] = _HUGE;
				d2rhodh2[i][j]  = _HUGE;
				d2rhodhdp[i][j] = _HUGE;
				d2rhodp2[i][j]  = _HUGE;


			}
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout << elap << " to build single phase table with p,h" << std::endl;

	// Update the boundaries of the points within the single-phase regions
	update_saturation_boundary_indices();
	
	// Update the cell validity for all cells
	update_cell_validity();

	return elap;
}
double TTSESinglePhaseTableClass::build_Trho(double Tmin, double Tmax, double rhomin, double rhomax, TTSETwoPhaseTableClass *SatL, TTSETwoPhaseTableClass *SatV)
{
	bool SinglePhase = false;

	if (Tmin < 0 && Tmax < 0 && rhomin < 0 && rhomax < 0)
	{
		rhomin = 9e9;
		rhomax = 0;
		Tmin = 9e9;
		Tmax = 0;
		// Use single-phase table to figure out the range for T,rho
		for (unsigned int i = 0; i<Nh; i++)
		{
			for (unsigned int j = 0; j<Np; j++)
			{
				if (ValidNumber(rho[i][j]) && rho[i][j] > rhomax){
					rhomax = rho[i][j];
				}
				if (ValidNumber(rho[i][j]) && rho[i][j] < rhomin){
					rhomin = rho[i][j];
				}
				if (ValidNumber(T[i][j]) && T[i][j] > Tmax){
					Tmax = T[i][j];
				}
				if (ValidNumber(T[i][j]) && T[i][j] < Tmin){
					Tmin = T[i][j];
				}
			}
		}
	}
	if (Tmin < pFluid->limits.Tmin){ Tmin = pFluid->limits.Tmin; }
	this->Tmin = Tmin;
	this->Tmax = Tmax;
	this->rhomin = rhomin;
	this->rhomax = rhomax;
	this->logrhomin = log(rhomin);

	rhoratio = pow(rhomax/rhomin,1/((double)Nrho-1));
	logrhoratio = log(rhoratio);

	CoolPropStateClassSI CPS(pFluid);

	long iFluid = get_Fluid_index(pFluid->get_name());

	double dT = (Tmax - Tmin)/((double)NT - 1);

	clock_t t1,t2;
	t1 = clock();

	for (unsigned int i = 0; i<NT; i++)
	{
		double Tval = Tmin + i*dT;
		T_Trho[i] = Tval;
		for (unsigned int j = 0; j<Nrho; j++)
		{
			double rhoval = rhomin*pow(rhoratio,(int)j);
			rho_Trho[j] = rhoval;
			
			// Check whether the point is single phase
			// If pressure between Ttriple point and Tcrit, might be two-phase or single phase, otherwise definitely single phase
			if (Tval <= pFluid->crit.T && Tval >= pFluid->params.Ttriple)
			{
				if (SatL == NULL || SatV == NULL){
					// Not using TTSE method, use saturation (slow...)
					CPS.update(iT,Tval,iQ,0.5);
					SinglePhase = (rhoval < CPS.rhoV() || rhoval > CPS.rhoL());
				}
				else{
					// Using the TTSE method, nice and fast
					double psatV = SatV->evaluate_T(Tval);
					double psatL = SatL->evaluate_T(Tval);
					SinglePhase = (rhoval < SatV->evaluate(iD,psatV)  || rhoval > SatL->evaluate(iD,psatL));
				}
			}
			else
			{
				SinglePhase = true;
			}
			
			// If enthalpy is outside the saturation region, continue and do the calculation as a function of T,rho
			if (SinglePhase)
			{
				CPS.update(iT,Tval,iD,rhoval);

				s_Trho[i][j] = CPS.s();
				dsdT_Trho[i][j] = CPS.dsdT_constrho();
				dsdrho_Trho[i][j] = CPS.dsdrho_constT();
				d2sdT2_Trho[i][j] = CPS.d2sdT2_constrho();
				d2sdTdrho_Trho[i][j] = CPS.d2sdrhodT();
				d2sdrho2_Trho[i][j] = CPS.d2sdrho2_constT();

				h_Trho[i][j] = CPS.h();
				dhdT_Trho[i][j] = CPS.dhdT_constrho();
				dhdrho_Trho[i][j] = CPS.dhdrho_constT();
				d2hdT2_Trho[i][j] = CPS.d2hdT2_constrho();
				d2hdTdrho_Trho[i][j] = CPS.d2hdrhodT();
				d2hdrho2_Trho[i][j] = CPS.d2hdrho2_constT();

				p_Trho[i][j] = CPS.p();
				dpdT_Trho[i][j] = CPS.dpdT_constrho();
				dpdrho_Trho[i][j] = CPS.dpdrho_constT();
				d2pdT2_Trho[i][j] = CPS.d2pdT2_constrho();
				d2pdTdrho_Trho[i][j] = CPS.d2pdrhodT();
				d2pdrho2_Trho[i][j] = CPS.d2pdrho2_constT();

				
				/// Transport properties
				///
				// Using second-order centered finite differences to calculate the transport property derivatives
				double deltaT = 1e-3, deltarho = 1e-4;

				if (deltarho > rhoval)
					deltarho = rhoval/100;

				// The viscosity values
				double muval =             IPropsSI(iV,iT,Tval,       iD,rhoval,         iFluid);
				double muplusrho =         IPropsSI(iV,iT,Tval,       iD,rhoval+deltarho,iFluid);
				double muminusrho =        IPropsSI(iV,iT,Tval,       iD,rhoval-deltarho,iFluid);
				double muplusT =           IPropsSI(iV,iT,Tval+deltaT,iD,rhoval,         iFluid);
				double muminusT =          IPropsSI(iV,iT,Tval-deltaT,iD,rhoval,         iFluid);
				double muplusT_plusrho =   IPropsSI(iV,iT,Tval+deltaT,iD,rhoval+deltarho,iFluid);
				double muplusT_minusrho =  IPropsSI(iV,iT,Tval+deltaT,iD,rhoval-deltarho,iFluid);
				double muminusT_plusrho =  IPropsSI(iV,iT,Tval-deltaT,iD,rhoval+deltarho,iFluid);
				double muminusT_minusrho = IPropsSI(iV,iT,Tval-deltaT,iD,rhoval-deltarho,iFluid);

				mu_Trho[i][j] = muval;
				dmudT_Trho[i][j] = (-muminusT + muplusT)/(2*deltaT);
				dmudrho_Trho[i][j] = (-muminusrho + muplusrho)/(2*deltarho);
				d2mudT2_Trho[i][j] = (muminusT - 2*muval + muplusT)/(deltaT*deltaT);
				d2mudrho2_Trho[i][j] = (muminusrho - 2*muval + muplusrho)/(deltarho*deltarho);
				d2mudTdrho_Trho[i][j] = (muplusT_plusrho - muplusT_minusrho - muminusT_plusrho + muminusT_minusrho)/(2*deltaT*deltarho);

				// The thermal conductivity values
				double kval =             IPropsSI(iL,iT,Tval,       iD,rhoval,         iFluid);
				double kplusrho =         IPropsSI(iL,iT,Tval,       iD,rhoval+deltarho,iFluid);
				double kminusrho =        IPropsSI(iL,iT,Tval,       iD,rhoval-deltarho,iFluid);
				double kplusT =           IPropsSI(iL,iT,Tval+deltaT,iD,rhoval,         iFluid);
				double kminusT =          IPropsSI(iL,iT,Tval-deltaT,iD,rhoval,         iFluid);
				double kplusT_plusrho =   IPropsSI(iL,iT,Tval+deltaT,iD,rhoval+deltarho,iFluid);
				double kplusT_minusrho =  IPropsSI(iL,iT,Tval+deltaT,iD,rhoval-deltarho,iFluid);
				double kminusT_plusrho =  IPropsSI(iL,iT,Tval-deltaT,iD,rhoval+deltarho,iFluid);
				double kminusT_minusrho = IPropsSI(iL,iT,Tval-deltaT,iD,rhoval-deltarho,iFluid);

				k_Trho[i][j] = kval;
				dkdT_Trho[i][j] = (-kminusT + kplusT)/(2*deltaT);
				dkdrho_Trho[i][j] = (-kminusrho + kplusrho)/(2*deltarho);
				d2kdT2_Trho[i][j] = (kminusT - 2*kval + kplusT)/(deltaT*deltaT);
				d2kdrho2_Trho[i][j] = (kminusrho - 2*kval + kplusrho)/(deltarho*deltarho);
				d2kdTdrho_Trho[i][j] = (kplusT_plusrho - kplusT_minusrho - kminusT_plusrho + kminusT_minusrho)/(2*deltaT*deltarho);
			}
			else
			{
				s_Trho[i][j] = _HUGE;
				dsdT_Trho[i][j] = _HUGE;
				dsdrho_Trho[i][j] = _HUGE;
				d2sdT2_Trho[i][j] = _HUGE;
				d2sdTdrho_Trho[i][j] = _HUGE;
				d2sdrho2_Trho[i][j] = _HUGE;

				h_Trho[i][j] = _HUGE;
				dhdT_Trho[i][j] = _HUGE;
				dhdrho_Trho[i][j] = _HUGE;
				d2hdT2_Trho[i][j] = _HUGE;
				d2hdTdrho_Trho[i][j] = _HUGE;
				d2hdrho2_Trho[i][j] = _HUGE;

				p_Trho[i][j] = _HUGE;
				dpdT_Trho[i][j] = _HUGE;
				dpdrho_Trho[i][j] = _HUGE;
				d2pdT2_Trho[i][j] = _HUGE;
				d2pdTdrho_Trho[i][j] = _HUGE;
				d2pdrho2_Trho[i][j] = _HUGE;

				mu_Trho[i][j] = _HUGE;
				dmudT_Trho[i][j] = _HUGE;
				dmudrho_Trho[i][j] = _HUGE;
				d2mudT2_Trho[i][j] = _HUGE;
				d2mudTdrho_Trho[i][j] = _HUGE;
				d2mudrho2_Trho[i][j] = _HUGE;

				k_Trho[i][j] = _HUGE;
				dkdT_Trho[i][j] = _HUGE;
				dkdrho_Trho[i][j] = _HUGE;
				d2kdT2_Trho[i][j] = _HUGE;
				d2kdTdrho_Trho[i][j] = _HUGE;
				d2kdrho2_Trho[i][j] = _HUGE;
			}
			//std::cout << format("%d %d\n",i,j);
		}
	}
	t2 = clock();
	double elap = (double)(t2-t1)/CLOCKS_PER_SEC;
	std::cout << elap << " to build single phase table for T,rho" << std::endl;

	// Update the boundaries of the points within the single-phase regions
	update_saturation_boundary_indices();

	update_cell_validity();

	return elap;
}

void TTSESinglePhaseTableClass::update_cell_validity()
{
	for (unsigned int i = 0; i < NT; i++)
	{
		for (unsigned int j = 0; j < Nrho; j++)
		{
			if (i < NT-1 && j < Nrho-1)
			{ 
				if(ValidNumber(s_Trho[i][j]) && ValidNumber(s_Trho[i][j+1]) && ValidNumber(s_Trho[i+1][j]) && ValidNumber(s_Trho[i+1][j+1]))
				{
					bicubic_cells.cells[i][j].valid_Trho = true;
				}
				else
				{
					bicubic_cells.cells[i][j].valid_Trho = false;
				}
			}
		}
	}

	for (unsigned int i = 0; i < Nh; i++)
	{
		for (unsigned int j = 0; j < Np; j++)
		{
			if (i < Nh-1 && j < Np-1)
			{ 
				if(ValidNumber(s[i][j]) && ValidNumber(s[i][j+1]) && ValidNumber(s[i+1][j]) && ValidNumber(s[i+1][j+1]))
				{
					bicubic_cells.cells[i][j].valid_hp = true;
				}
				else
				{
					bicubic_cells.cells[i][j].valid_hp = false;
				}
			}
		}
	}
}
//void TTSESinglePhaseTableClass::update_Trho_map()
//{
//	int ii,jj;
//	double Tmin,Tmax,rhomin,rhomax,rhoL,rhoV,TsatL,TsatV,dummy;
//	// Get the bounding values
//	Tmin = T[0][0];
//	rhomax = rho[0][0];
//	Tmax = T[Nh-1][Np-1];
//	rhomin = rho[Nh-1][0];
//	// Resize the arrays, using the same sizes as the base matrices
//	T_Trho.resize(Nh);
//	rho_Trho.resize(Np);
//	i_Trho.resize(Nh, std::vector<int>(Np, -1));
//	j_Trho.resize(Nh, std::vector<int>(Np, -1));
//	
//	for (unsigned int i = 0; i < Nh; i++)
//	{
//		double T = (Tmax-Tmin)/(Nh-1)*i+Tmin;
//		T_Trho[i] = T;
//
//		for (unsigned int j = 0; j < Np; j++)
//		{
//			double rho = (rhomax-rhomin)/(Np-1)*j+rhomin;
//			rho_Trho[j] = rho;
//
//			// T,rho --> p,h
//			double p = pFluid->pressure_Trho(T,rho);
//			double h = pFluid->enthalpy_Trho(T,rho);
//			double rhooV = pFluid->rhosatV(T);
//			double rhooL = pFluid->rhosatL(T);
//			double pV = pFluid->psatV_anc(T);
//			double pp = pFluid->pressure_Trho(T,rhooV);
//
//			// Find i,j from p,h
//			ii = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
//			jj = (int)round((log(p)-logpmin)/logpratio);
//
//			// Only keep values that are within the range for the table
//			if ( ii>=0 && ii < (int)Nh && jj>=0 && jj< (int)Np)
//			{
//				i_Trho[i][j] = ii;
//				j_Trho[i][j] = jj;
//			}
//			else
//			{
//				i_Trho[i][j] = -1;
//				j_Trho[i][j] = -1;
//			}
//		}
//	}
//}
void TTSESinglePhaseTableClass::update_saturation_boundary_indices()
{
	// Store some information about the phase boundaries so that we 
	// can use other inputs than p,h more easily

	for (unsigned int j = 0; j < Np; j++)
	{
		if (p[j] < pFluid->reduce.p.Pa)
		{
			IL[j] = -1;
			// Sweep left to right to find a phase boundary, use the first one that fails in the saturation region
			for (unsigned int i = 0; i < Nh; i++)
			{
				if (!ValidNumber(T[i][j]))
				{
					IL[j] = i;
					break;
				}
			}
			IV[j] = -1;
			// Sweep right to left to find a phase boundary, use the first one that fails in the saturation region
			for (int i = Nh-1; i > 0; i--)
			{
				if (!ValidNumber(T[i][j]))
				{
					IV[j] = i;
					break;	
				}
			}
			if (SatL != NULL && SatV != NULL)
			{
				TL[j] = SatL->evaluate(iT,p[j]);
				SL[j] = SatL->evaluate(iS,p[j]);
				DL[j] = SatL->evaluate(iD,p[j]);
				TV[j] = SatV->evaluate(iT,p[j]);
				SV[j] = SatV->evaluate(iS,p[j]);
				DV[j] = SatV->evaluate(iD,p[j]);
			}
			else
			{
				throw ValueError("SatL and SatV must be provided");
			}
		}
		else
		{
			IL[j] = -1;
			IV[j] = -1;
			TL[j] = _HUGE;
			SL[j] = _HUGE;
			DL[j] = _HUGE;
			TV[j] = _HUGE;
			SV[j] = _HUGE;
			DV[j] = _HUGE;
		}
	}
}
void TTSESinglePhaseTableClass::write_dotdrawing_tofile(char fName[])
{
	FILE *fp;
	fp = fopen(fName,"w");
	for (int j = Np-1; j>=0; j--)
	{
		for (unsigned int i = 0; i<Nh; i++)
		{
			if (ValidNumber(rho[i][j]))
			{
				fprintf(fp,".");
			}
			else
			{
				fprintf(fp,"X");
			}
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

double TTSESinglePhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *h, std::vector<double> *p, std::vector<double> *EOSv, std::vector<double> *TTSE)
{	
	double val=0;
	h->resize(N);
	p->resize(N);
	EOSv->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClassSI CPS(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;
		
		CPS.update(iH,h1,iP,p1);
		double sEOS = CPS.s();
		double cpEOS = CPS.cp();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*h)[i] = h1;
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,p1,CPS._logp,h1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOSv)[i] = sEOS; break;
		case iT:
			(*EOSv)[i] = TEOS; break;
		case iC:
			(*EOSv)[i] = cpEOS; break;
		case iD:
			(*EOSv)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g %g (h,p,EOS,TTSE,diff)\n",h1,p1,(*EOSv)[i],(*TTSE)[i],(*EOSv)[i]-(*TTSE)[i]).c_str();
	}
	return val;
}

double TTSESinglePhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		double h1 = ((double)rand()/(double)RAND_MAX)*(hmax-hmin)+hmin;

		if (p1 > pFluid->TTSESatL.pmax || h1 > pFluid->TTSESatV.evaluate(iH,p1) || h1 < pFluid->TTSESatL.evaluate(iH,p1))
		{
			// Get the value from TTSE
			evaluate(iParam,p1,log(p1),h1); 
		}
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}

double TTSESinglePhaseTableClass::evaluate(long iParam, double p, double logp, double h)
{
	// Use Bicubic interpolation if requested
	if (mode == TTSE_MODE_BICUBIC){ 
		return interpolate_bicubic_ph(iParam,p,logp,h);
	}
	int i = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
	int j = (int)round((logp-logpmin)/logpratio);
	
	if (i<0 || i>(int)Nh-1 || j<0 || j>(int)Np-1)
	{
		throw ValueError(format("Input to TTSE [p = %0.16g, h = %0.16g] is out of range",p,h));
	}

	// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
	// might be in the two-phase region which is not defined for single-phase table.  
	// Therefore, search around its neighbors for a better choice
	if (!ValidNumber(T[i][j])){
		nearest_good_neighbor(&i,&j);
	}

	// Distances from the node
	double deltap = p-this->p[j];
	double deltah = h-this->h[i];
	
	switch (iParam)
	{
	case iS:
		return s[i][j]+deltah*dsdh[i][j]+deltap*dsdp[i][j]+0.5*deltah*deltah*d2sdh2[i][j]+0.5*deltap*deltap*d2sdp2[i][j]+deltap*deltah*d2sdhdp[i][j]; break;
	case iT:
		return T[i][j]+deltah*dTdh[i][j]+deltap*dTdp[i][j]+0.5*deltah*deltah*d2Tdh2[i][j]+0.5*deltap*deltap*d2Tdp2[i][j]+deltap*deltah*d2Tdhdp[i][j]; break;
	case iD:
		return rho[i][j]+deltah*drhodh[i][j]+deltap*drhodp[i][j]+0.5*deltah*deltah*d2rhodh2[i][j]+0.5*deltap*deltap*d2rhodp2[i][j]+deltap*deltah*d2rhodhdp[i][j]; break;
	default:
		throw ValueError(format("Output key value [%d] to evaluate is invalid",iParam));
	}
}

void TTSESinglePhaseTableClass::bicubic_cell_coordinates_Trho(double Tval, double rhoval, double logrhoval, int *i, int *j)
{
	*i = (int)round((Tval-Tmin)/(Tmax-Tmin)*(NT-1));
	*j = (int)round((logrhoval-logrhomin)/logrhoratio);

	if (*i<0 || *i>(int)NT-1 || *j<0 || *j>(int)Nrho-1)
	{
		throw ValueError(format("Input to TTSE [T = %0.16g, logrho = %0.16g] is out of range",Tval,logrhoval));
	}
	
	if (Tval < this->T_Trho[*i])
	{
		*i -= 1;
	}
	if (rhoval < this->rho_Trho[*j])
	{
		*j -= 1;
	}

	nearest_good_neighbor_Trho_interpolate(i, j);
}

void TTSESinglePhaseTableClass::bicubic_cell_coordinates_ph(double hval, double pval, double logpval, int *i, int *j)
{
	*i = (int)round(((hval-hmin)/(hmax-hmin)*(Nh-1)));
	*j = (int)round((logpval-logpmin)/logpratio);

	if (*i<0 || *i>(int)Nh-1 || *j<0 || *j>(int)Np-1)
	{
		throw ValueError(format("Input to TTSE [p = %0.16g, h = %0.16g] is out of range",pval,hval));
	}
	
	if (hval < this->h[*i])
	{
		*i -= 1;
	}
	if (pval < p[*j])
	{
		*j -= 1;
	}

	nearest_good_neighbor_ph_interpolate(i, j);
}

std::vector<double> * TTSESinglePhaseTableClass::bicubic_cell_coeffs_ph(long iParam, int i, int j)
{
	std::vector<double> *alpha = NULL;
	std::vector< std::vector<double> > *f = NULL, *dfdp = NULL, *dfdh = NULL, *d2fdhdp = NULL;
	
	switch(iParam)
	{
	case iS:
		if (!bicubic_cells.cells[i][j].alpha_s_hp.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_s_hp;
			break;
		}
		else
		{
			f = &s; dfdh = &dsdh; dfdp = &dsdp; d2fdhdp = &d2sdhdp; break;
		}
	case iT:
		if (!bicubic_cells.cells[i][j].alpha_T_hp.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_T_hp;
			break;
		}
		else
		{
			f = &T; dfdh = &dTdh; dfdp = &dTdp; d2fdhdp = &d2Tdhdp; break;
		}
	case iD:
		if (!bicubic_cells.cells[i][j].alpha_rho_hp.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_rho_hp;
			break;
		}
		else
		{
			f = &rho; dfdh = &drhodh; dfdp = &drhodp; d2fdhdp = &d2rhodhdp; break;
		}
	default:
		throw ValueError("iParam is invalid");
	}

	if (alpha == NULL)
	{
		double dhdx = (h[i+1]-h[i]);
		double dpdy = (p[j+1]-p[j]);

		z_bicubic[0] = (*f)[i][j];
		z_bicubic[1] = (*f)[i+1][j];
		z_bicubic[2] = (*f)[i][j+1];
		z_bicubic[3] = (*f)[i+1][j+1];
		
		z_bicubic[4] = (*dfdh)[i][j]*dhdx;
		z_bicubic[5] = (*dfdh)[i+1][j]*dhdx;
		z_bicubic[6] = (*dfdh)[i][j+1]*dhdx;
		z_bicubic[7] = (*dfdh)[i+1][j+1]*dhdx;
		
		z_bicubic[8] = (*dfdp)[i][j]*dpdy;
		z_bicubic[9] = (*dfdp)[i+1][j]*dpdy;
		z_bicubic[10] = (*dfdp)[i][j+1]*dpdy;
		z_bicubic[11] = (*dfdp)[i+1][j+1]*dpdy;

		z_bicubic[12] = (*d2fdhdp)[i][j]*dpdy*dhdx;
		z_bicubic[13] = (*d2fdhdp)[i+1][j]*dpdy*dhdx;
		z_bicubic[14] = (*d2fdhdp)[i][j+1]*dpdy*dhdx;
		z_bicubic[15] = (*d2fdhdp)[i+1][j+1]*dpdy*dhdx;

		// Find the alpha values by doing the matrix operation alpha = Ainv*z
		matrix_vector_product(&z_bicubic, &alpha_bicubic);
		// Set the pointer in this function
		alpha = &alpha_bicubic;
		
		// Store this array of bicubic coefficients
		switch (iParam)
		{
		case iS:	
			bicubic_cells.cells[i][j].alpha_s_hp = alpha_bicubic; break;
		case iT:
			bicubic_cells.cells[i][j].alpha_T_hp = alpha_bicubic; break;
		case iD:
			bicubic_cells.cells[i][j].alpha_rho_hp = alpha_bicubic; break;
		default:
			throw ValueError("iParam is invalid");
		}
	}
	return alpha;
}

std::vector<double> * TTSESinglePhaseTableClass::bicubic_cell_coeffs_Trho(long iParam, int i, int j)
{

	std::vector<double> *alpha = NULL;
	std::vector< std::vector<double> > *f = NULL, *dfdT = NULL, *dfdrho = NULL, *d2fdTdrho = NULL;

	// Find the coefficients for the bicubic function
	switch(iParam)
	{
	case iS:
		if (!bicubic_cells.cells[i][j].alpha_s_Trho.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_s_Trho;
			break;
		}
		else
		{
			f = &s_Trho; dfdT = &dsdT_Trho; dfdrho = &dsdrho_Trho; d2fdTdrho = &d2sdTdrho_Trho; break;
		}
	case iH:
		if (!bicubic_cells.cells[i][j].alpha_h_Trho.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_h_Trho;
			break;
		}
		else
		{
			f = &h_Trho; dfdT = &dhdT_Trho; dfdrho = &dhdrho_Trho; d2fdTdrho = &d2hdTdrho_Trho; break;
		}
	case iP:
		if (!bicubic_cells.cells[i][j].alpha_p_Trho.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_p_Trho;
			break;
		}
		else
		{
			f = &p_Trho; dfdT = &dpdT_Trho; dfdrho = &dpdrho_Trho; d2fdTdrho = &d2pdTdrho_Trho; break;
		}
	case iV:
		if (!bicubic_cells.cells[i][j].alpha_mu_Trho.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_mu_Trho;
			break;
		}
		else
		{
			f = &mu_Trho; dfdT = &dmudT_Trho; dfdrho = &dmudrho_Trho; d2fdTdrho = &d2mudTdrho_Trho; break;
		}
	case iL:
		if (!bicubic_cells.cells[i][j].alpha_k_Trho.empty())
		{
			// Make alpha point to the pre-calculated values
			alpha = &bicubic_cells.cells[i][j].alpha_k_Trho;
			break;
		}
		else
		{
			f = &k_Trho; dfdT = &dkdT_Trho; dfdrho = &dkdrho_Trho; d2fdTdrho = &d2kdTdrho_Trho; break;
		}
	default:
		throw ValueError(format("iParam [%d] is invalid",iParam).c_str());
	}

	if (alpha == NULL)
	{
		double dTdx = (T_Trho[i+1]-T_Trho[i]);
		double drhody = (rho_Trho[j+1]-rho_Trho[j]);

		z_bicubic[0] = (*f)[i][j];
		z_bicubic[1] = (*f)[i+1][j];
		z_bicubic[2] = (*f)[i][j+1];
		z_bicubic[3] = (*f)[i+1][j+1];
		
		z_bicubic[4] = (*dfdT)[i][j]*dTdx;
		z_bicubic[5] = (*dfdT)[i+1][j]*dTdx;
		z_bicubic[6] = (*dfdT)[i][j+1]*dTdx;
		z_bicubic[7] = (*dfdT)[i+1][j+1]*dTdx;
		
		z_bicubic[8] = (*dfdrho)[i][j]*drhody;
		z_bicubic[9] = (*dfdrho)[i+1][j]*drhody;
		z_bicubic[10] = (*dfdrho)[i][j+1]*drhody;
		z_bicubic[11] = (*dfdrho)[i+1][j+1]*drhody;

		z_bicubic[12] = (*d2fdTdrho)[i][j]*drhody*dTdx;
		z_bicubic[13] = (*d2fdTdrho)[i+1][j]*drhody*dTdx;
		z_bicubic[14] = (*d2fdTdrho)[i][j+1]*drhody*dTdx;
		z_bicubic[15] = (*d2fdTdrho)[i+1][j+1]*drhody*dTdx;

		// Find the alpha values by doing the matrix operation alpha = Ainv*z
		matrix_vector_product(&z_bicubic, &alpha_bicubic);
		// Set the pointer in this function
		alpha = &alpha_bicubic;
		
		// Store this array of bicubic coefficients
		switch (iParam)
		{
		case iS:	
			bicubic_cells.cells[i][j].alpha_s_Trho = alpha_bicubic; break;
		case iH:
			bicubic_cells.cells[i][j].alpha_h_Trho = alpha_bicubic; break;
		case iP:
			bicubic_cells.cells[i][j].alpha_p_Trho = alpha_bicubic; break;
		case iV:
			bicubic_cells.cells[i][j].alpha_mu_Trho = alpha_bicubic; break;
		case iL:
			bicubic_cells.cells[i][j].alpha_k_Trho = alpha_bicubic; break;
		default:
			throw ValueError(format("iParam [%d] is invalid",iParam).c_str());
		}
	}

	return alpha;
}

/*!

\f[
A^{-1} =  \left[ \begin{array}{*{16}c} 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & -3 & 3 & 0 & 0 & -2 & -1 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 2 & -2 & 0 & 0 & 1 & 1 & 0 & 0 \\
 -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & -3 & 0 & 3 & 0 & 0 & 0 & 0 & 0 & -2 & 0 & -1 & 0 \\
 9 & -9 & -9 & 9 & 6 & 3 & -6 & -3 & 6 & -6 & 3 & -3 & 4 & 2 & 2 & 1 \\
 -6 & 6 & 6 & -6 & -3 & -3 & 3 & 3 & -4 & 4 & -2 & 2 & -2 & -2 & -1 & -1 \\
 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 & 1 & 0 & 1 & 0 \\
 -6 & 6 & 6 & -6 & -4 & -2 & 4 & 2 & -3 & 3 & -3 & 3 & -2 & -1 & -2 & -1 \\
 4 & -4 & -4 & 4 & 2 & 2 & -2 & -2 & 2 & -2 & 2 & -2 & 1 & 1 & 1 & 1 \end{array} \right]
 \f]

 \f[
 x = \frac{h-h_i}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial x}{\partial h} = \frac{1}{h_{i+1}-h_i}
 \f]
 \f[
 \frac{\partial h}{\partial x} = h_{i+1}-h_i
 \f]

 \f[
 y = \frac{p-p_j}{p_{j+1}-p_j}
 \f]
 \f[
 \frac{\partial y}{\partial p} = \frac{1}{p_{j+1}-p_j}
 \f]
  \f[
 \frac{\partial p}{\partial y} = p_{j+1}-p_j
 \f]

 \f[
 \frac{\partial f}{\partial x} = \frac{\partial f}{\partial h}\cdot\frac{\partial h}{\partial x}
 \f]
 \
*/
double TTSESinglePhaseTableClass::interpolate_bicubic_ph(long iParam, double pval, double logp, double hval)
{
	std::vector<double> *alpha = NULL;
	int i, j;

	// Get the i,j coordinates of the cell
	this->bicubic_cell_coordinates_ph(hval, pval,logp, &i, &j);

	double dhdx = (h[i+1]-h[i]);
	double dpdy = (p[j+1]-p[j]);
	double x = (hval-h[i])/dhdx;
	double y = (pval-p[j])/dpdy;

	// Find the coefficients for the cubic
	alpha = this->bicubic_cell_coeffs_ph(iParam, i, j);

	double x_0 = 1, x_1 = x, x_2 = x*x, x_3 = x*x*x;
	double y_0 = 1, y_1 = y, y_2 = y*y, y_3 = y*y*y;

	double summer;
	// Old version was "summer += (*alpha)[ii+jj*4]*x^i*y^j" for i in [0,3] and j in [0,3]
	summer = (*alpha)[0+0*4]*x_0*y_0+(*alpha)[0+1*4]*x_0*y_1+(*alpha)[0+2*4]*x_0*y_2+(*alpha)[0+3*4]*x_0*y_3
		    +(*alpha)[1+0*4]*x_1*y_0+(*alpha)[1+1*4]*x_1*y_1+(*alpha)[1+2*4]*x_1*y_2+(*alpha)[1+3*4]*x_1*y_3
			+(*alpha)[2+0*4]*x_2*y_0+(*alpha)[2+1*4]*x_2*y_1+(*alpha)[2+2*4]*x_2*y_2+(*alpha)[2+3*4]*x_2*y_3
			+(*alpha)[3+0*4]*x_3*y_0+(*alpha)[3+1*4]*x_3*y_1+(*alpha)[3+2*4]*x_3*y_2+(*alpha)[3+3*4]*x_3*y_3;
	
	return summer;
}
double TTSESinglePhaseTableClass::bicubic_evaluate_first_derivative_ph(long iOF, long iWRT, long iCONSTANT, double p, double logp, double h)
{
	std::vector<double> *alpha = NULL;
	int i, j;

	// Get the i,j coordinates of the cell
	this->bicubic_cell_coordinates_ph(h, p, logp, &i, &j);

	double dhdx = (this->h[i+1]-this->h[i]);
	double dpdy = (this->p[j+1]-this->p[j]);
	double x = (h-this->h[i])/dhdx;
	double y = (p-this->p[j])/dpdy;

	// Find the coefficients for the cubic
	alpha = this->bicubic_cell_coeffs_ph(iOF, i, j);

	double x_0 = 1, x_1 = x, x_2 = x*x, x_3 = x*x*x;
	double y_0 = 1, y_1 = y, y_2 = y*y, y_3 = y*y*y;

	double summer;
	if (iWRT == iH)
	{
		// Old version was "summer += (*alpha)[ii+jj*4]*i*x^(i-1)*y^j" for i in [0,3] and j in [0,3]
		double dsummer_dx = +(*alpha)[1+0*4]*1*x_0*y_0+(*alpha)[1+1*4]*1*x_0*y_1+(*alpha)[1+2*4]*1*x_0*y_2+(*alpha)[1+3*4]*1*x_0*y_3
				            +(*alpha)[2+0*4]*2*x_1*y_0+(*alpha)[2+1*4]*2*x_1*y_1+(*alpha)[2+2*4]*2*x_1*y_2+(*alpha)[2+3*4]*2*x_1*y_3
				            +(*alpha)[3+0*4]*3*x_2*y_0+(*alpha)[3+1*4]*3*x_2*y_1+(*alpha)[3+2*4]*3*x_2*y_2+(*alpha)[3+3*4]*3*x_2*y_3;
		return dsummer_dx/dhdx;
	}
	else if (iWRT == iP)
	{
		// Old version was "summer += (*alpha)[ii+jj*4]*j*x^i*y^(j-1)" for i in [0,3] and j in [0,3]
		double dsummer_dy = (*alpha)[0+1*4]*1*x_0*y_0+(*alpha)[0+2*4]*2*x_0*y_1+(*alpha)[0+3*4]*3*x_0*y_2
					       +(*alpha)[1+1*4]*1*x_1*y_0+(*alpha)[1+2*4]*2*x_1*y_1+(*alpha)[1+3*4]*3*x_1*y_2
					       +(*alpha)[2+1*4]*1*x_2*y_0+(*alpha)[2+2*4]*2*x_2*y_1+(*alpha)[2+3*4]*3*x_2*y_2
					       +(*alpha)[3+1*4]*1*x_3*y_0+(*alpha)[3+2*4]*2*x_3*y_1+(*alpha)[3+3*4]*3*x_3*y_2;
		return dsummer_dy/dpdy;
	}
	else
	{
		throw ValueError(format("iWRT [%d] is invalid", iWRT).c_str());
	}
	
	return summer;
}

double TTSESinglePhaseTableClass::interpolate_bicubic_Trho(long iParam, double Tval, double rhoval, double logrho)
{
	// A pointer to the values of the coefficients for the bicubic interpolation
	std::vector<double> *alpha = NULL;
	int i,j;
	std::vector< std::vector<double> > *f = NULL, *dfdT = NULL, *dfdrho = NULL, *d2fdTdrho = NULL;
	
	bicubic_cell_coordinates_Trho(Tval, rhoval, logrho, &i, &j);

	double dTdx = (T_Trho[i+1]-T_Trho[i]);
	double drhody = (rho_Trho[j+1]-rho_Trho[j]);
	double x = (Tval-T_Trho[i])/dTdx;
	double y = (rhoval-rho_Trho[j])/drhody;

	// Find the coefficients for the cubic
	alpha = this->bicubic_cell_coeffs_Trho(iParam, i, j);

	double x_0 = 1, x_1 = x, x_2 = x*x, x_3 = x*x*x;
	double y_0 = 1, y_1 = y, y_2 = y*y, y_3 = y*y*y;

	double summer;
	// Old version was "summer += (*alpha)[ii+jj*4]*x^i*y^j" for i in [0,3] and j in [0,3]
	summer = (*alpha)[0+0*4]*x_0*y_0+(*alpha)[0+1*4]*x_0*y_1+(*alpha)[0+2*4]*x_0*y_2+(*alpha)[0+3*4]*x_0*y_3
		    +(*alpha)[1+0*4]*x_1*y_0+(*alpha)[1+1*4]*x_1*y_1+(*alpha)[1+2*4]*x_1*y_2+(*alpha)[1+3*4]*x_1*y_3
			+(*alpha)[2+0*4]*x_2*y_0+(*alpha)[2+1*4]*x_2*y_1+(*alpha)[2+2*4]*x_2*y_2+(*alpha)[2+3*4]*x_2*y_3
			+(*alpha)[3+0*4]*x_3*y_0+(*alpha)[3+1*4]*x_3*y_1+(*alpha)[3+2*4]*x_3*y_2+(*alpha)[3+3*4]*x_3*y_3;
	
	return summer;
}

bool isbetween(double x1, double x2, double x)
{
	return ((x>=x1 && x <= x2) || (x>=x2 && x <= x1));
}

double TTSESinglePhaseTableClass::bicubic_evaluate_one_other_input(long iInput1, double Input1, long iOther, double Other)
{
	// My goal here is to find deltah
	int L,R,M,i,j,phase;
	double p,TL, TV, DL, DV, SL, SV, HL, HV, Q, ValV, ValL;
	std::vector<std::vector<double> > *mat;
	
	// Connect a pointer to the array of interest
	switch (iOther)
	{
	case iT:
		mat = &T; break;
	case iS:
		mat = &s; break;
	case iD:
		mat = &rho; break;
	}

	// One is pressure, we are getting enthalpy
	if (iInput1 == iP)
	{
		p = Input1;
		
		j = (int)round((log(p)-logpmin)/logpratio);

		if (p > pFluid->reduce.p.Pa)
		{
			// It's supercritical, we just need to iterate on the other parameter

			// Do interval halving over the whole range since
			// there can't be any saturation curve
			L = 0; R = Nh-1; M = (L+R)/2;
			while (R-L>1){
				if (isbetween((*mat)[M][j],(*mat)[R][j],Other)){
					L=M; M=(L+R)/2; continue;
				}
				else{
					R=M; M=(L+R)/2; continue;
				}
			}
			i = L;
		}
		else
		{
			// Get the saturation states
			switch (iOther)
			{
			case iT:
				TL = SatL->evaluate(iT,p); TV = SatV->evaluate(iT,p); ValL = TL; ValV = TV;
				if (Other > TV){ phase = iPHASE_GAS;} else if  (Other < TL){ phase = iPHASE_LIQUID;} else {phase = iPHASE_TWOPHASE; throw ValueError("Two-phase T,P inputs invalid for TTSE");}
				break;
			case iD:
				DL = SatL->evaluate(iD,p); DV = SatV->evaluate(iD,p); ValL = DL; ValV = DV;
				if (Other < DV){ phase = iPHASE_GAS;} else if  (Other > DL){ phase = iPHASE_LIQUID;} else {phase = iPHASE_TWOPHASE; Q = (1/Other-1/DL)/(1/DV-1/DL);}
				break;
			case iS:
				SL = SatL->evaluate(iS,p); SV = SatV->evaluate(iS,p); ValL = SL; ValV = SV;
				if (Other > SV){ phase = iPHASE_GAS;} else if  (Other < SL){ phase = iPHASE_LIQUID;} else {phase = iPHASE_TWOPHASE; Q = (Other-SL)/(SV-SL);}
				break;
			}

			if (phase == iPHASE_GAS)
			{
				L = IV[j]+1; R = Nh - 1; M = (L+R)/2;
				// Check if caught between saturation curve and first point
				if (isbetween((*mat)[L][j],ValV,Other))
				{
					// If caught, use gas value
					L = R;
				}
			}
			else if (phase == iPHASE_LIQUID)
			{
				L = 0; R = IL[j]-1; M = (L+R)/2;
				// Check if caught between saturation curve and first point
				if (isbetween((*mat)[R][j],ValL,Other))
				{
					// If caught, use liquid value
					R = L;
				}
			}
			else if (phase == iPHASE_TWOPHASE)
			{
				// Two-phase, finished
				HL = SatL->evaluate(iH,p); HV = SatV->evaluate(iH,p);
				return HV*Q + (1-Q)*HL;
			}
			else
			{
				throw ValueError("TTSE BC error on phase");
			}

			// Interval bisection
			while (R-L>1)
			{
				if (isbetween((*mat)[M][j],(*mat)[R][j],Other))
				{ 
					L=M; M=(L+R)/2; continue;
				}
				else
				{ 
					R=M; M=(L+R)/2; continue;
				}
			}
			i = L;
		}
		// Get the nearest full cell to this point
		nearest_good_neighbor_ph_interpolate(&i, &j);
		
		// Invert the bicubic to find h
		// Old version was "summer += (*alpha)[ii+jj*4]*x^i*y^j" for i in [0,3] and j in [0,3]
		/* 
		Summation function is 
		summer = (*alpha)[0+0*4]*x_0*y_0+(*alpha)[0+1*4]*x_0*y_1+(*alpha)[0+2*4]*x_0*y_2+(*alpha)[0+3*4]*x_0*y_3
		  		+(*alpha)[1+0*4]*x_1*y_0+(*alpha)[1+1*4]*x_1*y_1+(*alpha)[1+2*4]*x_1*y_2+(*alpha)[1+3*4]*x_1*y_3
				+(*alpha)[2+0*4]*x_2*y_0+(*alpha)[2+1*4]*x_2*y_1+(*alpha)[2+2*4]*x_2*y_2+(*alpha)[2+3*4]*x_2*y_3
				+(*alpha)[3+0*4]*x_3*y_0+(*alpha)[3+1*4]*x_3*y_1+(*alpha)[3+2*4]*x_3*y_2+(*alpha)[3+3*4]*x_3*y_3;

		but in this case we know y, which is simply

		y = (p-this->p[j])/(this->p[j+1]-this->p[j]);

		and we solve for x, which is cubic.  Convert to a form like
		
		0 = ax^3 + b*x^2 + c*x + d

		for use in solve_cubic where 

		x = (h-this->h[i])/(this->h[i+1]-this->h[i]);
		
		*/
		double y = (p-this->p[j])/(this->p[j+1]-this->p[j]);
		double y_0 = 1, y_1 = y, y_2 = y*y, y_3 = y*y*y;

		// Find the coefficients for the cubic
		std::vector<double> *alpha = NULL;
		
		// Get bicubic coefficients
		alpha = this->bicubic_cell_coeffs_ph(iOther, i, j);

		double a = (*alpha)[3+0*4]*y_0+(*alpha)[3+1*4]*y_1+(*alpha)[3+2*4]*y_2+(*alpha)[3+3*4]*y_3; // factors of x^3
		double b = (*alpha)[2+0*4]*y_0+(*alpha)[2+1*4]*y_1+(*alpha)[2+2*4]*y_2+(*alpha)[2+3*4]*y_3; // factors of x^2
		double c = (*alpha)[1+0*4]*y_0+(*alpha)[1+1*4]*y_1+(*alpha)[1+2*4]*y_2+(*alpha)[1+3*4]*y_3; // factors of x
		double d = (*alpha)[0+0*4]*y_0+(*alpha)[0+1*4]*y_1+(*alpha)[0+2*4]*y_2+(*alpha)[0+3*4]*y_3 - Other; // constant factors
		double x0,x1,x2;
		solve_cubic(a,b,c,d,&x0,&x1,&x2);
		
		double x;
		// Only one solution
		if (fabs(x0-x1) < 10*DBL_EPSILON &&  fabs(x0-x2) < 10*DBL_EPSILON && fabs(x1-x2) < 10*DBL_EPSILON){
			x = x0;
		}
		// See if first is good solution 0 < x < 1
		else if ( 0 <= x0 && x0 <= 1){
			x = x0;
		}
		else if ( 0 <= x1 && x1 <= 1){
			x = x1;
		}
		else if ( 0 <= x2 && x2 <= 1){
			x = x2;
		}
		else
		{
			throw ValueError("Multiple solutions found for cubic in TTSE-BICUBIC-P+Other");
		}
		return x*(this->h[i+1]-this->h[i])+this->h[i];
	}
	// One is enthalpy, we are getting pressure
	else if (iInput1 == iH)
	{
		throw ValueError("Sorry enthalpy and something else other than p is not valid for TTSE");
	}
	// Oops, neither enthalpy or pressure provided
	else
	{
		throw ValueError("Neither enthalpy nor pressure provided as the first parameter");
	}
}

bool TTSESinglePhaseTableClass::within_range_Trho(long iInput1, double Input1, long iOther, double Other)
{
	int i = (int)round((Input1-Tmin)/(Tmax-Tmin)*(NT-1));
	int j = (int)round((log(Other)-logrhomin)/logrhoratio);
	return (0 <= i && i <= (int)NT && 0 <= (int)j && j <= (int)Nrho);
}
bool TTSESinglePhaseTableClass::within_range_one_other_input(long iInput1, double Input1, long iOther, double Other)
{
	if (iInput1 != iP)
	{
		throw ValueError("iInput1 must be iP");
	}
	// Check if pressure is in range
	if (Input1 > pmax || Input1 < pmin){ 
		return false;
	}

	int j = (int)round((log(Input1)-logpmin)/logpratio);

	if (iOther == iT || iOther == iS)
	{
		double right, left;
		if (iOther == iT) 
			{ right = this->T[Nh-1][j]; left = this->T[0][j]; }
		else 
			{ right = this->s[Nh-1][j]; left = this->s[0][j]; }

		if (ValidNumber(right) && ValidNumber(left))
		{
			if (Other > right || Other < left)
				{ return false;	}
			else
				{ return true; }
		}
		else
		{
			return false;
		}
	}
	else if (iOther == iD)
	{
		double right = this->rho[Nh-1][j], left = this->rho[0][j];

		if (ValidNumber(right) && ValidNumber(left))
		{
			if (Other < right || Other > left)
				{ return false;	}
			else
				{ return true; }
		}
		else
		{
			return false;
		}
	}

	return true;
}
double TTSESinglePhaseTableClass::evaluate_one_other_input(long iInput1, double Input1, long iOther, double Other)
{
	// Use Bicubic interpolation if requested
	if (mode == TTSE_MODE_BICUBIC){ 
		return bicubic_evaluate_one_other_input(iInput1,Input1,iOther,Other);
	}
	// My goal here is to find deltah
	int L,R,M,i;
	double p,dh1,dh2,a,b,c;
	std::vector<std::vector<double> > *mat;
	
	// Connect a pointer to the array of interest
	switch (iOther)
	{
	case iT:
		mat = &T; break;
	case iS:
		mat = &s; break;
	case iD:
		mat = &rho; break;
	}

	// One is pressure, we are getting enthalpy
	// within_TTSE_range() guarantees that the pressure is within range
	if (iInput1 == iP)
	{
		p = Input1;
		
		int j = (int)round((log(p)-logpmin)/logpratio);
		double deltap = p-this->p[j];

		// Check the bounding values for the other input to see if it is within range

		if (j >= jpcrit_ceil) // Is is either supercritical pressure or just a little bit below critical pressure
		{
			// Very close to the boundary of the LUT, not 1-1 relationship between p-h and other
			// sets of inputs, need to allow for a bit of raggedness here
			if (   (iOther == iT && Other < 1.1*this->T[Np-1][j] && Other > this->T[Np-1][j])
			    || (iOther == iS && Other < this->s[Np-1][j]+0.1 && Other > this->s[Np-1][j])
			    || (iOther == iD && Other > 0.85*this->rho[Np-1][j] && Other < this->rho[Np-1][j])
			    )
			{
				i = Np-1;
			}
			else
			{
				// Do interval halving over the whole range since 
				// there can't be any saturation curve
				L = 0; R = Nh-1; M = (L+R)/2;
				while (R-L>1){
					if (isbetween((*mat)[M][j],(*mat)[R][j],Other)){ 
						L=M; M=(L+R)/2; continue;
					}
					else{ 
						R=M; M=(L+R)/2; continue;
					}
				}
				// Find which one of the bounds is closer
				if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other)){
					i = L;
				}
				else{
					i = R;
				}
			}
		}
		else
		{
			if (   (iOther == iT && Other > SatV->evaluate(iT,p))
				|| (iOther == iS && Other > SatV->evaluate(iS,p))
				|| (iOther == iD && Other < SatV->evaluate(iD,p))
				)
			{
				// SUPERHEATED!!
				//
				// If it is within between the saturation curve and the first point in the SH region,
				// just use the first point in the superheated region
				if (   (iOther == iT && Other < this->T[IV[j]+1][j])
				    || (iOther == iS && Other < this->s[IV[j]+1][j])
					|| (iOther == iD && Other > this->rho[IV[j]+1][j])
					)
				{
					i = IV[j]+1;
				}
				// Very close to the boundary of the LUT, not 1-1 relationship between p-h and other
				// sets of inputs, need to allow for a bit of raggedness here
				else if (   (iOther == iT && Other < 1.1*this->T[Np-1][j] && Other > this->T[Np-1][j])
				         || (iOther == iS && Other < this->s[Np-1][j]+0.1 && Other > this->s[Np-1][j])
					     || (iOther == iD && Other > 0.85*this->rho[Np-1][j] && Other < this->rho[Np-1][j])
					     )
				{
					i = Np-1;
				}
				else
				{
					// SUPERHEATED!! (and away from saturation)
					// Make sure it is in the bounds
					switch (iOther)
					{
					case iT:
						if (Other > this->T[Nh-1][j]) 
							throw ValueError(format("Input T [%g] is greater than max T [%g] at LUT pressure [%g]",Other,this->T[Np-1][j],this->p[j]));
						break;
					case iD:
						if (Other < this->rho[Nh-1][j])
							throw ValueError(format("Input rho [%g] is less than minimum rho [%g] at LUT pressure [%g]",Other,this->rho[Np-1][j],this->p[j]));
						break;
					case iS:
						if (Other > this->s[Nh-1][j]) 
							throw ValueError(format("Input s [%g] is greater than max s [%g] at LUT pressure [%g]",Other,this->s[Np-1][j],this->p[j]));
						break;
					}
					
					L = IV[j]+1; R = Np-1; M = (L+R)/2;
					while (R-L>1)
					{
						if (isbetween((*mat)[M][j],(*mat)[R][j],Other))
						{ 
							L=M; M=(L+R)/2; continue;
						}
						else
						{ 
							R=M; M=(L+R)/2; continue;
						}
					}
					// Find which one of the bounds is closer
					if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other))
						i = L;
					else
						i = R;
				}
			}
			else if (IL[j] < 2)
			{
				// We are at low pressure, so we don't know how to calculate, going to just use the i==1 element
				// if it is valid, or the i = 0 if not, otherwise, there are no values left and we have to fail
				if (ValidNumber(this->T[1][j])){
					i = 1;
				}
				else if (ValidNumber(this->T[0][j])){
					i = 0;
				}
				else{
					throw ValueError(format("Your inputs [%g,%g] do not yield a valid TTSE node",Input1,Other));
				}
			}
			
			else if (   (iOther == iT && Other < SatL->evaluate(iT,p))
				     || (iOther == iS && Other < SatL->evaluate(iS,p))
					 || (iOther == iD && Other > SatL->evaluate(iD,p))
					 )
			{
				// SUBCOOLED!!
				//
				// If it is within one spacing of the outlet variable of the saturation curve, 
				// just use the first point in the subcooled region
				if (   (iOther == iT && Other > this->T[IL[j]-1][j])
				    || (iOther == iS && Other > this->s[IL[j]-1][j])
					|| (iOther == iD && Other < this->rho[IL[j]-1][j])
					)
				{
					i = IL[j]-1;
				}
				else{
					// Make sure it is in the bounds of the LUT
					switch (iOther)
					{
					case iT:
						if (Other < this->T[0][j]) 
							throw ValueError(format("Input T [%g] is less than min T [%g] at LUT pressure [%g]",Other,this->T[0][j],this->p[j]));
						break;
					case iD:
						if (Other > this->rho[0][j]) 
							throw ValueError(format("Input rho [%g] is greater than max rho [%g] at LUT pressure [%g]",Other,this->rho[0][j],this->p[j]));
						break;
					case iS:
						if (Other < this->s[0][j])
							throw ValueError(format("Input s [%g] is less than min s [%g] at LUT pressure [%g]",Other,this->s[0][j],this->p[j]));
						break;
					}
					
					L = 0; R = IL[j]-1; M = (L+R)/2;
					// Its subcooled
					while (R-L>1)
					{
						if (isbetween((*mat)[M][j],(*mat)[R][j],Other))
						{ 
							L=M; M=(L+R)/2; continue;
						}
						else
						{ 
							R=M; M=(L+R)/2; continue;
						}
					}
					// Find which one of the bounds is closer
					if (fabs((*mat)[L][j]-Other)<fabs((*mat)[R][j]-Other))
						i = L;
					else
						i = R;
				}
			}
			else
			{
				// It's two-phase
				throw ValueError(format("It's two phase input"));
			}
		}

		// Now we calculate deltah
		switch (iOther)
		{
		// Quadratic in deltah
		// 0 = 0.5*deltah*deltah*d2Tdh2[i][j]+deltah*dTdh[i][j]+deltap*deltah*d2Tdhdp[i][j]+T[i][j]-T+deltap*dTdp[i][j]+0.5*deltap*deltap*d2Tdp2[i][j];
		// 0 = a*deltah^2+b*deltah+c
		case iT:
			a = 0.5*d2Tdh2[i][j];
			b = dTdh[i][j]+deltap*d2Tdhdp[i][j];
			c = T[i][j]-Other+deltap*dTdp[i][j]+0.5*deltap*deltap*d2Tdp2[i][j];
			break;
		case iS:
			a = 0.5*d2sdh2[i][j];
			b = dsdh[i][j]+deltap*d2sdhdp[i][j];
			c = s[i][j]-Other+deltap*dsdp[i][j]+0.5*deltap*deltap*d2sdp2[i][j];
			break;
		case iD:
			a = 0.5*d2rhodh2[i][j];
			b = drhodh[i][j]+deltap*d2rhodhdp[i][j];
			c = rho[i][j]-Other+deltap*drhodp[i][j]+0.5*deltap*deltap*d2rhodp2[i][j];
			break;
		}
		// Solutions from quadratic equation
		dh1 = (-b+sqrt(b*b-4*a*c))/(2*a);
		dh2 = (-b-sqrt(b*b-4*a*c))/(2*a);

		double hspacing = (hmax-hmin)/((double)Nh-1);
		// If only one is less than a multiple of enthalpy spacing, thats your solution
		if (fabs(dh1) < 10*hspacing && !(fabs(dh2) < 10*hspacing) )
			return dh1+this->h[i];
		else if (fabs(dh2) < 10*hspacing && !(fabs(dh1) < 10*hspacing) )
			return dh2+this->h[i];
		else{
			// Need to figure out which is the correct solution, try just the smaller one
			if (fabs(dh1)<fabs(dh2))
				return dh1+this->h[i];
			else
				return dh2+this->h[i];
		}
	}
	// One is enthalpy, we are getting pressure
	else if (iInput1 == iH)
	{
		throw ValueError("Sorry enthalpy and something else other than p is not valid for TTSE");
	}
	// Oops, neither enthalpy or pressure provided
	else
	{
		throw ValueError("Neither enthalpy nor pressure provided as the first parameter");
	}
}

double TTSESinglePhaseTableClass::evaluate_Trho(long iOutput, double T, double rho, double logrho)
{
	// Use Bicubic interpolation if requested
	if (this->mode == TTSE_MODE_BICUBIC){
		return interpolate_bicubic_Trho(iOutput,T,rho,logrho);
	}

	int i = (int)round((T-Tmin)/(Tmax-Tmin)*(NT-1));
	int j = (int)round((logrho-logrhomin)/logrhoratio);
	
	if (i<0 || i>(int)NT-1 || j<0 || j>(int)Nrho-1)
	{
		throw ValueError(format("Input to TTSE [T = %0.16g, rho = %0.16g] is out of range",T,rho));
	}

	// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
	// might be in the two-phase region which is not defined for single-phase table.  
	// Therefore, search around its neighbors for a better choice
	if (!ValidNumber(mu_Trho[i][j])){
		nearest_good_neighbor_Trho(&i,&j);
	}

	// Distances from the node
	double deltaT = T - this->T_Trho[i];
	double deltarho = rho - this->rho_Trho[j];
	
	switch (iOutput)
	{
	case iS:
		return s_Trho[i][j]+deltaT*dsdT_Trho[i][j]+deltarho*dsdrho_Trho[i][j]+0.5*deltaT*deltaT*d2sdT2_Trho[i][j]+0.5*deltarho*deltarho*d2sdrho2_Trho[i][j]+deltaT*deltarho*d2sdTdrho_Trho[i][j]; break;
	case iP:
		return p_Trho[i][j]+deltaT*dpdT_Trho[i][j]+deltarho*dpdrho_Trho[i][j]+0.5*deltaT*deltaT*d2pdT2_Trho[i][j]+0.5*deltarho*deltarho*d2pdrho2_Trho[i][j]+deltaT*deltarho*d2pdTdrho_Trho[i][j]; break;
	case iH:
		return h_Trho[i][j]+deltaT*dhdT_Trho[i][j]+deltarho*dhdrho_Trho[i][j]+0.5*deltaT*deltaT*d2hdT2_Trho[i][j]+0.5*deltarho*deltarho*d2hdrho2_Trho[i][j]+deltaT*deltarho*d2hdTdrho_Trho[i][j]; break;
	case iV:
		return mu_Trho[i][j]+deltaT*dmudT_Trho[i][j]+deltarho*dmudrho_Trho[i][j]+0.5*deltaT*deltaT*d2mudT2_Trho[i][j]+0.5*deltarho*deltarho*d2mudrho2_Trho[i][j]+deltaT*deltarho*d2mudTdrho_Trho[i][j]; break;
	case iL:
		return k_Trho[i][j]+deltaT*dkdT_Trho[i][j]+deltarho*dkdrho_Trho[i][j]+0.5*deltaT*deltaT*d2kdT2_Trho[i][j]+0.5*deltarho*deltarho*d2kdrho2_Trho[i][j]+deltaT*deltarho*d2kdTdrho_Trho[i][j]; break;
	default:
		throw ValueError(format("Output key value [%d] to evaluate is invalid",iOutput));
	}
}

double TTSESinglePhaseTableClass::evaluate_first_derivative(long iOF, long iWRT, long iCONSTANT, double p, double logp, double h)
{
	// Use Bicubic interpolation if requested
	if (mode == TTSE_MODE_BICUBIC){ 
		return bicubic_evaluate_first_derivative_ph(iOF, iWRT, iCONSTANT, p, logp, h);
	}

	int i = (int)round(((h-hmin)/(hmax-hmin)*(Nh-1)));
	int j = (int)round((logp-logpmin)/logpratio);

	if (i<0 || i>(int)Nh-1 || j<0 || j>(int)Np-1)
	{
		throw ValueError(format("Input to TTSE deriv [p = %0.16g, h = %0.16g] is out of range",p,h));
	}
	
	// If the value at i,j is too close to the saturation boundary, the nearest point i,j 
	// might be in the two-phase region which is not defined for single-phase table.  
	// Therefore, search around its neighbors for a better choice
	if (!ValidNumber(T[i][j]))
	{
		nearest_good_neighbor(&i,&j);
	}

	// Distances from the node
	double deltah = h-this->h[i];
	double deltap = p-this->p[j];
	
	// This is a first-order expansion of the derivative around the node point.
	//
	// Derivatives for constant p
	if (iOF == iT && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of T w.r.t. h for p constant (for cp, the constant-pressure specific heat)
		return dTdh[i][j]+deltah*d2Tdh2[i][j]+deltap*d2Tdhdp[i][j];
	}
	else if (iOF == iS && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of s w.r.t. h for p constant
		return dsdh[i][j]+deltah*d2sdh2[i][j]+deltap*d2sdhdp[i][j];
	}
	else if (iOF == iD && iWRT == iH && iCONSTANT == iP)
	{
		// Derivative of density w.r.t. h for p constant
		return drhodh[i][j]+deltah*d2rhodh2[i][j]+deltap*d2rhodhdp[i][j];
	}

	// Derivatives for constant h
	else if (iOF == iT && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of T w.r.t. p for h constant
		return dTdp[i][j]+deltap*d2Tdp2[i][j]+deltah*d2Tdhdp[i][j];
	}
	else if (iOF == iS && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of s w.r.t. p for h constant
		return dsdp[i][j]+deltap*d2sdp2[i][j]+deltah*d2sdhdp[i][j];
	}
	else if (iOF == iD && iWRT == iP && iCONSTANT == iH)
	{
		// Derivative of density w.r.t. p for h constant
		return drhodp[i][j]+deltap*d2rhodp2[i][j]+deltah*d2rhodhdp[i][j];
	}
	else{
		throw ValueError(format("Your inputs [%d,%d,%d,%g,%g] are invalid to evaluate_first_derivative",iOF,iWRT,iCONSTANT,p,h));
	}
}

TTSETwoPhaseTableClass::TTSETwoPhaseTableClass(Fluid *pFluid, double Q)
{
	this->pFluid = pFluid;
	this->Q = Q;
}
void TTSETwoPhaseTableClass::set_size(unsigned int N)
{
	this->N = N;
	
	// Seed the generator
	srand((unsigned int)time(NULL));

	// Resize all the arrays
	h.resize(N);
	p.resize(N);
	logp.resize(N);
	T.resize(N);
	dTdp.resize(N);
	d2Tdp2.resize(N);
	rho.resize(N);
	logrho.resize(N);
	drhodp.resize(N);
	d2rhodp2.resize(N);
	s.resize(N);
	dsdp.resize(N);
	d2sdp2.resize(N);
	h.resize(N);
	dhdp.resize(N);
	d2hdp2.resize(N);
}

double TTSETwoPhaseTableClass::build(double pmin, double pmax, TTSETwoPhaseTableClass *other)
{
	CoolPropStateClassSI CPS(pFluid);

	this->pmin = pmin;
	this->pmax = pmax;
	this->logpmin = log(pmin);
	this->logpmax = log(pmax);
	this->pratio = pow(pmax/pmin,1/((double)N-1));
	this->logpratio = log(pratio);

	double dlogp = (logpmax-logpmin)/(N-1);
	clock_t t1,t2;
	t1 = clock();
	// Logarithmic distribution of pressures
	for (unsigned int i = 0; i < N-1; i++)
	{
		// Calculate the pressure
		p[i] = exp(logpmin + i*dlogp);
		logp[i] = log(p[i]);
		// Update the class
		CPS.update(iP,p[i],iQ,Q);

		// Set the variables
		T[i] = CPS.T();
		dTdp[i] = CPS.dTdp_along_sat();
		d2Tdp2[i] = CPS.d2Tdp2_along_sat();
		h[i] = CPS.h();
		dhdp[i] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
		d2hdp2[i] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
		s[i] = CPS.s();
		dsdp[i] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
		d2sdp2[i] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
		rho[i] = CPS.rho();
		logrho[i] = log(CPS.rho());
		drhodp[i] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
		d2rhodp2[i] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

		// If other is provided
		if (other != NULL)
		{
			other->pmin = pmin;
			other->pmax = pmax;
			other->logpmin = log(pmin);
			other->logpmax = log(pmax);

			other->p[i] = this->p[i];
			other->logp[i] = log(this->p[i]);

			// Set the variables
			other->T[i] = CPS.T();
			other->dTdp[i] = CPS.dTdp_along_sat();
			other->d2Tdp2[i] = CPS.d2Tdp2_along_sat();
			other->h[i] = (other->Q>0.5) ? CPS.hV() : CPS.hL();
			other->dhdp[i] = (other->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
			other->d2hdp2[i] = (other->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
			other->s[i] = (other->Q>0.5) ? CPS.sV() : CPS.sL();
			other->dsdp[i] = (other->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
			other->d2sdp2[i] = (other->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
			other->rho[i] = (other->Q>0.5) ? CPS.rhoV() : CPS.rhoL();
			other->logrho[i] = log(other->rho[i]);
			other->drhodp[i] = (other->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
			other->d2rhodp2[i] = (other->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();
		}
	}
	// At the last point (the critical point)
	CPS.flag_SinglePhase = true; // Don't have it check the state or do a saturation call
	CPS.update(iT,CPS.pFluid->reduce.T+1e-12,iD,CPS.pFluid->reduce.rho+1e-12);
	p[N-1] = CPS.p();
	logp[N-1] = log(p[N-1]);
	T[N-1] = CPS.T();
	dTdp[N-1] = 1/CPS.dpdT_constrho();
	//d2Tdp2[N-1] = CPS.d2Tdp2_along_sat();
	h[N-1] = CPS.h();
	//dhdp[N-1] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
	//d2hdp2[N-1] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
	s[N-1] = CPS.s();
	//dsdp[N-1] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
	//d2sdp2[N-1] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
	rho[N-1] = CPS.rho();
	logrho[N-1] = log(CPS.rho());
	//drhodp[N-1] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
	//d2rhodp2[N-1] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();

	// If other is provided
	if (other != NULL)
	{
		other->p[N-1] = CPS.p();
		other->logp[N-1] = log(p[N-1]);
		other->T[N-1] = CPS.T();
		other->dTdp[N-1] = 1/CPS.dpdT_constrho();
		//other->d2Tdp2[N-1] = CPS.d2Tdp2_along_sat();
		other->h[N-1] = CPS.h();
		//other->dhdp[N-1] = (this->Q>0.5) ? CPS.dhdp_along_sat_vapor() : CPS.dhdp_along_sat_liquid();
		//other->d2hdp2[N-1] = (this->Q>0.5) ? CPS.d2hdp2_along_sat_vapor() : CPS.d2hdp2_along_sat_liquid();
		other->s[N-1] = CPS.s();
		//other->dsdp[N-1] = (this->Q>0.5) ? CPS.dsdp_along_sat_vapor() : CPS.dsdp_along_sat_liquid();
		//other->d2sdp2[N-1] = (this->Q>0.5) ? CPS.d2sdp2_along_sat_vapor() : CPS.d2sdp2_along_sat_liquid();
		other->rho[N-1] = CPS.rho();
		other->logrho[N-1] = log(CPS.rho());
		//other->drhodp[N-1] = (this->Q>0.5) ? CPS.drhodp_along_sat_vapor() : CPS.drhodp_along_sat_liquid();
		//other->d2rhodp2[N-1] = (this->Q>0.5) ? CPS.d2rhodp2_along_sat_vapor() : CPS.d2rhodp2_along_sat_liquid();
	}

	t2 = clock();
	std::cout << double(t2-t1)/CLOCKS_PER_SEC << " to build both two phase tables" << std::endl;
	return double(t2-t1)/CLOCKS_PER_SEC;
}

double TTSETwoPhaseTableClass::evaluate(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit below the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	// If the pressure is just barely above the critical pressure
	// or between the critical pressure and the next lowest point,
	// just just the next lowest point to avoid some of the 
	// derivatives that are infinite at the critical point
	if (i == (int)N || i == (int)N-1 ) i = N-2;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range[%g,%g], yielded index of: %d",p,pmin,pmax,i));
	}

	if (i == (int)N-2)
	{
		double y1, yc, dydp1;
		switch (iParam)
		{
		case iS:
			y1 = s[i]; yc = s[i+1]; dydp1 = dsdp[i]; break;
		case iT:
			y1 = T[i]; yc = T[i+1]; dydp1 = dTdp[i]; break;
		case iH:
			y1 = h[i]; yc = h[i+1]; dydp1 = dhdp[i]; break;
		case iD:
			y1 = rho[i]; yc = rho[i+1]; dydp1 = drhodp[i]; break;
		default:
			throw ValueError();
		}

		// Here we use interpolation based on a form y = a*x^(1/n) 
		// where we develop shifted coordinates
		// Y = y - yc (where y is s, h, rho, etc.)
		// X = pc - p

		double X = this->p[i+1] - p;
		double Y1 = y1 - yc;
		double X1 = this->p[i+1] - this->p[i];
		double dYdX1 = -dydp1;
		double n = Y1/(dYdX1*X1);
		double a = Y1/pow(X1,1/n);
		double Y = a*pow(X,1.0/n);
		return Y + yc;
	}
	else
	{
		// Spline interpolation http://en.wikipedia.org/wiki/Spline_interpolation since we
		// know the derivatives and the values at the bounding elements
		// Independent variable is logp
		// Dependent variable is varied (entropy, enthalpy, etc...)
		double x1 = this->logp[i];
		double x2 = this->logp[i+1];
		double t = (logp-x1)/(x2-x1);
		
		double y1, y2, k1, k2;
		switch (iParam)
		{
		case iS:
			y1 = this->s[i]; y2 = this->s[i+1]; k1 = this->p[i]*this->dsdp[i]; k2 = this->p[i+1]*this->dsdp[i+1]; break;
		case iT:
			y1 = this->T[i]; y2 = this->T[i+1]; k1 = this->p[i]*this->dTdp[i]; k2 = this->p[i+1]*this->dTdp[i+1]; break;
		case iH:
			y1 = this->h[i]; y2 = this->h[i+1]; k1 = this->p[i]*this->dhdp[i]; k2 = this->p[i+1]*this->dhdp[i+1]; break;
		case iD:
			/// log(p) v. log(rho) gives close to a line for most of the curve
			/// d(log(p))/d(log(p)) = 1/rho*d(rho)/d(log(p)) = p/rho*drho/dp
			y1 = this->logrho[i]; y2 = this->logrho[i+1]; k1 = this->p[i]/this->rho[i]*this->drhodp[i]; k2 = this->p[i+1]/this->rho[i+1]*this->drhodp[i+1]; break;
		default:
			throw ValueError();
		}
		
		double a = k1*(x2-x1)-(y2-y1);
		double b = -k2*(x2-x1)+(y2-y1);
		double y = (1-t)*y1+t*y2+t*(1-t)*(a*(1-t)+b*t);
		if (iParam == iD){ 
			return exp(y);
		}
		else{
			return y;
		}
	}
}
double TTSETwoPhaseTableClass::evaluate_T(double T)
{
	int L,R,M;
	double logp_spacing;

	logp_spacing = this->logp[2]-this->logp[1];

	// Do interval halving over the whole range to find the nearest temperature
	L = 0; R = N - 2; M = (L+R)/2;
	if (isbetween(this->T[N-2],pFluid->reduce.T,T))
	{
		// According to Matthis Thorade, dTdP|sat at the critical point is equal to dT/dP|rho evaluated at the 
		// critical temperature and density
		L = N-2;
		// Spline interpolation http://en.wikipedia.org/wiki/Spline_interpolation since we
		// know the derivatives and the values at the bounding elements
		// Independent variable is T
		// Dependent variable is logp
		double t = (T-this->T[L])/(this->T[L+1]-this->T[L]);
		double x1 = this->T[L];
		double x2 = this->T[L+1];
		double y1 = this->logp[L];
		double y2 = this->logp[L+1];
		// y is log(p); d(log(p))/dT = 1/p*(dp/dT) = 1/p/(dT/dp)
		double k1 = 1/this->p[L]/this->dTdp[L];
		double k2 = 1/this->p[L+1]/this->dTdp[L+1];
		double a = k1*(x2-x1)-(y2-y1);
		double b = -k2*(x2-x1)+(y2-y1);
		double logp = (1-t)*y1+t*y2+t*(1-t)*(a*(1-t)+b*t);
		return exp(logp);
	}
	else
	{
		while (R - L > 1)
		{
			if (T > this->T[M]){ 
				L=M; M=(L+R)/2; continue;
			}
			else{ 
				R=M; M=(L+R)/2; continue;
			}
		}
		// Spline interpolation http://en.wikipedia.org/wiki/Spline_interpolation since we
		// know the derivatives and the values at the bounding elements
		// Independent variable is T
		// Dependent variable is logp
		double t = (T-this->T[L])/(this->T[R]-this->T[L]);
		double x1 = this->T[L];
		double x2 = this->T[R];
		double y1 = this->logp[L];
		double y2 = this->logp[R];
		// y is log(p); d(log(p))/dT = 1/p*(dp/dT) = 1/p/(dT/dp)
		double k1 = 1/this->p[L]/this->dTdp[L];
		double k2 = 1/this->p[R]/this->dTdp[R];
		double a = k1*(x2-x1)-(y2-y1);
		double b = -k2*(x2-x1)+(y2-y1);
		double logp = (1-t)*y1+t*y2+t*(1-t)*(a*(1-t)+b*t);
		return exp(logp);
	}
	
}
double TTSETwoPhaseTableClass::evaluate_sat_derivative(long iParam, double p)
{
	double logp = log(p);
	int i = (int)round(((logp-logpmin)/(logpmax-logpmin)*(N-1)));
	// If the value is just a little bit out of the range, clip 
	// it back to the range of the LUT
	if (i == -1) i = 0;
	if (i == (int)N) i = N-1;
	// If it is really out of the range, throw an error
	if (i<0 || i>(int)N-1)
	{
		throw ValueError(format("p [%g] is out of range",p));
	}
		
	double pi = this->p[i];
	
	switch (iParam)
	{
	case iT:
		{
			/// First order expansion of dTdp around point of interest
			return dTdp[i]+(p-pi)*d2Tdp2[i];
		}
	case iH:
		{	
			/// First order expansion of dhdp around point of interest
			return dhdp[i]+(p-pi)*d2hdp2[i];
		}
	case iS:
		{	
			/// First order expansion of dsdp around point of interest
			return dsdp[i]+(p-pi)*d2sdp2[i];
		}
	case iD:
		{
			/// First order expansion of drhodp around point of interest
			return drhodp[i]+(p-pi)*d2rhodp2[i];
		}
	default:
		throw ValueError(format("Cannot use the key [%d] provided in evaluate_sat_derivative",iParam));
	}
}

double TTSETwoPhaseTableClass::evaluate_randomly(long iParam, unsigned int N)
{		
	clock_t t1,t2;
	t1 = clock();
	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;

		// Get the value from TTSE
		evaluate(iParam,p1);
	}
	t2 = clock();
	return (double)(t2-t1)/CLOCKS_PER_SEC/(double)N*1e6;
}


double TTSETwoPhaseTableClass::check_randomly(long iParam, unsigned int N, std::vector<double> *p, std::vector<double> *EOSv, std::vector<double> *TTSE)
{	
	double val=0;
	p->resize(N);
	EOSv->resize(N);
	TTSE->resize(N);
	
	CoolPropStateClassSI CPS(pFluid);

	for (unsigned int i = 0; i < N; i++)
	{
		double p1 = ((double)rand()/(double)RAND_MAX)*(pmax-pmin)+pmin;
		
		CPS.update(iP,p1,iQ,this->Q);
		double hEOS = CPS.h();
		double sEOS = CPS.s();
		double TEOS = CPS.T();
		double rhoEOS = CPS.rho();

		// Store the inputs
		(*p)[i] = p1;

		// Get the value from TTSE
		(*TTSE)[i] = evaluate(iParam,p1);
		
		// Get the value from EOS
		switch (iParam)
		{
		case iS: 
			(*EOSv)[i] = sEOS; break;
		case iT:
			(*EOSv)[i] = TEOS; break;
		case iH:
			(*EOSv)[i] = hEOS; break;
		case iD:
			(*EOSv)[i] = rhoEOS; break;
		default:
			throw ValueError();
		}
		
		std::cout << format("%g %g %g %g TTSE (p,EOS,TTSE, diff)\n",p1,(*EOSv)[i],(*TTSE)[i],((*EOSv)[i]-(*TTSE)[i])).c_str();
	}
	return val;
}

