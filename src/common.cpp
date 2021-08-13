//commmon.cpp
// Source file for common utilities, 
 
/*********************************************************************/
/* This source code file is part of PRIMoRDiA software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

//--------------------------------------------------------------------------------
//Including header from the c++
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <omp.h>
#include <fstream>
#include <experimental/filesystem>

#include "../include/common.h"

using std::string;
namespace fs = std::experimental::filesystem;

bool dos			= false;
bool extra_RD		= false;
bool pymol_script	= false;
double energy_crit	= 1;
bool M_R			= false;
bool comp_H			= false;
/*********************************************************************************/
unsigned int NP		= omp_get_max_threads();
Itimer chronometer;
std::unique_ptr<Ilog> m_log ( new Ilog() );
/*********************************************************************************/
bool M_verbose	= false;
bool M_logfile	= false;
/************************************************************************************************************/
string atomType[] = {
					 "H",																			      "He",
					 "Li","Be",											  		 "B", "C", "N", "O", "F", "Ne",
					 "Na","Mg",										  			 "Al","Si","P", "S", "Cl","Ar",
					 "K","Ca", "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
					 "Rb","Sr", "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te", "I","Xe",
					 "Cs","Ba",
								"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
								"Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
					 "Fr","Ra",
							  "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bh","Bk","Cf","Es","Fm","Md","No","Lr"
					 };
/************************************************************************************************************/
double atomMass[] = {
					1.00794,																				   		 				 			4.0026,
					6.941,9.012182,																	10.811,	12.0107,14.0067,15.9994,18.9984,	20.1797,
					22.9897,24.305,																	26.981,28.0855,30.973762,32.065,35.453,39.948,39.948,
					40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.6934,58.9332,58.693,63.546,65.39,69.723,72.64,74.9216,78.96,79.904,83.8,
					85.467,87.62l,88.9059,91.224,92.9064,95.94,98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,126.9045,127.6,131.293,
					132.9055,137.327,
									138.9055,140.116,140.9077,144.24,145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,174.967,
									178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,204.3833,207.2,208.9804,209,210,222,
					223,226,
							227,231.0359,232.0381,237,238.0289,243,244,247,247,251,252,257,258,259,261,262
					};					
/************************************************************************************************************/
double wdw_radius[] = {
						1.1,                                                                                  1.4,
						1.82,1.55,                                                  1.92,1.70,1.55,1.52,1.47,1.54,
						2.27,1.73,                                                  1.84,2.10,1.80,1.80,1.75,1.88,
						2.75,2.31,2.15,2.11,2.07,2.06,2.05,2.04,2.00,1.97,1.96,2.01,1.87,2.11,1.85,1.90,1.85,2.02,
						3.03,2.49,2.32,2.23,2.46,2.17,2.16,2.13,2.10,2.10,2.11,2.18,1.93,2.17,2.06,2.06,1.98,2.16,
						3.43,2.68,
								2.43,2.42,2.40,2.39,2.38,2.36,2.35,2.34,2.33,2.31,2.30,2.29,2.27,2.26,2.24,
								2.23,2.22,2.18,2.16,2.16,2.13,2.43,2.14,2.23,1.96,2.02,2.07,1.97,2.02,2.20,
						3.48,2.43,
								2.47,2.45,2.43,2.41,2.39,2.43,2.44,2.45,2.46,2.44,2.45,2.45,2.45,2.46,2.46,2.46
					};
					
double wdw_volume[] = {
						5.5758,																																	11.4952,																
						25.2549,15.6001,																				29.6507,20.5815,15.6001,14.7117,13.3071,15.3001,
						49.0014,21.6905,																				26.0966,38.7962,24.4314,24.4314,22.4515,27.8359,							
						87.1223,51.6377,41.6338,39.3531,37.1571,36.6212,36.0905,35.5649,33.5136,32.0280,31.5427,34.0188,27.3940,39.3531,26.5244,28.7337,26.5244,34.5291,
						116.5357,64.6739,52.3112,46.4564,62.3644,42.8066,42.2175,40.4827,38.7962,38.7962,39.3531,43.4011,30.1164,42.8066,36.6212,36.6212,32.5182,42.2175,
						169.0493,80.6372,									  
									60.1104,59.3714,57.9115,57.1906,56.4757,55.0639,54.3669,53.6758,52.9906,51.6377,50.9700,50.3081,49.0014,48.3567,47.0842,
									46.4564,45.8342,43.4011,42.2175,42.2175,40.4827,60.1104,41.0556,46.4564,31.5427,34.5291,37.1571,32.0280,34.5291,44.6066,
						176.5504,60.1104,									   
									63.1280,61.6069,60.1104,58.6384,57.1906,60.1104,60.8556,61.6069,62.3644,60.8556,61.6069,61.6069,61.6069,62.3644,62.3644,62.3644
					};				
					
/************************************************************************************************************/
std::string residue_3Lname[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
"PHE","PRO","SER","THR","TRP","TYR","VAL"};

std::string residue_1Lname[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};
/************************************************************************************************************/
//Zamyatnin, A.A., Protein volume in solution, Prog. Biophys. Mol. Biol., 24:107-123 (1972), PMID: 4566650.			
/************************************************************************************************************/
int factorial(int n){  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;}
/***********************************************************************************/
int doublefactorial(int n){  return (n == 1 || n == 0) ? 1 : factorial(n - 2) * n;}
/***********************************************************************************/
float get_atom_mass(string sym){
	int indx = 0;
	for (unsigned int i=0;i<103;i++) if ( sym == atomType[i] ) indx = i;
	return atomMass[indx];
}
/*********************************************************************************/
int get_atomic_number(string sym){
	int atomic_num;
	for (unsigned int i=0;i<103;i++) if ( sym == atomType[i] ) atomic_num = i+1;
	return atomic_num;
}
/*********************************************************************************/
string get_atomic_symbol(int i){ return atomType[i-1]; }
/*********************************************************************************/
double get_wdw_volume(int i){ return wdw_volume[i-1]; }
/*********************************************************************************/
bool IF_file(const char* name){
	fs::path file_name(name);
	return fs::exists(file_name);
}
/*********************************************************************************/
bool IF_file(fs::path& name){
	return fs::exists(name);
}
/*********************************************************************************/
bool check_file_ext(string ext,const char* file_name){
	fs::path f_name(file_name);
	if ( f_name.extension() == ext ) return true;
	else return false;
}
/*********************************************************************************/
string get_file_name(const char* path_name){
	fs::path f_name(path_name);
	string resultS( f_name.filename() );
	return resultS;
}
/*********************************************************************************/
string get_file_ext(const char* path_name){
	fs::path f_name(path_name);
	string resultS( f_name.extension() );
	return resultS;
}
/*********************************************************************************/
string remove_extension(const char* file_name){
	string file_name_string  = file_name;
	int point = 0;
	char dot = '.';
	for(unsigned int i = 0;i<file_name_string.size();i++){
		char character = file_name_string[i];
		if ( character == dot){	point = i; }
	}
	int pos = file_name_string.size()-point;
	return file_name_string.substr(0,file_name_string.size()-pos);
}
/*********************************************************************************/
string change_extension(const char* file_name,string new_ext){
	string name_wth_ext = remove_extension(file_name);
	name_wth_ext += new_ext;
	return name_wth_ext;
}
/*********************************************************************************/
void rename_file(const char* file_name,std::string new_file_name){
	fs::path f_name(file_name);
	fs::path nf_name( fs::current_path() / new_file_name );
	fs::rename(f_name,nf_name);
}
/********************************************************************************/
double D_E_conv(string sc_not){
	std::replace(sc_not.begin(),sc_not.end(),'D','E');
	double result = std::stod(sc_not);
	return result;
}
/********************************************************************************/
double mean_dvec(std::vector<double>& vec){
	double result = 0.0;
	for(unsigned int i=0;i<vec.size();i++){
		result += vec[i];
	}	
	return result/vec.size();
}
/********************************************************************************/
double max_dvec(std::vector<double>& vec){
	std::sort( vec.begin(), vec.end() );
	if ( vec.size() > 0) return vec[vec.size()-1];
	else return 0.0;
}
/********************************************************************************/
double min_dvec(std::vector<double>& vec){
	std::sort( vec.begin(), vec.end() );
	if ( vec.size() > 0) return vec[0];
	else return 0.0;
}
/********************************************************************************/
double sum_dvec(std::vector<double>& vec){
	double result = 0.0;
	for(unsigned int i=0;i<vec.size();i++){
		result += vec[i];
	}	
	return result;
}
/********************************************************************************/
std::vector<double> norm_dvec(std::vector<double>& vec, double size){
	std::vector<double> norma(vec);
	double max = sum_dvec(vec);
	for(unsigned int i=0; i<vec.size(); i++ ){
		vec[i] /= max;
		vec[i] *= size;
	}
}
/********************************************************************************/
string str_array(std::string& line, int in, int fin){
	string result = "";
	result.resize(fin-in+2);
	int cnt = 0;
	for(int i=in;i<fin;i++){ 
		result[cnt++] = line[i];
	}
	return result;
}
//================================================================================
//END OF FILE
//================================================================================