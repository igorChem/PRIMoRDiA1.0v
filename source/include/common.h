//commom.h
//Declarations for functions and source file with constant values. 

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

#ifndef COMMON
#define COMMON
//--------------------------------------------------------------------------------
//Including header from the c++
#include <iostream>
#include <string>
#include <vector> 
#include <memory>
#include <experimental/filesystem>
//--------------------------------------------------------------------------------
//Including header from PRIMoRDiA library. 
#include "log_class.h"
#include "Itimer.h"

//===========================================================
// GLOBAL VARIABLES for internal usage: DECLARION
//===========================================================
extern bool dos;
extern bool extra_RD;
extern bool pymol_script;
extern double energy_crit; 
extern bool M_R;

extern unsigned int NP; // global  holding the maximum number of openMP threads to be used.
extern Itimer chronometer; // global object that returns total wall time of  execution. 
extern std::unique_ptr<Ilog> m_log; // object that writes to a log file and/or outputs messages to the console.

//===========================================================
// GLOBAL VARIABLES for main function usage DECLARION
//===========================================================

extern bool M_verbose;
extern bool M_logfile;

//===========================================================
// Auxiliary Functions definitions
//------------------------------------------------------------------------------------------
/**
 * @fn get_atom_mass
 * @brief Function to get the atomic mass from array.
 * @param Element symbol.
 * @return Atomic mass value.
 */
float get_atom_mass(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_number
 * @brief Function to get the atomic number from array.
 * @param Element symbol.
 * @return Atomic number.
 */
int get_atomic_number(std::string sym);

//-----------------------------------------------------------------------------------------
/**
 * @fn get_atomic_symbol
 * @brief Function to get the atomic symbol from array.
 * @param Atomic nuber.
 * @return Atomic symbol.
 */
std::string get_atomic_symbol(int i);

//---------------------------------------------------------------------------------------
/**
 * @fn factorial
 * @brief Calculate the factorial of an integer value.
 * @param Value to get the factorial from.
 * @return Factorial result 
 */
int factorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn double factorial
 * @brief Calculate the semi-factorial/double factorial of an integer value.
 * @param Value to get the factorial from.
 * @return Semi-factorial of the integer passed. 
 */
int doublefactorial(int n);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to a std string object.
 * @param Value to be converted.
 * @return String object.
 */
std::string int_to_string(int val);

//---------------------------------------------------------------------------------------
/**
 * @fn int_to_string
 * @brief Convert integer value to a std string.
 * @param Integer to be converted to string.
 * @return String from integer passed.
 */
std::string double_to_string(double val);

//---------------------------------------------------------------------------------------
/**
 * @fn IF_file
 * @brief Test if file can be open.
 * @param Name of the file.
 * @return If the file can be open.
 */
bool IF_file(const char* name); 
bool IF_file(std::experimental::filesystem::path& name); 

//---------------------------------------------------------------------------------------
/**
 * @fn check_file_ext
 * @brief Test if file the file extension is the one that is required.
 * @param Extension to be tested.
 * @param Name of the file.
 * @return If the files has the right extension.
 */
bool check_file_ext(std::string ext,const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn remove_extension
 * @brief return a string with the filename without the extension part 
 * @param constant char* containing the file name.
 * @return string with the filename without the extension.
 */
std::string remove_extension(const char* file_name);

//---------------------------------------------------------------------------------------
/**
 * @fn change_extension
 * @brief Return a string with a file name with a new extension
 * @param char pointer constant with the file name
 * @param string with the new extension
 * @return string of the file name without the extension
 */
std::string change_extension(const char* file_name,std::string new_ext);

//---------------------------------------------------------------------------------------
/**
 * @brief Rename files in current directory.
 * @param Path to the file.
 * @param New name for the file.
 * @return None.
 */
void rename_file(const char* file_name,std::string new_file_name);

//---------------------------------------------------------------------------------------
/**
 * @brief Extract only the name component of a path
 * @param File path.
 * @return Name component of the file path. 
 */
std::string get_file_name(const char* path_name);
//---------------------------------------------------------------------------------------
/**
 * @brief Convert double in string format from file with "D" to "E"
 * @param Double to be converted
 * @return Double with the converted values
 */
double D_E_conv(std::string sc_not);

//---------------------------------------------------------------------------------------
/**
 * @brief  Get the average mean from a vector.
 * @param STL vector of doubles
 * @return Average mean
 */
double mean_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief Get the maximum value from a vector.
 * @param STL vector of doubles
 * @return Maximum value
 */
double max_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief Get the minimum from a vector.
 * @param STL vector of doubles
 * @return Minimum value
 */
double min_dvec(std::vector<double>& vec);
//---------------------------------------------------------------------------------------
/**
 * @brief Get the sum from a vector.
 * @param STL vector of doubles
 * @return Sum.
 */
double sum_dvec(std::vector<double>& vec);

//---------------------------------------------------------------------------------------
/**
 * @brief Substring a line.
 * @param Line to be subdued.
 * @param Initinal character index
 * @param Final character index
 * @return Result.
 */
std::string str_array(std::string& line, int in, int fin);

#endif
//================================================================================
//END OF FILE
//================================================================================