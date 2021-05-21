//orca_files.h
//header file for ORCA QM output file parser

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

#ifndef ORCAFILES
#define ORCAFILES

//Including c++ headers
#include <string>
#include <sstream>
#include <cstring>
#include <memory>
// PRIMoRDiA headers
#include "../include/common.h"

class Ibuffer; //foward declarations
class Imolecule;

//====================================================
/**
 * @class orca_files
 * @author Igor Barden Grillo
 * @date 10/03/20
 * @file orca_files.h
 * @brief Class to define functions to read and parse output files produced by ORCA calculations.
 */
class orca_files{
	public:
		//member variables
		const char* name_f;
		Imolecule molecule;
		bool is_open;
		//costructors/destructor
		orca_files();
		~orca_files();
		orca_files(const char* file_name);
		orca_files(const orca_files& rhs) = delete;
		orca_files& operator=(const orca_files& rhs) = delete;
		//member functions
		void parse_out();
		void get_overlap(int ov_in, int ov_fin);
		
};
//===================================================
/**
 * @class basis_orca
 * @author Igor Barden Grillo.
 * @date 10/03/20
 * @file orca_files.h
 * @brief Class to define storage objects to hold and deal with atomic orbitals extracted from ORCA output files.
 */
class basis_orca{
	public:
		//member variables.
		std::string element_type;
		std::vector<double> shell_size;
		std::vector<std::string> shell_sym;	
		std::vector<double> coefficients;
		std::vector<double> exp;
		// constructors/destructors
		basis_orca(){};
};

#endif

//================================================================================
//END OF FILE
//================================================================================