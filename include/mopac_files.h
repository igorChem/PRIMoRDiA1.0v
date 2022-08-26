//header file for QM output file parser
// mopac_files.h

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

#ifndef MOPACFILES
#define MOPACFILES

//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <memory>
//Including PRIMoRDiA headers
#include "../include/common.h"

class Imolecule; //foward declaration
class mopac_files{
	public:
		//member variables
		bool LMO;
		bool RHF;
		double f_chg;
		const char* name_f;
		bool is_open;
		std::string type;
		Imolecule molecule;
		//constructors/destructor
		mopac_files();
		mopac_files(const char* file_name);
		mopac_files(const mopac_files& rhs_mop) = delete;
		mopac_files& operator=(const mopac_files& rhs_mop) = delete;
		~mopac_files();
		//member functions
		void parse_aux();
		void parse_out();
		void parse_mgf();
		void get_overlap_m();
		void get_mo(bool beta);
		void get_mo_energies(bool beta);
		void get_mo_occupancies(bool beta);
};

#endif
//================================================================================
//END OF FILE
//================================================================================