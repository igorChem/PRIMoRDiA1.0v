//header file for GAUSSIAN QM output file parser
// gamess_files.h

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

#ifndef GAUSSFILES
#define GAUSSFILES

//Including c++ headers
#include <string>
#include <cstring>
#include <memory>
// PRIMoRDiA headers
#include "../include/common.h"

class Ibuffer;
class Imolecule;

//====================================================
/**
 * @class orca_files
 * @author igor
 * @date 10/03/20
 * @file orca_files.h
 * @brief 
 */
class gaussian_files{
	public:
		//member varibales
		const char* name_f;
		std::unique_ptr<Ibuffer> buffer;
		std::unique_ptr<Imolecule> molecule;
		bool is_open;
		//constructors/destructor
		gaussian_files();
		~gaussian_files();
		gaussian_files(const char* file_name);
		gaussian_files(const gaussian_files& rhs) = delete;
		gaussian_files& operator=(const gaussian_files& rhs) = delete;
		//member variables
		void parse_fchk();
		void get_overlap_m();
};

#endif
//================================================================================
//END OF FILE
//================================================================================