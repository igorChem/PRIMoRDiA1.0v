//header file for GAMESS QM output file parser
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

#ifndef GAMESSFILES
#define GAMESSFILES

//Including c++ headers
#include <string>
#include <sstream>
#include <cstring>
#include <memory>
// PRIMoRDiA headers
#include "../include/common.h"
//--------------------------
class Iline; //fowards declarations 
class Ibuffer;
class Imolecule;
//--------------------------
/**
 * @class gamess_files
 * @author Igor Barden Grillo
 * @date 10/03/20
 * @file gamess_files.h
 * @brief Class where functions are defined to read and parse output files produced by GAMESS calculations
 */
class gamess_files{
	public:
		//member variables
		const char* name_f;
		bool is_open;
		std::unique_ptr<Ibuffer> buffer;
		std::unique_ptr<Imolecule> molecule;
		// constructors/destructor
		gamess_files();
		gamess_files(const char* file_name);
		gamess_files( const gamess_files& rhs) = delete;
		gamess_files& operator=( const gamess_files& rhs) = delete;
		~gamess_files();
		//member function
		void parse_log();
};

#endif
//================================================================================
//END OF FILE
//================================================================================