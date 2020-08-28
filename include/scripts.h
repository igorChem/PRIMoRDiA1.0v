//scripts.h

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

#ifndef SCRIPTS
#define SCRIPTS
//--------------------------------------------------------------------------
#include <string>
//--------------------------------------------------------------------------
/**
 * @class scripts
 * @author Igor Barden Grillo
 * @date 15/03/20
 * @file scripts.h
 * @brief Class to handle and write scripts for pos-processment by other programs
 */
class scripts{
	public:
		//member variables
		const char* file_name;
		std::string s_type;
		//constructors/destructor
		scripts();
		scripts(const char* Nm);
		scripts(const scripts& rhs) = delete;
		scripts& operator=(const scripts& rhs) = delete;
		~scripts();
		//member functions
		void write_r_dos(std::vector<double>& energies);
		void write_pymol();
};
#endif
//================================================================================
//END OF FILE
//================================================================================