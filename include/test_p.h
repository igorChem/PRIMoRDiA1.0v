//test_p.h

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

#ifndef TEST_PRIMORDIA
#define TEST_PRIMORDIA

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "../include/common.h"

//==================================================================
/**
 * @class test_p
 * @author barden
 * @date 17/04/18
 * @file test_p.h
 * @brief Class to instantiate a object with a collection of memeber functions for PRIMoRDiA-libs unit testing.
 */
class test_p {
	public:	
		//member varuables
		std::vector<std::string> mopac_aux;
		std::vector<std::string> mopac_out;
		std::vector<std::string> mopac_mgf;
		std::vector<std::string> orca_out;
		std::vector<std::string> gauss_fchk;
		std::vector<std::string> gauss_log;
		std::vector<std::string> gamess_log;
		std::vector<std::string> molden;
		std::vector<std::string> pdb_files;
		//constructors/destructor
		test_p();
		~test_p();
		//member functions
		void test_primordia_1();
		void test_primordia_2();
		void test_primordia_3();
		void test_int_molden();
};

#endif
//================================================================================
//END OF FILE
//================================================================================