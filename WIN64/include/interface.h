//inerface.h

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

#ifndef INTERFACE 
#define INTERFACE

//--------------------------------------------
#include <iostream>
#include <cstring>
#include <string>
#include <vector>

//------------------------------------------
using std::cout;
using std::endl;
using std::stoi; 
using std::string;

//==================================================
/**
 * @class interface
 * @author Igor Barden Grillo
 * @date 14/03/20
 * @file interface.h
 * @brief Interface functions of the application.
 */
class interface{
	public:
		//member functions
		int m_argc;
		std::vector<std::string> m_argv;
		std::string runtyp;	
		//constructors/destructor
		interface();
		interface(int argc, char* argv[] );
		interface(const interface& rhs) = delete;
		interface& operator=(const interface& rhs) = delete;
		~interface();
		//member functions
		void run();
		void MO_cube();
		void ED_cube();
		void Comp_cube();
		void write_help();
		void test_run();
		void write_input();
		void print_options();
		void pos_traj_res(int* res);
};

#endif
//================================================================================
//END OF FILE
//================================================================================