//header file for QM output file parser
// QMparser.h

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

#ifndef QMPARSER
#define QMPARSER
//---------------------------------------------
//Including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <memory>
//---------------------------------------------
#include "../include/common.h"
//---------------------------------------------
class Imolecule;
/**
 * @class QMparser
 * @author igor Barden Grillo
 * @date 24/03/18
 * @file QMparser.h
 * @brief Class to read and stores files from output of quantum chemical calculations programs
 * to generate objects with info to calculate reactivity descriptors.
 * 
 * This class instatiates a Ibuffer object to load all the lines of a given file from a quantum
 * chemistry output calculation and parse. For now the member functions to parse file of gamess and mopac
 * are implemented, but the class make available the needed tools to parse formated file and store the chemical
 * relevant information. 
 */
class QMparser{
	public:
		//member variables
		std::string program;
		bool parsed;
		const char* name_f;
		//constructors/destructor
		QMparser();
		QMparser(const char* file_name, std::string Program);
		QMparser(const QMparser& rhs_QMp) = delete;
		QMparser& operator=(const QMparser& rhs_QMp) = delete;
		~QMparser();
		//member functions
		Imolecule get_molecule();
	
};
 
#endif

//================================================================================
//END OF FILE
//================================================================================