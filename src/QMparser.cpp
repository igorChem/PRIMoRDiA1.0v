//Including c++ headers
//QMparser.cpp

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

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <algorithm>  

//Including PRIMoRDiA headers
//-------------------------------------------------------
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Imolecule.h"
#include "../include/QMparser.h"
#include "../include/gamess_files.h"
#include "../include/orca_files.h"
#include "../include/gaussian_files.h"
#include "../include/mopac_files.h"
//-------------------------------------------------------
// Aliases for standard c++ scope functions
using std::cout;
using std::endl;
using std::string;
/************************************************************************************/
QMparser::QMparser()	:
	program("none")		,
	name_f("none")		,
	parsed(false)		{
}
/************************************************************************************/
QMparser::QMparser(const char* file_name,
								string Program			):
	program(Program)									,
	name_f(file_name)									{
}
/************************************************************************************/
Imolecule QMparser::get_molecule(){
	Imolecule empty_molecule;
	empty_molecule.name = "empty";
	double initi_time = omp_get_wtime();
	if ( program == "mopac" ){
		mopac_files file_obj(name_f);
		if ( file_obj.is_open ){
			if ( file_obj.type == "AUX" ){
				file_obj.parse_aux();
			}else if ( file_obj.type == "OUT" ){
				file_obj.parse_out();
				file_obj.parse_mgf();
			}
			else if ( file_obj.type == "MGF" ){
				file_obj.parse_mgf();
			}
			double fin_time = omp_get_wtime() - initi_time;
			m_log->input_message("Time to read QM file: ");
			m_log->input_message(fin_time);
			return file_obj.molecule;
		}else{
			m_log->write_error("File file not open! Ending without succes the parsing process!");
			return empty_molecule;
		}
	}else if (program == "gamess"){
		gamess_files file_obj(name_f);
		if ( file_obj.is_open ){
			file_obj.parse_log();
			return file_obj.molecule;
		}else{
			m_log->write_error("File not open! Ending without succes the parsing process!");
			return empty_molecule;
		}
	}else if (program == "orca"){
		orca_files file_obj(name_f);
		if ( file_obj.is_open ){
			file_obj.parse_out();
			return file_obj.molecule;
		}else{
			m_log->write_error("File not open! Ending without succes the parsing process!");
			return empty_molecule;
		}
	}else if (program == "gaussian"){
		gaussian_files file_obj(name_f);
		if ( file_obj.is_open ){
			file_obj.parse_fchk();
			file_obj.get_overlap_m();
			return file_obj.molecule;
		}else{
			m_log->write_error("File not open! Ending without succes the parsing process!");
			return empty_molecule;
		}
	}else{
		cout << "Program keyword not recognized!" << endl;
		m_log->write_warning("Parsing process end without valid molecular information stored!");
		return empty_molecule;
	}
	
}
/**********************************************************************************/
QMparser::~QMparser(){}
//================================================================================
//END OF FILE
//================================================================================