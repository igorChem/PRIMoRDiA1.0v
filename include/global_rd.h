//global_rd.h
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
// header file for global_rd.cpp source file 
// class for global reactivity descriptors 

#ifndef GLOBAL_RD
#define GLOBAL_RD
//include C++ headers
#include <iostream>
#include <string> 
#include <fstream>
//include our library header 
#include "../include/common.h"

//===================================================================================
/**
 * Class to represent objects that hold necessary information to and calculate the global reactivity
 * descriptors following the Conceptual Density Functional Theory.  
 * Contatc: barden.igor@gmail.com
 * @class global_rd
 * @author Igor barden Grillo
 * @date 20/03/18
 * @file global_rd.h
 * @brief Instantiate objects with quantum chemistry info extracted from other programs calculation to get the
 * global reactivity descriptors for the system. * 
 */
class global_rd{
	public:
		//member variables
		std::string name;
		std::vector<std::string> rd_names;
		std::vector<std::string> rd_abrev;
		std::vector<double> grds;
		bool KA;
		bool DF;
		//constructors/destructor
		global_rd();
		global_rd(const Imolecule& mol);
		global_rd(const Imolecule& mol_neutro,const Imolecule& mol_cation,const Imolecule& mol_anion);
		global_rd(const global_rd& rd_rhs);
		global_rd& operator=(const global_rd& rd_rhs);
		global_rd(global_rd&& rd_rhs) noexcept;
		global_rd& operator=(global_rd&& rd_rhs) noexcept;
		~global_rd();
		//member functions
		friend global_rd operator-(const global_rd& lhs_grd,const global_rd& rhs_grd);
		void calculate_rd();
		void print_rd();
		void write_rd();
};

#endif 
//================================================================================
//END OF FILE
//================================================================================