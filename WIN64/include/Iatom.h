//Iatom.h

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

#ifndef IATOM
#define IATOM
//----------------------------------------------------------
//including c++ headers
#include <string>
#include <vector>
//----------------------------------------------------------
//including PRIMoRDiA headers
#include "../include/common.h"
//----------------------------------------------------------
class Iaorbital; // foward declaration

//=================================================================================================
/**
 * This class is meant to hold values for atomic properties, also manipulate and modify, from quantum
 * chemistry programs outputs
 * 
 * @class Iatom
 * @author Igor Barden Grillo - barden.igor@gmail.com
 * @date 20/03/18
 * @file Iatom.h
 * @brief A Iatom Class to instatiate atom representation object to Quantum Mechanical properties calculations. 
 */
class Iatom {
	public:
		// member variables
		double xcoord; // x-axis atomic coordinate
		double ycoord;
		double zcoord;
		std::string element;
		float charge;
		float atomic_mass;
		unsigned int atomicN;
		unsigned int norb;
		double wdw_volume;
		std::vector<Iaorbital> orbitals;
		// constructors
		Iatom();
		Iatom(double x,double y,double z, std::string typ); 
		Iatom(const Iatom& rhs_atom); 
		Iatom& operator=(const Iatom& rhs_atom);
		Iatom(Iatom&& rhs_atom) noexcept;
		Iatom& operator=(Iatom&& rhs_atom) noexcept;
		~Iatom();
		// member functions 
		void print();
		void set_type(const std::string typ);
		void set_coord(double x, double y, double z);
		void add_orbital(int level,const std::string sym, double coef);
		void add_orbital(Iaorbital orb);
		void unit_test();
};

#endif
//================================================================================
//END OF FILE
//================================================================================