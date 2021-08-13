// header file for local_rd class
// local_rd.h

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

#ifndef LOCAL_RD
#define LOCAL_RD

//include statements from c++ library
#include <iostream>
#include <string> 
#include <vector>
//include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Icube.h"

//foward declarations 
class Iatom;
class Imolecule;
class Icube;
class global_rd;
//=================================================================================================
/**
 * This class is meant to represent a collection of local reactivity descriptors for a given molecular system,
 * in volumetric representation hold by Icube objects or in condensed form hold by STL vectors of doubles. 
 * The constructors of the current class instantiates objects passing necessary information to make the calculation
 * ans reporting of the local RD using the Koopman approximation or finite differences method.
 * @class local_rd
 * @author Igor Barden Grillo
 * @date 20/03/18
 * @file local_rd.h
 * @brief Class to calculate local reactivity descriptors from scalar fields and atomic charges from 
 * quantum chemistry info provided in the programs outputs. 
 */
class local_rd {
	public:
		//member variables
		std::string name;
		bool FD;
		bool LH;
		bool band;
		bool TFD;
		int charge;
		std::vector<Icube> lrds;
		std::vector<std::string> rd_names;
		//constructors/destructor
		local_rd();
		local_rd(const Icube& HOmo, const Icube& LUmo);
		local_rd(const Icube& elec_dens, const Icube& HOmo, const Icube& LUmo);
		local_rd(const Icube& elecDens, const Icube& cationDens, const Icube& anionDens,int chg);
		local_rd(const local_rd& lrd_rhs);
		local_rd& operator=(const local_rd& lrd_rhs);
		local_rd(local_rd&& lrd_rhs) noexcept;
		local_rd& operator=(local_rd&& lrd_rhs) noexcept;
		~local_rd();
		//member functions
		friend local_rd operator-(const local_rd& lrd_lhs,const local_rd& lrd_rhs);	
		void calculate_fukui_Band(const Icube& homo_b, const Icube& lumo_b);
		void calculate_RD(const global_rd& grd);
		void calculate_Fukui_potential();
		void calculate_hardness(const global_rd& grd);
		void calculate_MEP(const Imolecule& mol);
		void write_LRD();
};

#endif
//================================================================================
//END OF FILE
//================================================================================