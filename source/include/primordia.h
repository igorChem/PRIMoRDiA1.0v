//primordia.h

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

#ifndef PRIMORDIA 
#define PRIMORDIA

//C++ headers files
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

//foward declarations
class Imolecule;
class Icube;
class global_rd;
class local_rd;
class local_rd_cnd;

//============================================================================================================
/**
 * These class wraps all the functionalities of this libraries and make available in different ways of initializations
 * in respect of the type of the data passed, the type of approxiamtion calculus and representation of the reactivity
 * descriptors.
 * @class primordia
 * @author Igor Barden Grillo
 * @date 11/04/18
 * @file primordia.h
 * @brief Interface class for calculate the reactivity descriptors based on conceptual DFT.
 */
class primordia {
	public:
		//member varibles
		std::string name;
		unsigned int band; 
		std::unique_ptr<global_rd> grd;
		std::unique_ptr<local_rd> lrdVol;
		std::unique_ptr<local_rd_cnd> lrdCnd;
		//constructors/destructor
		primordia();
		primordia(const primordia& pr_rhs);
		primordia& operator=(const primordia& pr_rhs);
		primordia(primordia&& pr_rhs) noexcept;
		primordia& operator=(primordia&& pr_rhs) noexcept;
		~primordia();
		//member function
		friend primordia operator-(const primordia& pr_lhs,const primordia& pr_rhs);
		void init_FOA(const char* file_neutro,int gridN,std::string loc_hard,bool mep, std::string Program);
		void init_FD(const char* file_neutro,const char* file_cation,const char* file_anion, int grdN, int charge,bool mep,std::string loc_hard, std::string Program);			
		void init_protein_RD(const char* file_neutro,std::string locHardness,int grdN,int bandgap,double* ref_atom,int size,const char* pdb, bool mep , std::string bt, std::string Program);
		void init_QS_KA(Imolecule& mol, int gridN);
		void init_QS_FD( Imolecule& mol1, Imolecule& mol2, Imolecule& mol3, int charge ,int gridN);
		void write_pymol(const std::string& pdb, bool fixed);
};

#endif
//================================================================================
//END OF FILE
//================================================================================