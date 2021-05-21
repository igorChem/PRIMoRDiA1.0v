// header file for local_rd_cnd class
// local_rd_cnd.h

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
//-------------------------------------------------------------------------------
#ifndef LOCAL_RD_CND
#define LOCAL_RD_CND
//----------------------------------------------------------------------------------
#include <iostream>
#include <string> 
#include <vector>
#include <cmath>
//include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Imolecule.h"
//--------------------------------------------------------------------------------------
class Iatom;
//class Imolecule;
class QMdriver;
class global_rd;
class Iprotein;
class protein_lrd;

//===============================================================
class local_rd_cnd{
	public:
		//member variables
		std::string name;
		bool ed;
		bool mep;
		bool finite_diff;
		int charge;
		Imolecule molecule;
		std::vector<double> EAS;
		std::vector<double> NAS;
		std::vector<double> RAS;
		std::vector<double> dual;
		std::vector<double> Softness_Dual;
		std::vector<double> Hyper_softness;
		std::vector<double> hardness_A;
		std::vector<double> hardness_B;
		std::vector<double> hardness_C;
		std::vector<double> hardness_D;
		std::vector<double> multifilic;
		std::vector<double> electrophilicity;
		std::vector<double> fukushima;
		std::vector<double> electron_density;
		//constructors/destructor
		local_rd_cnd();
		local_rd_cnd(Imolecule&& mol) noexcept;
		local_rd_cnd(const Imolecule& mol_neut,const Imolecule& mol_cation,const Imolecule& mol_anion);
		local_rd_cnd(const local_rd_cnd& lrd_rhs);
		local_rd_cnd& operator=(const local_rd_cnd& lrd_rhs);
		local_rd_cnd(local_rd_cnd&& lrd_rhs) noexcept;
		local_rd_cnd& operator=(local_rd_cnd&& lrd_rhs) noexcept;
		~local_rd_cnd(){};
		//member functions
		friend local_rd_cnd operator-(const local_rd_cnd& lrd_lhs,const local_rd_cnd& lrd_rhs);
		void calculate_Fukui();
		void calculate_Fukui_band(int band);
		void calculate_Fukui_EW(int band);
		void calculate_RD(const global_rd& grd);
		void calculate_Hardness(const global_rd& grd);
		protein_lrd rd_protein(const Iprotein& prot);
		void write_rd_protein_pdb(const Iprotein& protein);
		void write_rd_protein_reaction(const Iprotein& prot);
		void write_LRD();
};

#endif 
//================================================================================
//END OF FILE
//================================================================================