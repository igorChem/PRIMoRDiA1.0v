//residue.h

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
// include statements from c++ library

#ifndef RESIDUELRD
#define RESIDUELRD

#include <vector>
#include <string> 


class Iprotein;
//=============================================================================
/**
 * @class residue_lrd
 * @author Igor Barden Grillo
 * @date 05/07/20
 * @file residue_lrd.h
 * @brief Class to store reactivity descripors for a biological polymer residue, and its stats.
 */
class residue_lrd{
	public:
		//member variables
		std::vector<double> rd_sum;
		std::vector<double> rd_avg;
		//constructors/destructor
		residue_lrd();
		residue_lrd(const residue_lrd& rhs);
		residue_lrd(residue_lrd&& rhs) noexcept;
		residue_lrd& operator=(const residue_lrd& rhs);
		residue_lrd& operator=(residue_lrd&& rhs) noexcept;
		~residue_lrd();
		//member functions
};

//=============================================================================
/**
 * @class protein_lrd
 * @author Igor Barden Grillo
 * @date 05/07/20
 * @file residue_lrd.h
 * @brief Class to store objects with reactivity descriptors for all residues in a biological polymer structure.
 */
class protein_lrd{
	public:
		//member variables
		std::vector<residue_lrd> residues_rd;
		std::vector<double>	protein_sum_avg;
		std::vector<double>	protein_avg_avg;
		std::vector<double>	protein_max;
		std::vector<double>	protein_min;
		std::vector<std::string> labels;
		std::vector<int> hydrophobicity;
		//constructors/destructor
		protein_lrd();
		protein_lrd(std::vector<residue_lrd> res_rd);
		protein_lrd(const protein_lrd& rhs);
		protein_lrd(protein_lrd&& rhs) noexcept;
		protein_lrd& operator=(const protein_lrd& rhs);
		protein_lrd& operator=(protein_lrd&& rhs) noexcept;
		~protein_lrd();
		//member function
		void write_protein_lrd(const Iprotein& prot);
		void determine_hydrophilicity();
		void recalculate_stats();
};

#endif

//================================================================================
//END OF FILE
//================================================================================