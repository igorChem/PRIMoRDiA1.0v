//comp_hardness.h

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
#include <vector>

class global_rd;
class local_rd_cnd;
class protein_lrd;
class Iprotein;

class comp_hard{
	public:
		double g_comp_hard;
		std::vector< std::vector<double> > l_comp_hard;
		std::vector< std::vector<double> > l_comp_hard_bio;
		comp_hard();
		comp_hard(const global_rd& grd, double vm);
		comp_hard(const global_rd& grd, const local_rd_cnd& lrd, double vm);
		~comp_hard();
		comp_hard(const comp_hard& rhs);
		comp_hard& operator=(const comp_hard& rhs);
		comp_hard( comp_hard&& rhs ) noexcept;
		comp_hard& operator=( comp_hard&& rhs ) noexcept;
		void calculate_protein(const protein_lrd& lrd, const Iprotein& prot);
		void write_comp_hardness(const char* name);
};