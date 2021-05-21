//reaction_coord.h

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
#include <string>

class primordia;
//-------------------------------------------------
class RC{
	public:
		unsigned int nsteps;
		unsigned int d_indx;
		unsigned int n_atoms;
		std::string rc_label;
		std::vector<int> atom_indx;
		std::vector<double> crd;
		RC();
		~RC();
		RC( std::vector<primordia>& trj_info, std::vector<int>& atoms);
		RC(const RC& rhs);
		RC& operator=(const RC& rhs);
		RC(RC&& rhs) noexcept;
		RC& operator=(RC&& rhs) noexcept;
};