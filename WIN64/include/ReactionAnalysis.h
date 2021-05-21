//ReactionAnalysis.h

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

#ifndef REACANALYSIS
#define REACANALYSIS

using std::vector;
using std::string;

class RC;

class ReactionAnalysis{
	public:
		std::vector<RC> RCs;
		std::vector<int> mnt_atoms;
		std::vector<int> mnt_residues;
		std::vector<int> rc1_indxs; 
		std::vector<int> rc2_indxs; 
		unsigned int ndim;
		unsigned int nrcs;
		unsigned int dimX;
		unsigned int dimY;
		ReactionAnalysis();
		~ReactionAnalysis();
		ReactionAnalysis(const ReactionAnalysis& rhs);
		ReactionAnalysis& operator=(const ReactionAnalysis& rhs);
		ReactionAnalysis(ReactionAnalysis&& rhs) noexcept;
		ReactionAnalysis& operator=(ReactionAnalysis&& rhs) noexcept;
		void update();
};

#endif