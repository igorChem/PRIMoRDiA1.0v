//pair.h

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

#include <string>
#include <vector>

class primordia;
class pair_rd{
	public:
		std::vector<int> a_indices;
		std::string label;
		std::vector<std::string> rd_labels;
		std::vector<double> CTP;	// charge transfer pair
		std::vector<double> HPI_A;	// hard-hard pair interactions
		std::vector<double> HPI_B;	// hard-hard pair interactions
		std::vector<double> HPI_C;	// hard-hard pair interactions
		std::vector<double> HPI_D;	// hard-hard pair interactions
		std::vector<double> SPI;	// softness pair interactions
		std::vector<double> EEP;	// Estabilization energy interactions
		pair_rd();
		pair_rd(std::vector<primordia>& lrd, int ax, int ay);
		~pair_rd();
		pair_rd(const pair_rd& rhs);
		pair_rd& operator=(const pair_rd& rhs);
		pair_rd(pair_rd&& rhs) noexcept;
		pair_rd& operator=(pair_rd&& rhs) noexcept;
};