//QMdriver.h

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

#ifndef QMDRIVER
#define QMDRIVER
//C++ headers
#include <iostream>
#include <memory>
//PRIMoRDiA headers
#include "../include/common.h"
//-------------------------
//foward declarations 
class Iaorbital;
class Imolecule;
//======================================================================================
/**
 * @class QMdriver
 * @author Igor Barden Grillo
 * @date 21/09/18
 * @file QMdriver.h
 * @brief This class is meant to manipulate quantum chemical coefficients from output programs 
 * to calculate properties for each atoms
 */
class QMdriver{
	public:	
		//member varibles
		Imolecule molecule;
		bool data_ok;
		//constructors/destructor
		QMdriver();
		QMdriver(Imolecule&& mol) noexcept;
		QMdriver(const QMdriver& qmd_rhs) = delete;
		QMdriver& operator=(const QMdriver& qmd_rhs) = delete;
		~QMdriver();
		//member functions
		double MO_in_atoms(int atom, int MO,bool beta);
		double EAS_in_atoms(int atom,int band);
		double EAS_in_atoms_EW(int atom);
		double NAS_in_atoms(int atom,int band);
		double NAS_in_atoms_EW(int atom);	
		double fukushima(int atom,int band);
		double density_in_atoms(int atom);
};

#endif 
//================================================================================
//END OF FILE
//================================================================================