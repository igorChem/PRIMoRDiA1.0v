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

#include <string>
#include <vector>
#include "../include/ReactionAnalysis.h"
#include "../include/reaction_coord.h"


using std::vector;
using std::string;

/********************************************************************/
ReactionAnalysis::ReactionAnalysis():
	ndim(1)							,
	nrcs(1)							,
	dimX(0)							,
	dimY(0)							{	
}
/********************************************************************/
ReactionAnalysis::~ReactionAnalysis(){}
/********************************************************************/
ReactionAnalysis::ReactionAnalysis(const ReactionAnalysis& rhs)	:
	RCs(rhs.RCs)												,
	mnt_atoms(rhs.mnt_atoms)									,
	mnt_residues(rhs.mnt_residues)								,
	rc1_indxs(rhs.rc1_indxs)									,
	rc2_indxs(rhs.rc2_indxs)									,
	ndim(rhs.ndim)												,
	nrcs(rhs.nrcs)												,	
	dimX(rhs.dimX)												,
	dimY(rhs.dimY)												{	
}
/********************************************************************/
ReactionAnalysis& ReactionAnalysis::operator=(const ReactionAnalysis& rhs){
	if ( this != &rhs ){
		RCs			= rhs.RCs;
		mnt_atoms	= rhs.mnt_atoms;
		mnt_residues= rhs.mnt_residues;
		rc1_indxs	= rhs.rc1_indxs;
		rc2_indxs	= rhs.rc2_indxs;
		ndim		= rhs.ndim;
		nrcs		= rhs.nrcs;
		dimX		= rhs.dimX;
		dimY		= rhs.dimY;
	}
	return *this;
}
/********************************************************************/
ReactionAnalysis::ReactionAnalysis(ReactionAnalysis&& rhs) noexcept	:
	RCs( move(rhs.RCs) )											,
	mnt_atoms( move(rhs.mnt_atoms) )								,
	mnt_residues( move(rhs.mnt_residues) )							,
	rc1_indxs( move(rhs.rc1_indxs) )								,
	rc2_indxs( move(rhs.rc2_indxs) )								,
	ndim(rhs.ndim)													,
	nrcs(rhs.nrcs)													,
	dimX(rhs.dimX)													,
	dimY(rhs.dimY)													{	
}
/********************************************************************/
ReactionAnalysis& ReactionAnalysis::operator=(ReactionAnalysis&& rhs) noexcept{
	if ( this != &rhs ){
		RCs			= move(rhs.RCs);
		mnt_atoms	= move(rhs.mnt_atoms);
		mnt_residues= move(rhs.mnt_residues);
		rc1_indxs	= move(rhs.rc1_indxs);
		rc2_indxs	= move(rhs.rc2_indxs);
		ndim		= rhs.ndim;
		nrcs		= rhs.nrcs;
		dimX		= rhs.dimX;
		dimY		= rhs.dimY;
	}
	return *this;
}
/********************************************************************/
void ReactionAnalysis::update(){
	nrcs = RCs.size();
	if ( rc2_indxs.size() > 2 ){
		ndim = 2;
	}
}
/////////////////////////////////////////////////////////////////////