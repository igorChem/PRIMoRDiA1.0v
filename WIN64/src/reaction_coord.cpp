//reaction_coord.cpp

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

#include "../include/reaction_coord.h"
#include "../include/primordia.h"
#include "../include/Imolecule.h"
#include "../include/Iatom.h"

using std::vector;

/*************************************************************/
RC::RC()			:
	nsteps(0)		,
	d_indx(0)		,
	n_atoms(0)		,
	rc_label("rc0")	{
}
/*************************************************************/
RC::~RC(){}
/*************************************************************/
RC::RC( std::vector<primordia>& trj_info,
		vector<int>& atoms)				:
		nsteps(0)						,
		d_indx(0)						,
		n_atoms( atoms.size() )			,
		rc_label("RC")					,
		atom_indx(atoms)				{
	
			
	double tmp_double  = 0.000;
	for( int i =0;i<trj_info.size();i++){
		if ( n_atoms == 2 ){
			crd.push_back( trj_info[i].mol_info.calc_dist( atom_indx[0]-1, atom_indx[1]-1 ) );
		}else if ( n_atoms == 3 ){
			tmp_double = trj_info[i].mol_info.calc_dist( atom_indx[0]-1, atom_indx[1]-1 );
			crd.push_back( tmp_double - trj_info[i].mol_info.calc_dist( atom_indx[1]-1, atom_indx[2]-1 ) );
		}
	}
	tmp_double = mean_dvec(crd);
	for( int i=0; i<crd.size(); i++){
		crd[i] -= tmp_double;
	} 
	
	rc_label += " (";
	rc_label += trj_info[0].mol_info.atoms[ atom_indx[0]-1 ].element;
	rc_label += std::to_string(atom_indx[0]);
	rc_label += "--";
	rc_label += trj_info[0].mol_info.atoms[ atom_indx[1]-1 ].element;
	rc_label += std::to_string(atom_indx[1]);
	if ( atom_indx.size() == 3 ){
		rc_label += "--";
		rc_label += trj_info[0].mol_info.atoms[ atom_indx[2]-1 ].element;
		rc_label += std::to_string(atom_indx[2]);
	}
	rc_label += ")";

}
/*************************************************************/
RC::RC(const RC& rhs)		:
	nsteps(rhs.nsteps)		,
	d_indx(rhs.d_indx)		,
	n_atoms(rhs.n_atoms)	,
	rc_label(rhs.rc_label)	,
	atom_indx(rhs.atom_indx),
	crd(rhs.crd)			{
}
/*************************************************************/
RC& RC::operator=(const RC& rhs){
	if ( this != &rhs ){
		nsteps		= rhs.nsteps;
		d_indx		= rhs.d_indx;
		n_atoms		= rhs.n_atoms;
		rc_label	= rhs.rc_label;
		atom_indx	= rhs.atom_indx;
		crd			= rhs.crd;
	}
	return *this;
}
/*************************************************************/
RC::RC(RC&& rhs) noexcept			:
	nsteps(rhs.nsteps)				,
	d_indx(rhs.d_indx)				,
	n_atoms(rhs.n_atoms)			,
	rc_label( move(rhs.rc_label) )	,
	atom_indx( move(rhs.atom_indx) ),
	crd( move(rhs.crd) )			{
	
}
/*************************************************************/
RC& RC::operator=(RC&& rhs) noexcept{
	if ( this != &rhs ){
		nsteps		= rhs.nsteps;
		d_indx		= rhs.d_indx;
		n_atoms		= rhs.n_atoms;
		rc_label	= rhs.rc_label;
		atom_indx	= rhs.atom_indx;
		crd			= rhs.crd;
	}
	return *this;
}


///////////////////////////////////////////////////////////////