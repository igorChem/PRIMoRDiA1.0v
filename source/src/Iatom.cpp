// Iatom.cpp 
//--------------------------------------------------------------------------
/** source codes for classes that will represent molecules and hold chemical and geometrical 
information 
*/

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

//including c++ headerrs
#include <iostream> 
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
//--------------------------------------------------------------------------
//including PRIMoRDiA headers
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"

using std::string;
using std::move;

//=========================================================================================
// Classes member functions definitions
//------------------------------------------------------------------------------------------
/*******************************************************************************************/
Iatom::Iatom()			: 
	xcoord(0)			, 
	ycoord(0)			,
	zcoord(0)			,
	element("H")		,
	charge(0)			,
	atomic_mass(1.00794),
	atomicN(1)			,
	norb(0)				{
}
/*******************************************************************************************/
Iatom::Iatom(double x	, 
					double y	,
					double z	,
					string typ)	:
	xcoord(x)					,
	ycoord(y)					,
	zcoord(z)					, 
	element( move(typ) )		,
	charge(0)					,
	atomic_mass(0)				,
	atomicN(0)					,
	norb(0)						{
		
	atomicN 	= get_atomic_number(element); 
	atomic_mass = get_atom_mass(element);
}
/*******************************************************************************************/
Iatom::Iatom(const Iatom& rhs_atom)		: 
	xcoord(rhs_atom.xcoord)				,
	ycoord(rhs_atom.ycoord)				,
	zcoord(rhs_atom.zcoord)				,
	element(rhs_atom.element)			,
	atomicN(rhs_atom.atomicN)			,
	atomic_mass(rhs_atom.atomic_mass)	,
	charge(rhs_atom.charge)				,
	orbitals(rhs_atom.orbitals)			,
	norb(rhs_atom.norb)					{ 
}
/*******************************************************************************************/
Iatom& Iatom::operator=(const Iatom& rhs_atom) {
	if(this != &rhs_atom){
		xcoord		= rhs_atom.xcoord; 
		ycoord		= rhs_atom.ycoord; 
		zcoord		= rhs_atom.zcoord; 
		element		= rhs_atom.element;
		charge		= rhs_atom.charge;
		atomicN		= rhs_atom.atomicN;
		atomic_mass	= rhs_atom.atomic_mass;
		norb		= rhs_atom.norb;
		orbitals	= rhs_atom.orbitals;
	}
	return *this;
}
/*******************************************************************************************/
Iatom::Iatom(Iatom&& rhs_atom) noexcept	:
	xcoord(rhs_atom.xcoord)				,
	ycoord(rhs_atom.ycoord)				,
	zcoord(rhs_atom.zcoord)				,
	element( move(rhs_atom.element) )	,
	atomicN(rhs_atom.atomicN)			,
	atomic_mass(rhs_atom.atomic_mass)	,
	charge(rhs_atom.charge)				,
	orbitals( move(rhs_atom.orbitals) )	,
	norb(rhs_atom.norb)					{
}
/*******************************************************************************************/
Iatom& Iatom::operator=(Iatom&& rhs_atom) noexcept {
	if( this != &rhs_atom ){
		xcoord		= rhs_atom.xcoord; 
		ycoord		= rhs_atom.ycoord; 
		zcoord		= rhs_atom.zcoord; 
		element		= move(rhs_atom.element);
		charge		= rhs_atom.charge;
		atomicN		= rhs_atom.atomicN;
		atomic_mass = rhs_atom.atomic_mass;
		norb		= rhs_atom.norb;
		orbitals	= move(rhs_atom.orbitals);
	}
	return *this;
}
/*******************************************************************************************/
void Iatom::print() {
	std::cout << "x coordinate: "		<< xcoord		<< "\n" 
		 << "y coordinate: "			<< ycoord		<< "\n"
		 << "x coordinate: "			<< zcoord		<< "\n"
		 << "The atom is a "			<< element		<< " element " << "of atomic number " << atomicN << "\n"
		 << "The partial charge is: "	<< charge		<< "\n"
		 << "The atomic mass is: "		<< atomic_mass	<< "\n"
		 << "The atomic number is: "	<< atomicN		<< "\n"
		 << "Number of orbitals :"		<< norb			<< "\n"
		 << std::endl;
	std::cout << "Printing the atomic orbitals " << std::endl;
	for (int i=0;i<norb;i++) orbitals[i].print();
}
/*******************************************************************************************/
void Iatom::set_type(string typ) {
	element= move(typ);
	atomicN= get_atomic_number(element);
	atomic_mass = get_atom_mass(element);
}
/*******************************************************************************************/
void Iatom::set_coord(double x, double y, double z) {
	xcoord = x;
	ycoord = y;
	zcoord = z;
}
/*******************************************************************************************/
void Iatom::add_orbital(int level, string sym, double coef){
	orbitals.emplace_back(level,sym,coef);		
	norb++;
}
/*******************************************************************************************/
void Iatom::add_orbital(Iaorbital orb){
	orbitals.emplace_back( move(orb) );
	norb++;
}
/*******************************************************************************************/
Iatom::~Iatom(){}
//================================================================================
//END OF FILE
//================================================================================