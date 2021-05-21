// source file for global_rd.cpp source file 
// class for global reactivity descriptors 

//-------------------------------------------------------------------
/* Class to provide objects to store global reactivity descriptor 
 * information and methods to calculate them, also to write and print
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
//including c++ header
#include <iostream>
#include <string> 
#include <fstream>
#include <cmath>
#include <iomanip>

//including PRIMoRDiA header
#include "../include/common.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/global_rd.h"

using std::move;

//================================================================================================
//Class member functions definitions
/*****************************************************************************************/
global_rd::global_rd()		:
	name("nonamed")			,
	homo_en(0.0)			,
	lumo_en(0.0)			,
	energ_tot(0.0)			,
	energ_cat(0.0)			,
	energ_an(0.0)			,
	hardness(0.0)			,
	chemical_pot(0.0)		,
	softness(0.0)			,
	IP(0.0)					,
	EA(0.0)					,
	Electrophilicity(0.0)	,
	gap(0.0)				,
	deltaN(0.0)				,
	nMax(0.0)				,
	HOF(0.0)				,
	KA(false)				,
	DF(false)				{
}
/*****************************************************************************************/
global_rd::global_rd(const Imolecule& mol)		:
	name(mol.name)								,
	homo_en(mol.homo_energy)					,
	lumo_en(mol.lumo_energy)					,
	energ_tot(mol.energy_tot)					,
	energ_cat(0.0)								,
	energ_an(0.0)								,
	hardness(0.0)								,
	chemical_pot(0.0)							,
	softness(0.0)								,
	IP(-mol.homo_energy)						,
	EA(-mol.lumo_energy)						,
	Electrophilicity(0.0)						,
	gap( mol.homo_energy -mol.lumo_energy )		,
	deltaN(0.0)									,
	nMax(0.0)									,
	HOF(mol.heat_of_formation)					,
	KA(true)									,
	DF(false)									{
	
}
/*****************************************************************************************/
global_rd::global_rd(const Imolecule& mol_neutro			,
					 const Imolecule& mol_cation			,
					 const Imolecule& mol_anion)			:
	name(mol_neutro.name)									,
	homo_en(mol_neutro.homo_energy)							,
	lumo_en(mol_neutro.lumo_energy)							,
	energ_tot(mol_neutro.energy_tot)						,
	energ_cat(mol_cation.energy_tot)						,
	energ_an(mol_anion.energy_tot)							,
	hardness(0.0)											,
	chemical_pot(0.0)										,
	softness(0.0)											,
	IP(energ_cat - energ_tot)								,
	EA(energ_tot - energ_an)								,
	Electrophilicity(0.0)									,
	gap( mol_neutro.homo_energy-mol_neutro.lumo_energy )	,
	deltaN(0.0)												,
	nMax(0.0)												,
	HOF(mol_neutro.heat_of_formation)						,
	KA(false)												,
	DF(true)												{
	
}
/*****************************************************************************************/
global_rd::global_rd(const global_rd& rd_rhs)		:
		name(rd_rhs.name)							,
		homo_en(rd_rhs.homo_en)						,
		lumo_en(rd_rhs.lumo_en)						,
		energ_tot(rd_rhs.energ_tot)					,
		energ_an(rd_rhs.energ_an)					,
		energ_cat(rd_rhs.energ_cat)					,
		hardness(rd_rhs.hardness)					,
		softness(rd_rhs.softness)					,
		chemical_pot(rd_rhs.chemical_pot)			,
		IP(rd_rhs.IP)								,
		EA(rd_rhs.EA)								,
		Electrophilicity(rd_rhs.Electrophilicity)	,
		gap(rd_rhs.gap)								,
		deltaN(rd_rhs.deltaN)						,
		DF(rd_rhs.DF)								,
		KA(rd_rhs.KA)								,
		nMax(rd_rhs.nMax)							,
		HOF(rd_rhs.HOF)								{

}
/*****************************************************************************************/
global_rd& global_rd::operator=(const global_rd& rd_rhs){
	if (this!=&rd_rhs){
		name			= rd_rhs.name;
		homo_en			= rd_rhs.homo_en;
		lumo_en			= rd_rhs.lumo_en;
		energ_tot		= rd_rhs.energ_tot;
		energ_an		= rd_rhs.energ_an;
		energ_cat		= rd_rhs.energ_cat;
		hardness		= rd_rhs.hardness;
		softness		= rd_rhs.softness;
		chemical_pot	= rd_rhs.chemical_pot;
		IP				= rd_rhs.IP;
		EA				= rd_rhs.EA;
		Electrophilicity= rd_rhs.Electrophilicity;
		gap				= rd_rhs.gap;
		deltaN			= rd_rhs.deltaN;
		nMax			= rd_rhs.nMax;
		HOF				= rd_rhs.HOF;
		KA				= rd_rhs.KA;
		DF				= rd_rhs.DF;
	}
	return *this;
}
/*****************************************************************************************/
global_rd::global_rd(global_rd&& rd_rhs)  noexcept :
		name( std::move(rd_rhs.name) )				,
		homo_en(rd_rhs.homo_en)						,
		lumo_en(rd_rhs.lumo_en)						,
		energ_tot(rd_rhs.energ_tot)					,
		energ_an(rd_rhs.energ_an)					,
		energ_cat(rd_rhs.energ_cat)					,
		hardness(rd_rhs.hardness)					,
		softness(rd_rhs.softness)					,
		chemical_pot(rd_rhs.chemical_pot)			,
		IP(rd_rhs.IP)								,
		EA(rd_rhs.EA)								,
		Electrophilicity(rd_rhs.Electrophilicity)	,
		gap(rd_rhs.gap)								,
		deltaN(rd_rhs.deltaN)						,
		DF(rd_rhs.DF)								,
		KA(rd_rhs.KA)								,
		nMax(rd_rhs.nMax)							,
		HOF(rd_rhs.HOF)								{
}
/*****************************************************************************************/
global_rd& global_rd::operator=(global_rd&& rd_rhs) noexcept {
	if (this!=&rd_rhs){
		name			= std::move (rd_rhs.name);
		homo_en			= rd_rhs.homo_en;
		lumo_en			= rd_rhs.lumo_en;
		energ_tot		= rd_rhs.energ_tot;
		energ_an		= rd_rhs.energ_an;
		energ_cat		= rd_rhs.energ_cat;
		softness		= rd_rhs.softness;
		hardness		= rd_rhs.hardness;
		chemical_pot	= rd_rhs.chemical_pot;
		IP				= rd_rhs.IP;
		EA				= rd_rhs.EA;
		Electrophilicity= rd_rhs.Electrophilicity;
		gap				= rd_rhs.gap;
		deltaN			= rd_rhs.deltaN;
		nMax			= rd_rhs.nMax;
		KA				= rd_rhs.KA;
		DF				= rd_rhs.DF;
		HOF				= rd_rhs.HOF;
	}	
	return *this;
}
/*****************************************************************************************/
void global_rd::calculate_rd(){
		chemical_pot= -(IP+EA)/2;
		hardness= (IP-EA)/2;
		softness= 1/hardness;
		Electrophilicity = chemical_pot*chemical_pot*softness/2;
		nMax= -chemical_pot*softness*2;
		if ( gap < 0 ){
			gap*=-1;
		}
}
/*****************************************************************************************/
void global_rd::print_rd(){
	std::string typestr;
	std::string typestr2;
	if ( KA ) { 
		typestr  = "Koopman approximation \n";
		typestr2 = "KA";
	}
	if ( DF ) {
		typestr  = "Finite differences approximation\n";
		typestr2 = "FD";
	}	
	std::cout << "Molecule name "						<< name
			  << "\n Calculus method "					<< typestr
			  << "\n Printing reactivity descriptors"	<< std::endl;
	std::cout << "Ionization potential " 				<< IP << "\n" 
			  << "Electron affinity: "					<< EA << "\n"
			  << "Chemical potential: "					<< chemical_pot << "\n"
			  << "Hardness: "							<< hardness << "\n"
			  << "Softness: "							<< softness << "\n"
			  << "Electrophilicity: "					<< Electrophilicity << "\n"
			  << "Max electron recible: "				<< nMax << "\n"
			  << "HOMO-LUMO gap: "						<< gap << std::endl; 
}
/*****************************************************************************************/
global_rd operator-(const global_rd& lhs_grd, const global_rd& rhs_grd){
	global_rd Result(lhs_grd);	
	Result.hardness			= lhs_grd.hardness			- rhs_grd.hardness;
	Result.softness			= lhs_grd.softness			- rhs_grd.softness;
	Result.chemical_pot		= lhs_grd.chemical_pot		- rhs_grd.chemical_pot;
	Result.IP				= lhs_grd.IP				- rhs_grd.IP;
	Result.EA				= lhs_grd.EA				- rhs_grd.EA;
	Result.Electrophilicity	= lhs_grd.Electrophilicity	- rhs_grd.Electrophilicity;
	Result.gap				= lhs_grd.gap				- rhs_grd.gap;
	Result.nMax				= lhs_grd.nMax				- rhs_grd.nMax;
	Result.deltaN			= (lhs_grd.chemical_pot		- rhs_grd.chemical_pot) / (lhs_grd.hardness + rhs_grd.hardness);
	return Result;
} 
/*****************************************************************************************/
void global_rd::write_rd(){
	std::string typestr;
	std::string typestr2;
	if ( KA ) { 
		typestr  = "Koopman approximation \n";
		typestr2 = "KA";
	}
	if ( DF ) {
		typestr  = "Finite differences approximation\n";
		typestr2 = "FD";
	}	
	std::string file_name = name + typestr2 + ".GRD";
	std::ofstream file_grd;
	file_grd.open(file_name.c_str());
	
	file_grd << std::fixed;
	file_grd.precision(5);

	file_grd  << "Calculus method: "	<< typestr
			  << "name: "				<<  name << "\n"
			  << "HOF "					<< std::setw(10) << std::left
			  << "Energy "				<< std::setw(10) << std::right
			  << "ECation "				<< std::setw(10) << std::right
			  << "EAnion "				<< std::setw(10) << std::right
			  << "IP "					<< std::setw(10) << std::right
			  << "EA "					<< std::setw(10) << std::right
			  << "ECP "					<< std::setw(10) << std::right
			  << "Hardness "			<< std::setw(10) << std::right
			  << "Softness "			<< std::setw(10) << std::right
			  << "Electrophilicity "	<< std::setw(10) << std::right
			  << "nMax "				<< std::setw(10) << std::right
			  << "gap"					<< std::endl;
	
	file_grd << HOF				<< " " << std::setw(10) << std::left
			 << energ_tot		<< " " << std::setw(10) << std::right
			 << energ_cat		<< " " << std::setw(10) << std::right
			 << energ_an		<< " " << std::setw(10) << std::right
			 << IP				<< " " << std::setw(10) << std::right
			 << EA				<< " " << std::setw(10) << std::right
			 << chemical_pot	<< " " << std::setw(10) << std::right
			 << hardness		<< " " << std::setw(10) << std::right
			 << softness		<< " " << std::setw(10) << std::right
			 << Electrophilicity<< " " << std::setw(10) << std::right
			 << nMax<< " "		<< std::setw(10) << std::right
			 << gap				<< std::endl;
	file_grd.close();
}
/*****************************************************************************************/
global_rd::~global_rd(){}
//================================================================================
//END OF FILE
//================================================================================

