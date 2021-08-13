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
	KA(false)				,
	DF(false)				{
		
	rd_names = {"HOMO_ENERGY",				// 0
				"LUMO_ENERGY",				// 1
				"TOTAL_ENERGY",				// 2
				"ENERGY_CATION",			// 3
				"ENERGY_ANION",				// 4
				"IONIZATION_POTENTIAL",		// 5
				"ELECTRONIC_AFFINITY",		// 6
				"CHEMICAL_POTENTIAL",		// 7
				"HARDNESS",					// 8
				"SOFTNESS",					// 9
				"TOTAL_ELECTROPHILICITY",	// 10
				"GAP",						// 11
				"N_MAX",					// 12
				"HEAT_OF_FORMATION"			// 13
	};
	
	rd_abrev = {"HOMO_E", "LUMO_E", "T_ENERGY", "ENERGY_CAT","ENERGY_AN", "IP", "EA", "ECP", 
				"HARDNESS","SOFTNESS","Electrophilicity","GAP", "N_MAX", "HOF"};
		
	grds.resize( rd_names.size() );
}
/*****************************************************************************************/
global_rd::global_rd(const Imolecule& mol)		:
	name(mol.name)								,
	KA(true)									,
	DF(false)									{
	
	rd_names = {"HOMO_ENERGY",				// 0
				"LUMO_ENERGY",				// 1
				"TOTAL_ENERGY",				// 2
				"ENERGY_CATION",			// 3
				"ENERGY_ANION",				// 4
				"IONIZATION_POTENTIAL",		// 5
				"ELECTRONIC_AFFINITY",		// 6
				"CHEMICAL_POTENTIAL",		// 7
				"HARDNESS",					// 8
				"SOFTNESS",					// 9
				"TOTAL_ELECTROPHILICITY",	// 10
				"GAP",						// 11
				"N_MAX",					// 12
				"HEAT_OF_FORMATION"			// 13
	};
	
	rd_abrev = {"HOMO_E", "LUMO_E", "T_ENERGY", "ENERGY_CAT","ENERGY_AN", "IP", "EA", "ECP", 
				"HARDNESS","SOFTNESS","Electrophilicity","GAP", "N_MAX", "HOF"};
	
	grds.resize( rd_names.size() );
	
	grds[0] = mol.homo_energy;
	grds[1] = mol.lumo_energy;
	grds[2] = mol.energy_tot;
	grds[5] = -grds[0];
	grds[6] = -grds[1];
	grds[11]= grds[0]-grds[1];
	grds[14]= mol.heat_of_formation;
}
/*****************************************************************************************/
global_rd::global_rd(const Imolecule& mol_neutro			,
					 const Imolecule& mol_cation			,
					 const Imolecule& mol_anion)			:
	name(mol_neutro.name)									,
	KA(false)												,
	DF(true)												{
		
	rd_names = {"HOMO_ENERGY",				// 0
				"LUMO_ENERGY",				// 1
				"TOTAL_ENERGY",				// 2
				"ENERGY_CATION",			// 3
				"ENERGY_ANION",				// 4
				"IONIZATION_POTENTIAL",		// 5
				"ELECTRONIC_AFFINITY",		// 6
				"CHEMICAL_POTENTIAL",		// 7
				"HARDNESS",					// 8
				"SOFTNESS",					// 9
				"TOTAL_ELECTROPHILICITY",	// 10
				"GAP",						// 11
				"N_MAX",					// 12
				"HEAT_OF_FORMATION"			// 13
	};
	
	rd_abrev = {"HOMO_E", "LUMO_E", "T_ENERGY", "ENERGY_CAT","ENERGY_AN", "IP", "EA", "ECP", 
				"HARDNESS","SOFTNESS","Electrophilicity","GAP", "N_MAX", "HOF"};
	
	grds.resize( rd_names.size() );
	
	grds[0] = mol_neutro.homo_energy;
	grds[1] = mol_neutro.lumo_energy;
	grds[2] = mol_neutro.energy_tot;
	grds[5] = mol_neutro.energy_tot - mol_cation.energy_tot;
	grds[6] = mol_anion.energy_tot - mol_neutro.energy_tot;
	grds[11]= grds[0]-grds[1];
	grds[13]= mol_neutro.heat_of_formation;
	
}
/*****************************************************************************************/
global_rd::global_rd(const global_rd& rd_rhs)		:
		name(rd_rhs.name)							,
		rd_names(rd_rhs.rd_names)					,
		rd_abrev(rd_rhs.rd_abrev)					,
		grds(rd_rhs.grds)							,
		DF(rd_rhs.DF)								,
		KA(rd_rhs.KA)								{

}
/*****************************************************************************************/
global_rd& global_rd::operator=(const global_rd& rd_rhs){
	if (this!=&rd_rhs){
		name			= rd_rhs.name;
		rd_names		= rd_rhs.rd_names;
		rd_abrev		= rd_rhs.rd_abrev;
		grds			= rd_rhs.grds;
		KA				= rd_rhs.KA;
		DF				= rd_rhs.DF;
	}
	return *this;
}
/*****************************************************************************************/
global_rd::global_rd(global_rd&& rd_rhs)  noexcept :
		name( std::move(rd_rhs.name) )				,
		rd_names( std::move(rd_rhs.rd_names) )		,
		rd_abrev( std::move(rd_rhs.rd_abrev) )		,
		grds( std::move(rd_rhs.grds) )				,
		DF(rd_rhs.DF)								,
		KA(rd_rhs.KA)								{
}
/*****************************************************************************************/
global_rd& global_rd::operator=(global_rd&& rd_rhs) noexcept {
	if (this!=&rd_rhs){
		name			= std::move(rd_rhs.name);
		rd_names		= std::move(rd_rhs.rd_names);
		rd_abrev		= std::move(rd_rhs.rd_abrev);
		grds			= std::move(rd_rhs.grds);
		KA				= rd_rhs.KA;
		DF				= rd_rhs.DF;
	}	
	return *this;
}
/*****************************************************************************************/
void global_rd::calculate_rd(){
		grds[7] = -( grds[5] + grds[6] )/2; // electronic chemical potential
		grds[8] = ( grds[5] - grds[6] )/2; //hardness
		grds[9] = 1/grds[8]; // softness
		
		grds[10] = grds[7]*grds[7]*grds[9]/2; // total electrophilicity
		grds[12] = -grds[7]*grds[9]*2; // nmax of electrons that the system would receive
		if ( grds[11] < 0 ){
			grds[11] *=-1;
		}
}
/*****************************************************************************************/
void global_rd::print_rd(){
	std::string typestr;
	std::string typestr2;
	if ( KA ) { 
		typestr  = "Frozen Orbital Approximation \n";
		typestr2 = "FOA";
	}
	if ( DF ) {
		typestr  = "Finite Differences approximation\n";
		typestr2 = "FD";
	}	
	std::cout << "Molecule name "						<< name
			  << "\n Calculus method "					<< typestr
			  << "\n Printing reactivity descriptors"	<< std::endl;
	std::cout << "Ionization potential " 				<< grds[5] << "\n" 
			  << "Electron affinity: "					<< grds[6] << "\n"
			  << "Chemical potential: "					<< grds[7] << "\n"
			  << "Hardness: "							<< grds[8] << "\n"
			  << "Softness: "							<< grds[9] << "\n"
			  << "Electrophilicity: "					<< grds[10] << "\n"
			  << "Max electron recible: "				<< grds[12] << "\n"
			  << "HOMO-LUMO gap: "						<< grds[11] << std::endl; 
}
/*****************************************************************************************/
global_rd operator-(const global_rd& lhs_grd, const global_rd& rhs_grd){
	global_rd Result(lhs_grd);	
	for (unsigned i=0; i<lhs_grd.grds.size() ;i++ ){
		Result.grds[i] = lhs_grd.grds[i] - rhs_grd.grds[i];
	}
	return Result;
} 
/*****************************************************************************************/
void global_rd::write_rd(){
	std::string typestr;
	std::string typestr2;
	if ( KA ) { 
		typestr  = "Frozen Ortbial Approximation \n";
		typestr2 = "FOA";
	}
	if ( DF ) {
		typestr  = "Finite differences approximation\n";
		typestr2 = "FD";
	}	
	std::string file_name = name + typestr2 + ".GRD";
	std::ofstream file_grd;
	file_grd.open( file_name.c_str() );
	
	file_grd << std::fixed;
	file_grd.precision(5);

	file_grd  << "Calculus method: "	<< typestr
			  << "name: "				<<  name << "\n";
			  
	file_grd <<  std::endl;
	
	for( unsigned int i=0; i<rd_abrev.size(); i++){
		file_grd << rd_abrev[i] << std::setw(10) << std::left;
	}	
	file_grd <<  std::endl;
	
	for( unsigned int i=0; i<grds.size(); i++){
		file_grd << grds[i] << std::setw(10) << std::left;
	}	
	file_grd.close();
}
/*****************************************************************************************/
global_rd::~global_rd(){}
//================================================================================
//END OF FILE
//================================================================================

