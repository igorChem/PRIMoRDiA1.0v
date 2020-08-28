//autoprimordia.cpp

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

//=================================
//C++ Headers 
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <memory>

//=================================
//PRIMoRDiA headers
#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/Imolecule.h"
#include "../include/global_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/Icube.h"
#include "../include/gridgen.h"
#include "../include/QMparser.h"
#include "../include/primordia.h"
#include "../include/autoprimordia.h"
//===============================
// std functions alias
using std::unique_ptr;
using std::move;
using std::string;
using std::stoi;
using std::cout;
using std::endl;
/*************************************************************/
AutoPrimordia::AutoPrimordia(){
}
/*************************************************************/
AutoPrimordia::AutoPrimordia(const char* file_list)	:
	m_file_list(file_list)												{
	
	bool compare		= false;
	unsigned int i		= 0;
	unsigned int mode= 0;
	
	const char* neut	= ".";
	const char* cation	= ".";
	const char* anion 	= ".";
	
	const char* homo_cube	= "."; 
	const char* lumo_cube	= ".";
	const char* neut_cube	= ".";
	const char* cub_cation	= ".";
	const char* cub_anion	= ".";
	
	string program	= ".";
	string locHard 	= ".";
	string btm			= "BD";
	int gridsize		= 0;
	int bgap			= 0;
	bool LH				= false; 
	bool mep			= false;
	bool ed				= false;
	int charge			= 0;
	
	//----------------------------------------------------------------------
	m_log->input_message("Starting the descriptors calculation!\n");
	//----------------------------------------------------------------------
	unique_ptr<Ibuffer> list_f ( new Ibuffer(file_list,true) );
	//----------------------------------------------------------------------
	m_log->input_message("Starting to read input information given in the first input line.\n");
	//----------------------------------------------------------------------
	for (i=0;i<list_f->lines[0].words.size();i++ ){
		if			( list_f->lines[0].words[i]	== "compare") compare		= true; 
		else if	( list_f->lines[0].words[i]	== "dos"		) dos			= true;
		else if	( list_f->lines[0].words[i]	== "pymols"		) pymol_script	= true;
		else if	( list_f->lines[0].words[i]	== "extrard"	) extra_RD		= true;
		else if	( list_f->lines[0].words[i]	== "eband"		) energy_crit	= std::stod(list_f->lines[0].words[i+1]);
		else if	( list_f->lines[0].words[i]	== "Rscript"	) M_R			= true;
		else if	( list_f->lines[0].words[i]	== "traj"		) {
			atoms.resize( std::stoi(list_f->lines[0].words[i+1] ) );
			for(int j=0;j<atoms.size();j++){
				atoms[j] = std::stoi(list_f->lines[0].words[i+j+2] ) ;
			}
		}
	}
	//----------------------------------------------------------------------
	for (i=1;i<list_f->nLines;i++){
		if ( list_f->lines[i].words.size() <= 0 ){
			break;
		}else{
			 mode = stoi(list_f->lines[i].words[0]);
			 if ( mode == 1 ) {
				m_log->input_message("");
				neut	= list_f->lines[i].words[1].c_str();
				locHard	= list_f->lines[i].words[2];
				gridsize= stoi(list_f->lines[i].words[3]);
				program	= list_f->lines[i].words[4];
				for (int j=0;j<list_f->lines[i].words.size();j++ ){
					if	( list_f->lines[i].words[j]  == "mep" ) mep = true;
				}
				primordia rd;
				rd.init_FOA(neut,gridsize,locHard,mep,program);
				RDs.emplace_back( move(rd) ); 
			 }else if ( mode == 2 ){
				m_log->input_message("");
				neut	= list_f->lines[i].words[1].c_str();
				cation	= list_f->lines[i].words[2].c_str();
				anion	= list_f->lines[i].words[3].c_str();
				locHard	= list_f->lines[i].words[4];
				gridsize= stoi(list_f->lines[i].words[5]);
				charge	= stoi(list_f->lines[i].words[6]);
				program	= list_f->lines[i].words[7];
				for (int j=0;j<list_f->lines[i].words.size();j++ ){
					if  ( list_f->lines[i].words[j]  == "mep" ) mep = true;
				}
				primordia rd;
				rd.init_FD(neut,cation,anion,gridsize,charge,mep,locHard,program);
				RDs.emplace_back( move(rd) ); 
			}else if ( mode == 3 ){
				m_log->input_message("");
				neut	= list_f->lines[i].words[1].c_str();
				locHard	= list_f->lines[i].words[2];
				gridsize= stoi(list_f->lines[i].words[3]);
				bgap	= stoi(list_f->lines[i].words[4]);
				cation	= list_f->lines[i].words[5].c_str();
				string program = list_f->lines[i].words[6];
				double r_atom[3];
				r_atom[0] = stod(list_f->lines[i].words[7]);
				r_atom[1] = stod(list_f->lines[i].words[8]);
				r_atom[2] = stod(list_f->lines[i].words[9]);
				int sze   = stoi(list_f->lines[i].words[10]);
				for (int j=0;j<list_f->lines[i].words.size();j++ ){
					if  ( list_f->lines[i].words[j]		== "mep" )	mep	= true;
					else if ( list_f->lines[i].words[j]	== "EW" )	btm	= "EW";
				}
				primordia rd;
				rd.init_protein_RD(neut,locHard,gridsize,bgap,r_atom,sze,cation,mep,btm,program);
				RDs.emplace_back( move(rd) );
			}else{
				continue;
			}
		}
	}
	m_log->input_message("Ending reactivity descriptors calculations.");
	/*
	if ( compare ){
		for (i = 1;i<RDs.size();i++){
			std::unique_ptr<compPrimordia> RD_comp ( new compPrimordia(RDs[0],RDs[i]) );
			RD_comp->write_report();
		}
	}
	 */
}
/*************************************************************/
void AutoPrimordia::write_global(){
	
	string file_name = m_file_list;
	file_name 		+= ".global"; 
	
	std::ofstream file_grd( file_name.c_str() );
	file_grd << "Global RD for each system in the PRIMoRDiA input\n";
	file_grd << "Molecule Energy Ecat Ean HOF IP EA ECP Hardness Softness Electrophilicity nMax gap\n";
	file_grd << std::fixed;
	file_grd.precision(8);
	for(unsigned i = 0;i<RDs.size();i++){
		file_grd << RDs[i].grd->name				<< " "
					<< RDs[i].grd->energ_tot		<< " "
					<< RDs[i].grd->energ_cat		<< " "
					<< RDs[i].grd->energ_an 		<< " "
					<< RDs[i].grd->HOF				<< " "
					<< RDs[i].grd->IP				<< " "
					<< RDs[i].grd->EA				<< " "
					<< RDs[i].grd->chemical_pot		<< " "
					<< RDs[i].grd->hardness			<< " "
					<< RDs[i].grd->softness			<< " "
					<< RDs[i].grd->Electrophilicity	<< " "
					<< RDs[i].grd->nMax				<< " "
					<< RDs[i].grd->gap
					<< endl;
	} 
	file_grd.close();
}
/*************************************************************/
void AutoPrimordia::traj_atoms(){
	
	unsigned int i,j,k,l;
	
	if ( atoms.size() > 0 ){
		traj_atoms_rd.resize( atoms.size() );
		for( i=0;i<atoms.size();i++ ){
			traj_atoms_rd[i].resize(4);
			for( j=0;j<4;j++ ){
				traj_atoms_rd[i][j].resize( RDs.size() );
			}
		}
	
		for (i=0;i<traj_atoms_rd.size();i++){
			for( j=0;j<RDs.size();j++){
					traj_atoms_rd[i][0][j] = RDs[j].lrdCnd->EAS[atoms[i]-1];
					traj_atoms_rd[i][1][j] = RDs[j].lrdCnd->NAS[atoms[i]-1];
					traj_atoms_rd[i][2][j] = RDs[j].lrdCnd->hardness_B[atoms[i]-1];
					traj_atoms_rd[i][3][j] = RDs[j].lrdCnd->hardness_C[atoms[i]-1];
			}
		}
	
		string file_name = m_file_list;
		file_name 			+= ".atom_lrd"; 
	
		std::ofstream file_lrd( file_name.c_str() );
		file_lrd << std::fixed;
		file_lrd.precision(6);
	
		for (j=0;j<RDs.size();j++){
			for (i=0;i<traj_atoms_rd.size();i++){
				file_lrd<< traj_atoms_rd[i][0][j] << " " 
						<< traj_atoms_rd[i][1][j] << " " 
						<< traj_atoms_rd[i][2][j] << " "
						<< traj_atoms_rd[i][3][j] << " ";
			}
			file_lrd << endl;
		}
		file_lrd.close();
	}
}
/*************************************************************/
AutoPrimordia::~AutoPrimordia(){}
/*************************************************************/
//================================================================================
//END OF FILE
//================================================================================