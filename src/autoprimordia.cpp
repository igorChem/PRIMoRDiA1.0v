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
#include <algorithm>
#include <omp.h>
//=================================
//PRIMoRDiA headers
#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/global_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/gridgen.h"
#include "../include/pos_traj.h"
#include "../include/pair_rd.h"
#include "../include/reaction_coord.h"
#include "../include/primordia.h"
#include "../include/autoprimordia.h"
#include "../include/scripts.h"
#include "../include/ReactionAnalysis.h"
//===============================
// std functions alias
using std::move;
using std::string;
using std::stoi;
using std::cout;
using std::endl;
using std::vector;
using std::to_string;
/*************************************************************/
AutoPrimordia::AutoPrimordia(){}
/*************************************************************/
AutoPrimordia::AutoPrimordia(const char* file_list):
	m_file_list(file_list)						 {
	
	m_log->input_message("Starting the descriptors calculation!\n");
	m_log->input_message("Starting to process the input file parameters for reactivity descriptors!\n\n");
	
	Ibuffer list_f(file_list,true);
	for (unsigned i=0;i<list_f.nLines;i++ ){
		for (unsigned j=0;j<list_f.lines[i].words.size();j++ ){
			if ( list_f.lines[i].words[0] == "#RT" ){
				run_type = list_f.lines[i].words[1];
			}
			else if ( list_f.lines[i].words[0] == "#PR" ){
				if	( list_f.lines[i].words[j] == "eband" ){
					energy_crit = list_f.lines[i].get_int(j+1);
				}
				else if ( list_f.lines[i].words[j] == "dos" )		dos			= true;
				else if ( list_f.lines[i].words[j] == "extrard" )	extra_RD	= true;
				else if	( list_f.lines[i].words[j] == "Rscript" )	M_R			= true;
				else if	( list_f.lines[i].words[j] == "composite" )	comp_H		= true;
				else if	( list_f.lines[i].words[j] == "pymols" )	
					pymol_script = true;
			}
		}
	}
}
/*************************************************************/
void AutoPrimordia::init(){
	if		( run_type == "normal" ) {
		this->calculate_rd();
	}
	else if	( run_type == "trajectory" ){
		this->calculate_rd_from_traj();
		this->md_trajectory_analysis();
	}	
	else if	( run_type == "reaction" ){
		this->calculate_rd_from_traj();
		this->reaction_analysis();
	}
	this->write_global();
	
}
/*************************************************************/
void AutoPrimordia::calculate_rd(){
	
	//temporary variables
	unsigned int i		= 0;
	unsigned int mode	= 0;
	const char* neut	= ".";
	const char* cation	= ".";
	const char* anion 	= ".";
	string program		= ".";
	string locHard 		= ".";
	string btm			= "BD";
	int gridsize		= 0;
	int bgap			= 0;
	bool mep			= false;
	int charge			= 0;
	double dens_tmp		= 0.0;
	
	//------------------------------------------------------
	//opening the file again
	Ibuffer list_f (m_file_list,true);
	
	for ( i=1;i<list_f.nLines;i++ ){
		if ( list_f.lines[i].words.size() <= 0 ){
			m_log->input_message("There are no contents in the line! Verify your input file!\n");
			break;
		}
		else if ( list_f.lines[i].words[0][0] == '#' ){ continue; }	
		else{
			mode = list_f.lines[i].get_int(0);
			for ( unsigned j=0; j<list_f.lines[i].words.size(); j++ ){
				if		( list_f.lines[i].words[j]  == "mep") mep = true;
				else if ( list_f.lines[i].words[j]  == "vm" ) dens_tmp = list_f.lines[i].get_double(j+1);
				else if ( list_f.lines[i].words[j]	== "EW" ) btm = "EW";
				else if ( list_f.lines[i].words[j]	== "BD" ) btm = "BD";
			}
			m_log->inp_delim(2);
			m_log->input_message("Starting New Entry!\n");
			primordia rd;			
			switch ( mode ){
				case 1:
					neut	= list_f.lines[i].words[1].c_str();
					locHard	= list_f.lines[i].words[2];
					gridsize= list_f.lines[i].get_int(3);
					program	= list_f.lines[i].words[4];
					rd.init_FOA(neut,gridsize,locHard,mep,program,dens_tmp);
				break;
				case 2:
					neut	= list_f.lines[i].words[1].c_str();
					cation	= list_f.lines[i].words[2].c_str();
					anion	= list_f.lines[i].words[3].c_str();
					locHard	= list_f.lines[i].words[4];
					gridsize= list_f.lines[i].get_int(5);
					charge	= list_f.lines[i].get_int(6);
					program	= list_f.lines[i].words[7];
					rd.init_FD(neut,cation,anion,gridsize,charge,mep,locHard,program,dens_tmp);
				break;
				case 3:
					neut	= list_f.lines[i].words[1].c_str();
					locHard	= list_f.lines[i].words[2];
					gridsize= list_f.lines[i].get_int(3);
					bgap	= list_f.lines[i].get_int(4);
					cation	= list_f.lines[i].words[5].c_str();
					program = list_f.lines[i].words[6];
					double r_atom[3];
					r_atom[0] = list_f.lines[i].get_double(7);
					r_atom[1] = list_f.lines[i].get_double(8);
					r_atom[2] = list_f.lines[i].get_double(9);
					int sze   = list_f.lines[i].get_int(10);
					rd.init_protein_RD(neut,locHard,gridsize,bgap,r_atom,sze,cation,mep,btm,program);
				break;
			}
			RDs.emplace_back( move(rd) ); 
		}
	}	
}
/*************************************************************/
void AutoPrimordia::calculate_rd_from_traj(){
	
	unsigned int mode = 3;
	unsigned int np_ngrd = 1;
	string temp_name	= "";
	string temp_name2	= "";
	string temp_name3	= "";
	string prefix		= "";
	string out_ext		= "";
	string program		= "";
	string locHard 		= ".";
	string btm			= "BD";
	int gridsize		= 0;
	int bgap			= 0;
	bool mep			= false;
	int charge			= 0;
	double dens_tmp		= 0.0;
	double r_atom[3];
	int sze				= 0;
	int start 			= 0;
	string pdb_prefix	= ".";
	vector<string> neut;
	vector<string> anions;
	vector<string> cations;
	vector<string> pdbs;
	
	
	Ibuffer list_f (m_file_list,true);
	
	for ( unsigned i=0; i<list_f.nLines; i++ ){
		for ( unsigned j=0; j<list_f.lines[i].words.size(); j++ ){
			if ( list_f.lines[i].words[0] == "#Reaction" ){
				if ( list_f.lines[i].words[j] == "dimX" ){
					trj_info.dimX = list_f.lines[i].get_int(j+1);
				}
				if ( list_f.lines[i].words[j] == "dimY" ){
					trj_info.dimY = list_f.lines[i].get_int(j+1);
				}
				if ( list_f.lines[i].words[j] == "start" ){
					start = list_f.lines[i].get_int(j+1);
				}
			}else if ( list_f.lines[i].words[0] == "#TRJ" ){
				if ( list_f.lines[i].words[j] == "frames" ){
					trj_info.dimX = list_f.lines[i].get_int(j+1);
				}
				if ( list_f.lines[i].words[j] == "start" ){
					start = list_f.lines[i].get_int(j+1);
				}
				if ( list_f.lines[i].words[j] == "residues"){
					for ( unsigned k=j+1; k<list_f.lines[i].words.size(); k++ ){
						trj_info.mnt_residues.push_back( list_f.lines[i].get_int(k) );
					}
				}
			}
			else if  ( list_f.lines[i].words[0][0]	== '#' ) { 
				continue;
			}
			else{
				m_log->inp_delim(1);
				mode	= list_f.lines[i].get_int(0);
				prefix	= remove_extension( list_f.lines[i].words[1].c_str() );
				out_ext = get_file_ext( list_f.lines[i].words[1].c_str() );
				locHard	= list_f.lines[i].words[2];
				gridsize= list_f.lines[i].get_int(3);
				switch( mode ){
					case 1:
						program = list_f.lines[i].words[4];
					break;
					case 2:
						charge	= list_f.lines[i].get_int(4);
						program = list_f.lines[i].words[5];
					break;	
					case 3:
						bgap		= list_f.lines[i].get_int(4);
						pdb_prefix	= remove_extension( list_f.lines[i].words[5].c_str() );
						program 	= list_f.lines[i].words[6];
						r_atom[0]	= list_f.lines[i].get_double(7);
						r_atom[1]	= list_f.lines[i].get_double(8);
						r_atom[2]	= list_f.lines[i].get_double(9);
						sze			= list_f.lines[i].get_int(10);
					break;
				}
				for( unsigned k=0; k<list_f.lines[i].words.size(); k++ ){
					if ( list_f.lines[i].words[k]	== "EW" ){ 
						btm = "EW";
					}
					if ( list_f.lines[i].words[k]	== "vm" ) { 
						dens_tmp = list_f.lines[k].get_double(j+1);
					}	
					if ( list_f.lines[i].words[k] == "mep"){ mep = true; }
				}
			}
		}
	}
	m_log->input_message("Calculating Reactivity Descriptors for a Reaction Path Trajectory!\n");
	
	if ( trj_info.dimY == 0 ){
		trj_info.dimY = 1;
	}
	if ( trj_info.dimY > 1 ){
		trj_info.ndim = 2;
	}				
				
	trj_info.rc1_indxs.resize( trj_info.dimX*trj_info.dimY );
	trj_info.rc2_indxs.resize( trj_info.dimY*trj_info.dimY );
	
	RDs.resize( trj_info.rc1_indxs.size() );
	if ( gridsize == 0 ){
		np_ngrd = NP;
		NP = 1;
	}

	unsigned cnt = 0;
	for( unsigned x=0; x<trj_info.dimX; x++ ){
		for( unsigned y=0; y<trj_info.dimY; y++ ){
			trj_info.rc1_indxs[cnt] = x +start;
			if ( trj_info.ndim == 2 ){
				trj_info.rc2_indxs[cnt] = y+start;
			}
			cnt++;
		}
	}
	
	temp_name = prefix;
	temp_name2= prefix;
	temp_name3= prefix;
	for( unsigned i=0; i<RDs.size(); i++ ){
		temp_name += to_string( trj_info.rc1_indxs[i] );
		if ( trj_info.ndim == 2 ){
			temp_name += "_";
			temp_name += to_string( trj_info.rc2_indxs[i] );
		}
		temp_name += out_ext;
		neut.push_back(temp_name);
		if ( mode == 2 ){
			temp_name2 +="_cat";
			temp_name2 += to_string( trj_info.rc1_indxs[i]  );
			temp_name3 +="_an";
			temp_name3 += to_string( trj_info.rc1_indxs[i]  );
			if ( trj_info.ndim == 2 ){
				temp_name2 += "_";
				temp_name2 += to_string( trj_info.rc2_indxs[i] );
				temp_name3 += "_";
				temp_name3 += to_string( trj_info.rc2_indxs[i] );
				temp_name2= prefix;
				temp_name3= prefix;
			}
			temp_name2 += out_ext;
			cations.push_back(temp_name2);
			temp_name3 += out_ext;
			anions.push_back(temp_name3);
		}
		else if ( mode == 3 ){
			temp_name2 = pdb_prefix;
			temp_name2 += to_string( trj_info.rc1_indxs[i] );
			if ( trj_info.ndim == 2 ){
				temp_name2 += "_";
				temp_name2 += to_string( trj_info.rc2_indxs[i] );
			}
			temp_name2 += ".pdb";
			pdbs.push_back(temp_name2);
			temp_name2= pdb_prefix;
		}
		temp_name = prefix;
	}
	
	unsigned i;
	switch ( mode ){
		case 1:
			omp_set_num_threads(np_ngrd);
			#pragma omp parallel
			{
				#pragma omp for
				for( i=0; i<RDs.size(); i++ ){
					RDs[i].init_FOA(neut[i].c_str(),gridsize,locHard,mep,program,dens_tmp);
				}
			}
		break;
		case 2:
			omp_set_num_threads(np_ngrd);
			#pragma omp parallel
			{
				#pragma omp for 
				for( i=0; i<RDs.size(); i++ ){
					RDs[i].init_FD(neut[i].c_str(),cations[i].c_str(),anions[i].c_str(),gridsize,charge,mep,locHard,program,dens_tmp);
				}
			}
		break;
		case 3:
			omp_set_num_threads(np_ngrd);
			#pragma omp parallel
			{
				#pragma omp for 
				for( i=0; i<RDs.size(); i++ ){
					RDs[i].init_protein_RD(neut[i].c_str(),locHard,gridsize,bgap,r_atom,sze,pdbs[i].c_str(),mep,btm,program);
				}
			}
		break;
	}
}
/*************************************************************/
void AutoPrimordia::reaction_analysis(){

	//temporary vars
	unsigned int pr = 0;

	vector<int> at_1; // reaction coordinate index 
	vector<int> at_2; // reaction coordinate index 
	vector< vector<int> > pr_ind;
	//--------------------------------------------------------
	
	m_log->input_message("Starting to analyse the Descriptors for the Reaction Path Trajectory!\n");
	
	Ibuffer list_f (m_file_list,true);
	
	for ( unsigned i=0; i<list_f.nLines; i++ ){		
		for ( unsigned j=0; j<list_f.lines[i].words.size(); j++ ){
			if ( list_f.lines[i].words[0] == "#Reaction" ){
				if ( list_f.lines[i].words[j] == "RC1" ){
					for( unsigned l=2; l<list_f.lines[i].words.size(); l++ ){
						at_1.push_back( list_f.lines[i].get_int(l) );
					}
					trj_info.RCs.emplace_back( RDs, at_1 );
					trj_info.mnt_atoms = at_1;
				}else if ( list_f.lines[i].words[j] == "RC2" ){
					for( unsigned l=2; l<list_f.lines[i].words.size(); l++){
						at_2.push_back( list_f.lines[i].get_int(l) );
						trj_info.mnt_atoms.push_back( at_2[ at_2.size()-1 ] );
					}
					trj_info.RCs.emplace_back( RDs, at_2 );
				}else if ( list_f.lines[i].words[j] == "Pair" ){
					pr_ind.resize(pr+1);
					pr_ind[pr].push_back( list_f.lines[i].get_int(j+1) );
					pr_ind[pr++].push_back( list_f.lines[i].get_int(j+2) );
				}
			}
		}
	}
	
	std::sort( trj_info.mnt_atoms.begin(), trj_info.mnt_atoms.end() );
	trj_info.mnt_atoms.erase( std::unique(trj_info.mnt_atoms.begin(), trj_info.mnt_atoms.end() ), trj_info.mnt_atoms.end() );
	
	//--------------------------------------------------------
	traj_rd atoms_lrd( RDs, trj_info.mnt_atoms, trj_info.mnt_residues );
	vector<pair_rd> pairs_lrds;
	for( unsigned i=0; i<pr ; i++ ){
		pairs_lrds.emplace_back( RDs, pr_ind[i][0], pr_ind[i][1] );
	}
	//------------------------------------------------------------------------
	string file_name = RDs[0].mol_info.name;
	file_name += ".atom_lrd"; 
	std::ofstream file_lrd( file_name.c_str() );
	file_lrd << std::fixed;
	file_lrd.precision(6);
	
	file_lrd << "RC1 "; 
	if ( trj_info.RCs.size() == 2 ){
		file_lrd << "RC2 ";
		trj_info.nrcs = 2;
	} 
	file_lrd << "rc1 ";
	if ( trj_info.ndim == 2 ){
		file_lrd << "rc2 ";
	}
	
	for( unsigned i=0; i<atoms_lrd.rds_labels.size() ; i++ ){
		file_lrd << atoms_lrd.rds_labels[i] << " ";
	}
	
	for( unsigned i=0; i<pr; i++ ){
		for ( unsigned j=0; j<pairs_lrds[i].rd_labels.size(); j++){
			file_lrd << pairs_lrds[i].rd_labels[j] << " ";
		}
	}
		
	file_lrd << "Energy HOF ECP Hardness Softness Electrophilicity";
	
	file_lrd << endl;
	
	for ( unsigned i=0; i<RDs.size(); i++){
		for( unsigned j=0; j<trj_info.RCs.size(); j++){
			file_lrd << trj_info.RCs[j].crd[i]	<< " ";
		}
		
		file_lrd << trj_info.rc1_indxs[i] << " ";
		
		if ( trj_info.ndim == 2 ){
			file_lrd << trj_info.rc2_indxs[i] << " ";
		}
		
		for ( unsigned k=0; k<trj_info.mnt_atoms.size(); k++ ){
			for ( unsigned m=0; m<9; m++ ){
				file_lrd << atoms_lrd.atoms_rd[k][m][i] << " ";
			}
		}
		for ( unsigned n=0; n<pr; n++ ){
			file_lrd << pairs_lrds[n].CTP[i]	<< " ";
			file_lrd << pairs_lrds[n].HPI_A[i]	<< " ";
			file_lrd << pairs_lrds[n].HPI_B[i]	<< " ";
			file_lrd << pairs_lrds[n].HPI_C[i]	<< " ";
			file_lrd << pairs_lrds[n].SPI[i]	<< " ";
			file_lrd << pairs_lrds[n].EEP[i]	<< " ";
		}
		file_lrd	<< RDs[i].grd.grds[2]	- RDs[0].grd.grds[2]		<< " "
					<< RDs[i].grd.grds[13]	- RDs[0].grd.grds[13]		<< " "
					<< RDs[i].grd.grds[7]	- RDs[0].grd.grds[7]		<< " "
					<< RDs[i].grd.grds[8]	- RDs[0].grd.grds[8]		<< " "
					<< RDs[i].grd.grds[9]	- RDs[0].grd.grds[9]		<< " ";
		file_lrd << endl;
	}
	file_lrd.close();
	
	vector<string> pr_lbs;
	for( unsigned i=0; i<pr; i++){
		for( unsigned j=0; j<pairs_lrds[i].rd_labels.size(); j++){
			pr_lbs.push_back(pairs_lrds[i].rd_labels[j]);
		}
	}
		
	if ( M_R ){
		scripts r_analysis( RDs[0].mol_info.name.c_str(), "reaction_analsys" );
		r_analysis.write_r_reaction_analysis(atoms_lrd,pr_lbs,trj_info,file_name);
		if ( trj_info.mnt_residues.size() > 0 ){
			scripts r_residues_analysis( RDs[0].mol_info.name.c_str(), "residues_analysis" );
			r_residues_analysis.write_r_residuos_barplot();
		}
	}
	
	atoms_lrd.calculate_res_stats();
	atoms_lrd.write_residues_reports();
}
/*************************************************************/
void AutoPrimordia::md_trajectory_analysis(){
	traj_rd trajectories( RDs, trj_info.mnt_atoms, trj_info.mnt_residues );
	trajectories.calculate_res_stats();
	trajectories.write_residues_reports();
	if ( M_R ){
		scripts r_residues_analysis( RDs[0].mol_info.name.c_str(), "residues_analysis" );
		r_residues_analysis.write_r_residuos_barplot();
	}
}
/*************************************************************/
void AutoPrimordia::write_global(){

	
	string fn = change_extension( m_file_list, ".global");
	std::ofstream file_grd(fn.c_str());
	
	file_grd << "GRD ";
	for(unsigned j = 0; j<RDs[0].grd.grds.size(); j++ ){
		file_grd << RDs[0].grd.rd_abrev[j] << " ";
	}
	
	file_grd << endl;
	file_grd << std::fixed;
	file_grd.precision(8);
	for(unsigned i = 0; i<RDs.size(); i++ ){
		file_grd << RDs[i].grd.name				<< " ";
			for(unsigned j = 0; j<RDs[i].grd.grds.size(); j++ ){
				file_grd << RDs[i].grd.grds[j] << " ";
			}
			file_grd << endl;
	} 
	file_grd.close();
}
/*************************************************************/
AutoPrimordia::~AutoPrimordia(){
	m_log->input_message("Ending reactivity descriptors calculations.\n");
}
//================================================================================
//END OF FILE
//================================================================================
