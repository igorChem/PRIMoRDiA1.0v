//Including c++ headers

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

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <algorithm>  

//Including PRIMoRDiA headers
//-------------------------------------------------------
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/mopac_files.h"
//-------------------------------------------------------
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>
//-------------------------------------------------------
// Aliases for standard c++ scope functions
using std::move;
using std::vector;
using std::string;
using std::unique_ptr;
using std::cout;
using std::endl;
using std::stoi;
using std::stod;

//======================================================================
mopac_files::mopac_files()		:
	RHF(true)					,
	LMO(false)					,
	f_chg(0)					,
	name_f("no_name")			,
	is_open(false)				,
	type("no_type")				,
	buffer( new Ibuffer() )		,
	molecule( new Imolecule() ){
}
/***************************************************************************************/
mopac_files::mopac_files(const char* file_name):
	molecule( new Imolecule() )					,
	buffer( new Ibuffer() )						,
	LMO(false)									,
	is_open(false)								,
	RHF(true)									{
	
	unsigned int i,j;
	name_f = file_name;
	molecule->name = get_file_name(file_name);
	molecule->name = remove_extension( molecule->name.c_str() );
	
	if ( IF_file( file_name ) ){
		is_open = true;
		if ( check_file_ext(".aux",file_name ) ) {
			type = "AUX";
			buffer.reset( new Ibuffer(file_name,1,15) );
			for(i=0;i<buffer->nLines;i++){
				if ( buffer->lines[i].IF_word("KEYWORDS",0,8) ){
					for(j=0;j<buffer->lines[i].line_len;j++){
						if ( buffer->lines[i].words[j] == "mozyme" || buffer->lines[i].words[j] == "MOZYME" ){
							LMO = true;
							RHF = true;
						}
					}
				}
			}
			
			//====================
			//reading keywords
			buffer.reset( new Ibuffer(name_f,"NUM_ELECTRONS","ATOM_X_OPT:ANGSTROMS[") );
			for(i=0;i<buffer->nLines;i++){
				if ( i ==0 ) {
					string tmp_a(buffer->lines[0].words[0],14, buffer->lines[0].words[0].size()-14);
					molecule->num_of_electrons = stoi( tmp_a);
					//cout << molecule->num_of_electrons << endl;
				}else{
					if ( buffer->lines[i].IF_word("HEAT_OF_FORMATION:",0,18)) {
						string tmp_ab(buffer->lines[i].words[0],27,14 );
						molecule->heat_of_formation = D_E_conv(tmp_ab);
					}else{
						if ( buffer->lines[i].IF_word("NUM_ALPHA_ELECTRONS=",0,20) ){
							string tmp_b(buffer->lines[i].words[0],20,buffer->lines[i].words[0].size()-20);
							if ( stoi(tmp_b) > 0 ) {
								RHF = false;
								molecule->num_of_electrons = 0;
								molecule->num_of_electrons += stoi(tmp_b);
								molecule->betad = true;
							}
						}else if ( buffer->lines[i].IF_word("NUM_BETA_ELECTRONS=",0,19) ){
							string tmp_c(buffer->lines[i].words[0],19,buffer->lines[i].words[0].size()-19);
							if ( stoi(tmp_c) > 0 ) {
								molecule->num_of_electrons += stoi(tmp_c);
							}
						}else if ( buffer->lines[i].IF_word("ENERGY_ELECTRONIC:EV=",0,21) ){
							string word1 ( buffer->lines[i].words[0],21,14);
							molecule->energy_tot = D_E_conv(word1);
						}
					}
				}
			}
		}else if ( check_file_ext(".out",file_name ) ) {
			type = "OUT";
			buffer.reset( new Ibuffer(name_f,14,25) );
			for(i=0;i<buffer->nLines;i++){
				if 	( buffer->lines[i].IF_word("MOZYME",1,6) 	){
					LMO	= true;
					RHF = true;
				}
				else if ( buffer->lines[i].IF_word("SINGLET",1,7) )	RHF	=true; 
				else if ( buffer->lines[i].IF_word("DOUBLET",1,7) )	RHF	=false; 
				else if ( buffer->lines[i].IF_word("TRIPLET",1,7) )	RHF =false; 
				else if ( buffer->lines[i].IF_word("QUARTET",1,7) )	RHF =false; 
				else if ( buffer->lines[i].IF_word("CHARGE",1,6) )	f_chg = stoi(buffer->lines[i].words[5]); 
			}
		}else if (check_file_ext(".mgf",file_name) ){
			type = "MGF";
		}else{
			cout << "Warning! The file has imcompatible extension name with mopac files!" << endl;
			m_log->input_message("Warning! The file has imcompatible extension name with mopac files!");
			is_open = false;
		}
	}else{
		m_log->input_message("Error opening MOPAC file! Verify its presence in the current directory!");
		cout << "the file named " << file_name << " cannot be opened! " << endl;
	}
}
/***************************************************************************************/
void mopac_files::parse_aux(){
	
	m_log->input_message("Starting to parse MOPAC AUX file.");
	unsigned i,j,k;
	vector<int> aoidx;
	vector<double> zetas;
	vector<int> shells;
	
	buffer.reset( new Ibuffer(name_f,"ATOM_EL[","ATOM_CORE[") );
	for(i=1; i<buffer->nLines-1;i++){
		for (j=0;j<buffer->lines[i].line_len;j++){
			Iatom atom;
			atom.set_type( buffer->lines[i].get_string(j) );
			molecule->add_atom(atom);
		}
	}
	//molecule->print();
	if ( molecule->atoms.size() == 0 ){ cout << "Warning! Zero atoms read in aux file. Verify your file!!" << endl; } 

	m_log->input_message("Found number of atoms in the aux file: ");
	m_log->input_message( int(molecule->num_of_atoms) );
	
	m_log->input_message("Number of electron in the system: ");
	m_log->input_message( int(molecule->num_of_electrons) );
	
	int counter = 0;
	buffer.reset ( new Ibuffer(name_f,"ATOM_X:ANGSTROMS[","AO_ATOMINDEX[") );
	for (i=1;i<buffer->nLines-1;i++ ){
		molecule->atoms[counter].xcoord 		= buffer->lines[i].pop_double(0);
		molecule->atoms[counter].ycoord 		= buffer->lines[i].pop_double(0);
		molecule->atoms[counter++].zcoord 	= buffer->lines[i].pop_double(0);
	}
	//molecule->print_coordinates();
	
	buffer.reset( new Ibuffer(name_f,"AO_ATOMINDEX[","ATOM_SYMTYPE[") );
	for (i=1;i<buffer->nLines-1;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) aoidx.push_back( buffer->lines[i].pop_int(0) );
	}
	
	m_log->input_message("Number of atomic orbitals: ");
	m_log->input_message( int( aoidx.size() ) );
	
	counter = 0;
	buffer.reset ( new Ibuffer(name_f,"ATOM_SYMTYPE[","AO_ZETA[") );
	for (i=1; i<buffer->nLines-1;i++ ){
		for (j=0; j<buffer->lines[i].line_len;j++){
			Iaorbital aorb;
			aorb.symmetry = buffer->lines[i].pop_string(0);
			if 		( aorb.symmetry == "PX" ) aorb.powx	= 1;
			else if ( aorb.symmetry == "PY" ) aorb.powy	= 1;
			else if ( aorb.symmetry == "PZ" ) aorb.powz = 1;
			molecule->atoms[aoidx[counter++]-1].add_orbital(aorb);
		}
	}

	buffer.reset ( new Ibuffer(name_f,"AO_ZETA[","ATOM_PQN[") );
	for (i=1;i<buffer->nLines-1;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) zetas.push_back( buffer->lines[i].pop_double(0) );		
	}
	
	buffer.reset ( new Ibuffer(name_f,"ATOM_PQN[","NUM_ELECTRONS=") );
	for (i=1;i<buffer->nLines-1;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) shells.push_back( buffer->lines[i].pop_int(0) );		
	}

	counter = 0;
	buffer.reset ( new Ibuffer(name_f,"ATOM_CHARGES[","OVERLAP_MATRIX[") );
	for (i=1;i<buffer->nLines-1;i++){
		for (j=0;j<buffer->lines[i].line_len;j++) molecule->atoms[counter++].charge = buffer->lines[i].pop_double(0);
	}
	
	molecule->get_ao_number();
	
	omp_set_num_threads(NP);
	#pragma omp parallel
	{
		#pragma omp single nowait 
		{
			this->get_overlap_m();
		}
		#pragma omp single nowait 
		{
			this->get_mo(false);
		}
		#pragma omp single nowait 
		{
			this->get_mo_energies(false);
		}
		if ( !RHF) {
			#pragma omp single nowait
			{
			this->get_mo(true);
			}
			#pragma omp single nowait
			{
			this->get_mo_energies(true);
			}
		}
	}
	
	buffer.reset(nullptr);
	
	counter = 0;
	for (i=0;i<molecule->atoms.size();i++){
		for (j=0;j<molecule->atoms[i].norb;j++){
			molecule->atoms[i].orbitals[j].alpha = zetas[counter];
			molecule->atoms[i].orbitals[j].shell = shells[counter++];
		}
	}
	molecule->update();
	if ( !molecule->check() ) { cout << "Problems in reading the mopac aux file: " << name_f << endl;}
	m_log->input_message("HOMO energy: ");
	m_log->input_message( double_to_string(molecule->homo_energy) );
	m_log->input_message("HOMO number: ");
	m_log->input_message( double_to_string(molecule->homoN) );
	m_log->input_message("LUMO energy: ");
	m_log->input_message( double_to_string(molecule->lumo_energy) );
	m_log->input_message("LUMO number: ");
	m_log->input_message( double_to_string(molecule->lumoN) );
}
/***************************************************************************************/
void mopac_files::parse_out(){
	m_log->input_message("Starting to parse out file from MOPAC.");
	
	buffer.reset( new Ibuffer(name_f,true) );
	for (int i=0;i<buffer->nLines;i++){
		if ( buffer->lines[i].IF_line("HEAT",1,"FORMATION",3,10) ){ 
			 molecule->heat_of_formation = buffer->lines[i].pop_double(5);
		}
		else if ( buffer->lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,5) || buffer->lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,8) ){
			molecule->energy_tot = buffer->lines[i].pop_double(3);
			//cout << molecule->energy_tot << endl;
		}
		else if ( buffer->lines[i].IF_line("HOMO",0,"LUMO",1,7) ){
			molecule->homo_energy = buffer->lines[i].pop_double(5);
			molecule->lumo_energy = buffer->lines[i].pop_double(5);
		}
		else if ( buffer->lines[i].IF_line("SUM",0,5) ){
			for ( int j=0;j<3;j++){ molecule->dipole_moment[j] = buffer->lines[i].pop_double(1); }
			molecule->total_dipmoment = buffer->lines[i].pop_double(1);
		}
	}

	m_log->input_message("Total Energy ");
	m_log->input_message(double_to_string(molecule->energy_tot));
	m_log->input_message("HOMO energy ");
	m_log->input_message(double_to_string(molecule->homo_energy));
	m_log->input_message("LUMO energy ");
	m_log->input_message(double_to_string(molecule->lumo_energy));
	
}
/***************************************************************************************/
void mopac_files::parse_mgf(){
	vector<double>	zetasS;
	vector<double>	zetasP;
	vector<double>	zetasD;
	vector<double>	inv_mat;
	vector<int>		orbN;
	vector<int>		orbN_beta;
	
	unsigned int i,j	= 0;
	int inmat_in		= 0;
	int inmat_fin		= 0;
	int noa				= 0;
	double tot_charge	= 0;
 	
	m_log->input_message("Starting to parse MOPAC MGF file!");
	
	
	string mgf_name = change_extension(name_f,".mgf");
	
	if ( IF_file( mgf_name.c_str() ) ){
	
		buffer.reset( new Ibuffer(mgf_name.c_str(),true) );
		int counter = 0;
		for(i=0;i<buffer->nLines;i++){
			if( i == 0 ) noa = buffer->lines[i].pop_int(0);
			else if ( i>0 && i<=noa ) {
				Iatom atom;
				atom.set_type( get_atomic_symbol( buffer->lines[i].pop_int(0) ) );
				atom.xcoord  = buffer->lines[i].pop_double(0);
				atom.ycoord  = buffer->lines[i].pop_double(0);
				atom.zcoord  = buffer->lines[i].pop_double(0);
				atom.charge  = buffer->lines[i].pop_double(0);
				molecule->add_atom(atom);
			}
			else if( i>noa && i<=(noa*2) ){
				zetasS.push_back( buffer->lines[i].pop_double(0) );
				zetasP.push_back( buffer->lines[i].pop_double(0) );
				zetasD.push_back( buffer->lines[i].pop_double(0) );
			}
			else if( buffer->lines[i].IF_word("ORBITAL",0,7) && inmat_fin == 0 )
				orbN.push_back(i);
			else if( buffer->lines[i].IF_word("INVERSE_MATRIX[",0,15) )
				inmat_in = i;
			else if( buffer->lines[i].IF_word("Keywords:",0,9) ) 
				inmat_fin = i;
			else if( inmat_fin > 0 && buffer->lines[i].IF_word("ORBITAL",0,7) ){
				molecule->betad = true;
				orbN_beta.push_back(i);
			}
		}

		string temp = "noname";
		for(int i=0;i<orbN.size();i++){
			int fin_ind = 0;
			if ( i==orbN.size()-1 )	fin_ind = inmat_in;
			else fin_ind = orbN[i+1];
			for(int j=orbN[i];j<fin_ind;j++){
				if ( j == orbN[i] ){
					molecule->occupied.push_back( buffer->lines[j].pop_int(1) );
					molecule->orb_energies.push_back( buffer->lines[j].pop_double(2) );
					molecule->MOnmb++;
				}else{ 
					for(int k=0;k<buffer->lines[j].line_len;k++){
						if (k == 0 &&  buffer->lines[j].words[k][0] == '-' ){
							if ( buffer->lines[j].words[k].size() > 15 ){
								temp = buffer->lines[j].words[k].substr(0,15);
								molecule->coeff_MO.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 15,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 30,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule->coeff_MO.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr( 45,buffer->lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule->coeff_MO.push_back( D_E_conv(temp) );
											temp = buffer->lines[j].words[k].substr( 60,buffer->lines[j].words[k].size() );
											molecule->coeff_MO.push_back( D_E_conv(temp) );
										}else molecule->coeff_MO.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO.push_back( D_E_conv( buffer->lines[j].words[k]) );
						}else{
							if ( buffer->lines[j].words[k].size() > 14 ){
								temp = buffer->lines[j].words[k].substr(0,14);
								molecule->coeff_MO.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 14,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 29,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule->coeff_MO.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr( 44,buffer->lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule->coeff_MO.push_back( D_E_conv(temp) );
											temp = buffer->lines[j].words[k].substr(59,buffer->lines[j].words[k].size() );
											molecule->coeff_MO.push_back( D_E_conv(temp) );
										}else molecule->coeff_MO.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO.push_back( D_E_conv(buffer->lines[j].words[k]) );
						}	
					}
				}
			}
		}
		for(int i=0;i<noa;i++){
			int sh = 0;
			Iaorbital orbS;
			orbS.alpha = zetasS[i];
			orbS.shell = 1;
			if 		( molecule->atoms[i].atomicN > 1  && molecule->atoms[i].atomicN < 11 ) sh		= orbS.shell = 2;
			else if ( molecule->atoms[i].atomicN > 10 && molecule->atoms[i].atomicN < 19 ) sh	= orbS.shell = 3;
			else if ( molecule->atoms[i].atomicN > 18 && molecule->atoms[i].atomicN < 37 ) sh	= orbS.shell = 4;
			molecule->atoms[i].add_orbital(orbS);
			if ( zetasP[i] > 1e-05 ){
				Iaorbital orbPx, orbPy, orbPz;
				orbPx.alpha = orbPy.alpha = orbPz.alpha = zetasP[i];
				orbPx.shell = orbPy.shell = orbPz.shell = sh;
				orbPx.symmetry = "PX";
				orbPy.symmetry = "PY";
				orbPz.symmetry = "PZ";
				orbPx.powx  = 1;  
				orbPy.powy  = 1;  
				orbPz.powz  = 1;
				molecule->atoms[i].add_orbital(orbPx);
				molecule->atoms[i].add_orbital(orbPy);
				molecule->atoms[i].add_orbital(orbPz);
			}
			if ( zetasD[i] > 1e-05 ){
				Iaorbital orbDxx, orbDyy, orbDzz, orbDxy, orbDyz, orbDxz;
				orbDxx.alpha = orbDyy.alpha = orbDzz.alpha = orbDxy.alpha = orbDyz.alpha = orbDxz.alpha = zetasD[i];
				orbDxx.shell = orbDyy.shell = orbDzz.shell = orbDxy.shell = orbDyz.shell = orbDxz.shell = sh;
				orbDxx.powx  = 2;
				orbDxx.symmetry  = "XX";
				orbDyy.powy  = 2;
				orbDyy.symmetry  = "YY";
				orbDzz.powz  = 2;
				orbDzz.symmetry  = "zz";
				orbDxy.powx  = orbDxy.powy  = 1;
				orbDxy.symmetry  = "XY";
				orbDyz.powy  = orbDyz.powz  = 1;
				orbDyz.symmetry  = "YZ";
				orbDxz.powx  = orbDxz.powz  = 1;
				orbDxz.symmetry  = "XZ";
				molecule->atoms[i].add_orbital(orbDxx);
				molecule->atoms[i].add_orbital(orbDyy);
				molecule->atoms[i].add_orbital(orbDzz);
				molecule->atoms[i].add_orbital(orbDxy);
				molecule->atoms[i].add_orbital(orbDyz);
				molecule->atoms[i].add_orbital(orbDxz);
			}
		}
	
		for(int i=inmat_in+1;i<inmat_fin;i++){
			for(int k=0;k<buffer->lines[i].line_len;k++){
				temp = buffer->lines[i].words[k];
				if ( k == 0 && temp[0] == '-' ){
					if ( temp.size() > 15 ){ 
						temp = temp.substr(0,15);
						inv_mat.push_back( D_E_conv(temp) );
						temp = buffer->lines[i].words[k].substr( 15,buffer->lines[i].words[k].size() );
						if ( temp.size() > 15 ){
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = buffer->lines[i].words[k].substr( 30,buffer->lines[i].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = buffer->lines[i].words[k].substr( 45,buffer->lines[i].words[k].size() );
								if ( temp.size() > 15 ){
									temp = temp.substr(0,15);
									inv_mat.push_back( D_E_conv(temp) );
									temp = buffer->lines[i].words[k].substr( 60,buffer->lines[i].words[k].size() );
									inv_mat.push_back( D_E_conv(temp) );
								}else inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) );
				}else{
					if( temp.size() > 14 ){
						temp = temp.substr(0,14);
						inv_mat.push_back( D_E_conv(temp) );
						temp = buffer->lines[i].words[k].substr( 14,buffer->lines[i].words[k].size() );
						if ( temp.size() > 15 ){
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = buffer->lines[i].words[k].substr( 29,buffer->lines[i].words[k].size() );
							if ( temp.size() > 15 ) {
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = buffer->lines[i].words[k].substr( 44,buffer->lines[i].words[k].size() );
								if ( temp.size() > 15 ){
									temp = temp.substr(0,15);
									inv_mat.push_back( D_E_conv(temp) );
									temp = buffer->lines[i].words[k].substr( 59,buffer->lines[i].words[k].size() );
									inv_mat.push_back( D_E_conv(temp) );
								}else inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) ); 
				}
			}
		}	
		int nmo = molecule->MOnmb;
		int k   = 0;
	
	
		for(int i=0;i<orbN_beta.size();i++){
			int fin_ind = 0;
			if ( i==orbN_beta.size()-1 ) fin_ind = buffer->nLines;
			else fin_ind = orbN_beta[i+1];
			for(int j=orbN_beta[i];j<fin_ind;j++){
				if ( j == orbN_beta[i] ){
					molecule->occupied_beta.push_back( buffer->lines[j].pop_int(1) );
					molecule->orb_energies_beta.push_back( buffer->lines[j].pop_double(2) );
					molecule->MOnmb_beta++;
				}else{
					for(int k=0;k<buffer->lines[j].line_len;k++){
						if (k == 0 &&  buffer->lines[j].words[k][0] == '-' ){
							if ( buffer->lines[j].words[k].size() > 15 ){
								temp = buffer->lines[j].words[k].substr(0,15);
								molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 15,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 30,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr( 45,buffer->lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
											temp = buffer->lines[j].words[k].substr(60,buffer->lines[j].words[k].size() );
											molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO_beta.push_back( D_E_conv(buffer->lines[j].words[k]) );
						}else{
							if ( buffer->lines[j].words[k].size() > 14 ){
								temp = buffer->lines[j].words[k].substr(0,14);
								molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = buffer->lines[j].words[k].substr( 14,buffer->lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = buffer->lines[j].words[k].substr( 29,buffer->lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = buffer->lines[j].words[k].substr( 44,buffer->lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
											temp = buffer->lines[j].words[k].substr(59,buffer->lines[j].words[k].size() );
											molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
										}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule->coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule->coeff_MO_beta.push_back( D_E_conv(buffer->lines[j].words[k]) );
						}		
					}
				}
			}
		}
		buffer.reset();
	
		Eigen::MatrixXd inv(nmo,nmo);
		Eigen::MatrixXd coeff_corrected(nmo,nmo);

		int bnmo = molecule->MOnmb_beta;
 
		Eigen::MatrixXd coeff_corrected_beta(bnmo,bnmo);
	
		for (i=0;i<nmo;i++){
			for (j=0;j<nmo;j++){  coeff_corrected(i,j) = molecule->coeff_MO[i*nmo+j]; }
		}
	
		for (i=0;i<bnmo;i++){
			for (j=0;j<bnmo;j++) coeff_corrected_beta(i,j) = molecule->coeff_MO_beta[i*nmo+j];
		}
	
		for(i=0;i<nmo;i++){
			for(j=0;j<=i;j++){
				if (i==j)  inv(i,j) = inv_mat[k++];
				else       inv(j,i) = inv(i,j) = inv_mat[k++];
			}
		}
	
		coeff_corrected = coeff_corrected*inv;
		for (int i=0;i<nmo;i++){
			for (int j=0;j<nmo;j++)  molecule->coeff_MO[i*nmo+j] = coeff_corrected(i,j);
		}
		
		if ( molecule->betad ){
			coeff_corrected_beta = coeff_corrected_beta*inv; 
			for (i=0;i<nmo;i++){
				for (j=0;j<nmo;j++)	molecule->coeff_MO_beta[i*nmo+j] = coeff_corrected_beta(i,j);
			}
		}
		
		for(int i=0;i<molecule->occupied.size();i++){
			molecule->num_of_electrons+= molecule->occupied[i];
		}
		
		molecule->update();
		
		m_log->input_message("HOMO energy: ");
		m_log->input_message(double_to_string(molecule->homo_energy));
		m_log->input_message("LUMO energy: ");
		m_log->input_message(double_to_string(molecule->lumo_energy));
		m_log->input_message("Atomic orbitals: ");
		m_log->input_message(molecule->get_ao_number() );
	
	}else{
		m_log->input_message("MGF file does not exist, skipping its parsing process!");
	}
}
/***************************************************************************************/
void mopac_files::get_overlap_m(){
	
	string keyword = "OVERLAP_MATRIX[";
	
	unsigned int nAO	= molecule->num_of_ao;
	int in_indx			= -1;
	int fin_indx		= ( ( nAO*(nAO+1) )/2 ) ;
	int nLines 			= 0;
	double temp			= 0.0;
	
	char tmp_line[500];
	
	if ( IF_file(name_f) ){
		std::ifstream buf(name_f);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			Iline Line(tmp_line);
			if ( in_indx == -1) {
				if ( Line.IF_word( keyword,0,keyword.size() ) ) {
					in_indx		= nLines;
				}
			}else if ( in_indx >0 && molecule->m_overlap.size() < fin_indx) {
				std::stringstream ssline(tmp_line);
				while ( ssline >> temp ){
					molecule->m_overlap.push_back(temp);
				}
			}
			nLines++;
		}
		buf.close();
		m_log->input_message("Size of Overlap matrix:");
		m_log->input_message( int( molecule->m_overlap.size() ) );
	}else{
		m_log->input_message("Nothing to read!\n");
		nLines = 0;
	}
}
/***************************************************************************************/
void mopac_files::get_mo(bool beta){
	
	string keyword		= "EIGENVECTORS[";
	if (LMO)keyword		= "LMO_VECTORS[";
	if (!RHF){keyword	= "ALPHA_EIGENVECTORS[";}
	if (beta){keyword	= "BETA_EIGENVECTORS[";}
	
	unsigned int nLines	= 0;
	double temp 		= 0.0;
	unsigned int nMO	= 0;
	unsigned int nAO	= molecule->num_of_ao;
	int nHomo			= molecule->num_of_electrons/2 - 1;
	int nLumo			= nHomo+1;
	unsigned int nMO_out= nAO*nAO;
	int in_indx			= -1;
	int fin_indx		= nAO*nAO/10;
	std::vector<double> mo_c;

	char tmp_line[500];
	string tmpt;
	
	if ( IF_file(name_f) ){
		std::ifstream buf(name_f);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			Iline Line(tmp_line);
			if ( in_indx == -1 ){
				if ( Line.IF_word( keyword,0,keyword.size() ) ){
					in_indx = nLines;
				}
			}else if ( in_indx > 0 && nMO < nMO_out){
				std::stringstream ssline(tmp_line);
				while ( ssline >> temp ){
					mo_c.push_back(temp);
					nMO++;
				}
			}else{
				break;
			}
			nLines++;
		}
		buf.close();
		if ( !beta ) {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule->coeff_MO) );
			m_log->input_message("Number of MO vectors:");
			m_log->input_message( int( molecule->coeff_MO.size() ) );
		}
		else {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule->coeff_MO_beta) );
			m_log->input_message("Number of beta MO vectors:");
			m_log->input_message( int( molecule->coeff_MO_beta.size() ) );
		}
	}else{
		string message = "Not possible to open the file: ";
		message += name_f;
		message +="\n";
		cout << message << endl;
		m_log->input_message(message);
	}
}
/***************************************************************************************/
void mopac_files::get_mo_energies(bool beta){
	
	string keyword			= "EIGENVALUES[";
	if ( !RHF ) { keyword	= "ALPHA_EIGENVALUES[";}
	if ( beta ) { keyword 	= "BETA_EIGENVALUES[";}
	if ( LMO) { keyword		= "LMO_ENERGY_LEVELS[";}
	
	unsigned int nLines		= 0;
	double temp 			= 0.0;
	unsigned int nMO		= 0;
	unsigned int nAO		= molecule->num_of_ao;
	int in_indx				= -1;
	int fin_indx			= nAO/10;
	std::vector<double> mo_c;
	
	char tmp_line[500];
	string tmpt;
	
	if ( IF_file(name_f) ){
		std::ifstream buf(name_f);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			Iline Line(tmp_line);
			if ( in_indx == -1 ){
				if ( Line.IF_word( keyword,0,keyword.size() ) ) {
					in_indx = nLines;
					fin_indx +=  nLines +1;
					tmpt = Line.words[0].substr(keyword.size(),Line.words[0].size());
					tmpt = tmpt.substr(0,tmpt.size()-2);
					if ( nAO != stoi(tmpt)  ){
						m_log->input_message("The molecular orbitals set are imcomplete!");
					}
				}
			}else if ( in_indx > 0 && mo_c.size() < nAO ){
				std::stringstream ssline(tmp_line);
				while ( ssline >> temp ){
					mo_c.push_back(temp);
				}
				nMO++;
			}
			nLines++;
		}
		buf.close();
		if ( !beta ) {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule->orb_energies) );
			molecule->MOnmb = molecule->orb_energies.size();
			m_log->input_message("Number of MO energy levels:");
			m_log->input_message( int( molecule->orb_energies.size() ) );
			if ( molecule->MOnmb < molecule->num_of_atoms ){
				m_log->input_message("Problem in reading molecular energies!");
				m_log->input_message("Number of MO energy levels:");
				m_log->input_message( int( molecule->orb_energies.size() ) );
				if ( LMO ){
					m_log->input_message("mozyme run");
					m_log->input_message(keyword);
				}
			}
		}
		else {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule->orb_energies_beta) );
			m_log->input_message("Number of beta MO energy levels:");
			m_log->input_message( int( molecule->orb_energies_beta.size() ) );
			molecule->MOnmb_beta = molecule->orb_energies_beta.size();
		}
		
	}else{
		string message = "Not possible to open the file: ";
		message += name_f;
		message +="\n";
		cout << message << endl;
		m_log->input_message(message);
	}

}
/***************************************************************************************/
mopac_files::~mopac_files(){}
//================================================================================
//END OF FILE
//================================================================================