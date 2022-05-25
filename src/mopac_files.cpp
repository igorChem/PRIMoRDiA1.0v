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

/**************************************/
// Keywords	

string _keywords= "KEYWORDS";					//0
string _mozyme	= "mozyme";						//1
string _MOZYME	= "MOZYME";						//2
string _n_elec  = "NUM_ELECTRONS";				//3
string _overlap = "OVERLAP_MATRIX[";			//4
string _atom_xa	= "ATOM_X:ANGSTROMS[";			//5
string _hof		= "HEAT_OF_FORMATION:";			//6
string _num_a_el= "NUM_ALPHA_ELECTRONS=";		//7
string _num_b_el= "NUM_BETA_ELECTRONS=";		//8
string _elec_en = "ENERGY_ELECTRONIC:EV=";		//9
string _tot_en	= "TOTAL_ENERGY:EV=";			//9
string _chg		= "CHARGE";						//10
string _orbital = "ORBITAL";					//11
string _inv_mat = "INVERSE_MATRIX[";			//12
string _keywords2= "Keywords:";					//13
string _atom_el	= "ATOM_EL[";					//14
string _atom_cor= "ATOM_CORE[";					//15
string _atom_ao = "AO_ATOMINDEX[";				//16
string _sym_type= "ATOM_SYMTYPE";				//17
string _ao_zeta = "AO_ZETA[";					//18
string _ao_pqn  = "ATOM_PQN[";					//19
string _atom_chg = "ATOM_CHARGES[";				//10
string _gradients = "GRADIENTS:KCAL/MOL/ANGSTROM[";

vector<string> _states = {"SINGLET", "DOUBLET", "TRIPLET", "QUARTET","QUINTET"};
/**************************************/

//======================================================================
mopac_files::mopac_files()		:
	RHF(true)					,
	LMO(false)					,
	f_chg(0)					,
	name_f("no_name")			,
	is_open(false)				,
	type("no_type")				{
}
/***************************************************************************************/
mopac_files::mopac_files(const char* file_name):
	LMO(false)									,
	is_open(false)								,
	RHF(true)									{
	
	name_f = file_name;
	molecule.name = get_file_name(file_name);
	molecule.name = remove_extension( molecule.name.c_str() );
	
	
	if ( IF_file( file_name ) ){
		is_open = true;
		//Ibuffer Buffer(file_name,true);
		if ( check_file_ext(".aux",file_name ) ) {
			type = "AUX";
			Ibuffer Buffer(file_name,_keywords,_overlap);
			for( unsigned i=0; i<Buffer.lines.size(); i++ ){
				if ( Buffer.lines[i].IF_word(_keywords,0,8 ) ){
					for( unsigned j=0; j<Buffer.lines[i].words.size(); j++ ){
						if ( Buffer.lines[i].IF_word(_mozyme,j,6) ){
							LMO = true;
							RHF = true;
						}
						else if ( Buffer.lines[i].IF_word(_MOZYME,j,6) ){
							LMO = true;
							RHF = true;
						}
					}
				}
				else if ( Buffer.lines[i].IF_word(_n_elec,0,13) ){
					string tmp_a(Buffer.lines[i].words[0],14,Buffer.lines[0].words[0].size()-14);
					molecule.num_of_electrons = stoi(tmp_a);
				}
				else if ( Buffer.lines[i].IF_word(_hof,0,18) ){
					string tmp_ab(Buffer.lines[i].words[0],27,14 );
					molecule.heat_of_formation = D_E_conv(tmp_ab);
				}
				else if ( Buffer.lines[i].IF_word(_num_a_el,0,20) ){
					string tmp_b(Buffer.lines[i].words[0],20,Buffer.lines[i].words[0].size()-20);
					if ( stoi(tmp_b) > 0 ) {
						RHF = false;
						molecule.num_of_electrons = 0;
						molecule.num_of_electrons += stoi(tmp_b);
						molecule.betad = true;
					}
				}else if ( Buffer.lines[i].IF_word(_num_b_el,0,19) ){
					string tmp_c(Buffer.lines[i].words[0],19,Buffer.lines[i].words[0].size()-19);
					if ( stoi(tmp_c) > 0 ) {
						molecule.num_of_electrons += stoi(tmp_c);
					}
				}
				else if ( Buffer.lines[i].IF_word(_elec_en,0,21) ){
					string word1 ( Buffer.lines[i].words[0],21,14);
					molecule.elec_energy = D_E_conv(word1);
					molecule.elec_energy *= 0.0367493;
				}
				else if ( Buffer.lines[i].IF_word(_tot_en,0,16) ){
					string word2 ( Buffer.lines[i].words[0],16,14);
					molecule.energy_tot = D_E_conv(word2);
					molecule.energy_tot *= 0.0367493;
				}
			}
		}else if ( check_file_ext(".out",file_name ) ) {
			type = "OUT";
			Ibuffer Buffer_o(name_f,14,25);
			for( unsigned i=0; i<Buffer_o.nLines; i++ ){
				if 	( Buffer_o.lines[i].IF_word(_MOZYME,1) ){
					LMO	= true;
					RHF = true;
				}
				else if ( Buffer_o.lines[i].IF_word(_states[0],1,7) ) RHF =true; 
				else if ( Buffer_o.lines[i].IF_word(_states[1],1,7) ) RHF =false; 
				else if ( Buffer_o.lines[i].IF_word(_states[2],1,7) ) RHF =false; 
				else if ( Buffer_o.lines[i].IF_word(_states[3],1,7) ) RHF =false; 
				else if ( Buffer_o.lines[i].IF_word(_states[4],1,7) ) RHF =false; 
				else if ( Buffer_o.lines[i].IF_word(_chg,1,6) )	f_chg = stoi(Buffer_o.lines[i].words[5]); 
			}
		}else if (check_file_ext(".mgf",file_name) ){
			type = "MGF";
		}else{
			cout << "Warning! The file has imcompatible extension name with mopac files!" << endl;
			m_log->write_warning("Warning! The file has imcompatible extension name with mopac files!");
			is_open = false;
		}
	}else{
		m_log->write_error("Error opening MOPAC file! Verify its presence in the current directory!");
		cout << "the file named " << file_name << " cannot be opened! " << endl;
	}
}
/***************************************************************************************/
void mopac_files::parse_aux(){
	
	m_log->input_message("Starting to parse MOPAC AUX file.\n");
	vector<int> aoidx;
	vector<double> zetas;
	vector<int> shells;
	
	vector<unsigned> _in(7);
	vector<unsigned> _out(7);
	
	Ibuffer Buffer(name_f,_keywords,_overlap);	
	_out[6] = Buffer.nLines-1;
	for( unsigned i=0; i<Buffer.nLines-1;i++){
		if ( Buffer.lines[i].IF_word(_atom_el,0,8) ){
			_in[0] = i+1;
		}else if ( Buffer.lines[i].IF_word(_atom_cor,0,10) ){
			_out[0] = i;
		}else if ( Buffer.lines[i].IF_word(_atom_xa,0,17) ){
			_in[1] = i+1;
		}else if ( Buffer.lines[i].IF_word(_atom_ao,0,13) ){
			_out[1] = i;
			_in[2] = i+1;
		}else if ( Buffer.lines[i].IF_word(_sym_type,0,12) ){
			_out[2] = i;
			_in[3] = i+1;
		}else if ( Buffer.lines[i].IF_word(_ao_zeta,0,8) ){
			_out[3] = i;
			_in[4] = i+1;
		}else if ( Buffer.lines[i].IF_word(_ao_pqn,0,9) ){
			_out[4] = i;
			_in[5] = i+1;
		}else if ( Buffer.lines[i].IF_word(_n_elec,0,13) ){
			_out[5] = i;			
		}else if ( Buffer.lines[i].IF_word(_atom_chg,0,13) ){
			_in[6] = i+1;
		}
		else if ( Buffer.lines[i].IF_word(_gradients,0,28) ) {
			_out[6]=i+1;
		}
	}
	
	
	for( unsigned int i=_in[0]; i<_out[0]; i++ ){
		for ( unsigned j=0; j<Buffer.lines[i].words.size(); j++ ){
				Iatom atom;
				atom.set_type( Buffer.lines[i].words[j] );
				molecule.add_atom(atom);
		}
	}
	
	unsigned counter = 0;
	
	for (unsigned i=_out[0]+1; i<_in[1]-1; i++ ){
		for ( unsigned j=0;j<Buffer.lines[i].words.size();j++){
			molecule.atoms[counter++].atomicN  = Buffer.lines[i].get_double(0);
		}
	}
	
	counter = 0;
	
	for (unsigned i=_in[1]; i<_out[1]; i++ ){
		molecule.atoms[counter].xcoord  = Buffer.lines[i].get_double(0);
		molecule.atoms[counter].ycoord 	= Buffer.lines[i].get_double(1);
		molecule.atoms[counter++].zcoord = Buffer.lines[i].get_double(2);
	}
	
	for ( unsigned i=_in[2]; i<_out[2]; i++ ){
		for ( unsigned j=0; j<Buffer.lines[i].words.size(); j++ ){ 
			aoidx.push_back( Buffer.lines[i].get_int(j) );
		}
	}
	
	counter = 0;
	
	for (unsigned i=_in[3]; i<_out[3];i++ ){
		for ( unsigned j=0; j<Buffer.lines[i].words.size();j++){
			Iaorbital aorb;
			aorb.symmetry = Buffer.lines[i].get_string(j);
			if 		( aorb.symmetry == "PX" ) aorb.powx	= 1;
			else if ( aorb.symmetry == "PY" ) aorb.powy	= 1;
			else if ( aorb.symmetry == "PZ" ) aorb.powz = 1;
			molecule.atoms[aoidx[counter++]-1].add_orbital(aorb);
		}
	}

	for (unsigned i=_in[4]; i<_out[4]; i++){
		for ( unsigned j=0; j<Buffer.lines[i].words.size(); j++){
			zetas.push_back( Buffer.lines[i].get_double(j) );
		}
	}
	
	for (unsigned i=_in[5];i<_out[5];i++){
		for (unsigned j=0;j<Buffer.lines[i].words.size();j++){
			shells.push_back( Buffer.lines[i].get_int(j) );		
		}
	}
	
	counter = 0;
	
	for (unsigned i=_in[6]; i<_out[6]; i++){
		for (unsigned j=0;j<Buffer.lines[i].words.size();j++){
			molecule.atoms[counter++].charge = Buffer.lines[i].get_double(j);
		}
	}
		
	if ( molecule.atoms.size() == 0 ){ 
		cout << "Warning! Zero atoms read in aux file. Verify your file!!" << endl;
	} 

	m_log->input_message("Found number of atoms in the aux file: \n\t");
	m_log->input_message( int(molecule.num_of_atoms) );
	m_log->input_message("\n");
	
	m_log->input_message("Number of electron in the system: \n\t");
	m_log->input_message( int(molecule.num_of_electrons) );
	m_log->input_message("\n");
	
	m_log->input_message("Number of atomic orbitals: \n\t");
	m_log->input_message( int( aoidx.size() ) );
	m_log->input_message("\n"); 
	
	
	molecule.get_ao_number();
	
	this->get_overlap_m();
	this->get_mo(false);
	this->get_mo_energies(false);
	if ( !RHF ) {
		this->get_mo(true);
		this->get_mo_energies(true);
	}	
	
	counter = 0;
	
	for (unsigned i=0;i<molecule.atoms.size();i++){
		for ( unsigned j=0;j<molecule.atoms[i].norb;j++){
			molecule.atoms[i].orbitals[j].alpha = zetas[counter];
			molecule.atoms[i].orbitals[j].shell = shells[counter++];
		}
	}
	molecule.update();
	if ( !molecule.check() ) { 
		cout << "Problems in reading the mopac aux file: " << name_f << endl;
	}
	m_log->input_message("HOMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.homo_energy) );
	m_log->input_message("\nHOMO number: \n\t");
	m_log->input_message( std::to_string(molecule.homoN) );
	m_log->input_message("\nLUMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.lumo_energy) );
	m_log->input_message("\nLUMO number: \n\t");
	m_log->input_message( std::to_string(molecule.lumoN) );
	m_log->input_message("\n");
}
/***************************************************************************************/
void mopac_files::parse_out(){
	m_log->input_message("Starting to parse out file from MOPAC.\n");
	
	Ibuffer Buffer(name_f,true) ;
	for (unsigned  i=0;i<Buffer.nLines;i++){
		if ( Buffer.lines[i].IF_line("HEAT",1,"FORMATION",3,10) ){ 
			 molecule.heat_of_formation = Buffer.lines[i].pop_double(5);
		}
		else if ( Buffer.lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,5) || Buffer.lines[i].IF_line("ELECTRONIC",0,"ENERGY",1,8) ){
			molecule.energy_tot = Buffer.lines[i].pop_double(3);
			//cout << molecule.energy_tot << endl;
		}
		else if ( Buffer.lines[i].IF_line("HOMO",0,"LUMO",1,7) ){
			molecule.homo_energy = Buffer.lines[i].pop_double(5);
			molecule.lumo_energy = Buffer.lines[i].pop_double(5);
		}
		else if ( Buffer.lines[i].IF_line("SUM",0,5) ){
			for ( int j=0;j<3;j++){ molecule.dipole_moment[j] = Buffer.lines[i].pop_double(1); }
			molecule.total_dipmoment = Buffer.lines[i].pop_double(1);
		}
	}

	m_log->input_message("Total Energy: \n\t");
	m_log->input_message( std::to_string(molecule.energy_tot) );
	m_log->input_message("\nHOMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.homo_energy) );
	m_log->input_message("\nLUMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.lumo_energy) );
	m_log->input_message("\n");
	m_log->inp_delim(1);
}
/***************************************************************************************/
void mopac_files::parse_mgf(){
	vector<double>	zetasS;
	vector<double>	zetasP;
	vector<double>	zetasD;
	vector<double>	inv_mat;
	vector<int>		orbN;
	vector<int>		orbN_beta;
	
	int inmat_in		= 0;
	int inmat_fin		= 0;
	int noa				= 0;
	double tot_charge	= 0;
 	
	m_log->input_message("Starting to parse MOPAC MGF file!\n");
	
	string mgf_name = change_extension(name_f,".mgf");
	
	if ( IF_file( mgf_name.c_str() ) ){
	
		Ibuffer Buffer( mgf_name.c_str(),true );
		for( unsigned i=0; i<Buffer.nLines; i++ ){
			if( i == 0 ){
				noa = Buffer.lines[i].pop_int(0);
			}
			else if ( i>0 && i<=noa ) {
				Iatom atom;
				atom.set_type( get_atomic_symbol( Buffer.lines[i].pop_int(0) ) );
				atom.xcoord  = Buffer.lines[i].pop_double(0);
				atom.ycoord  = Buffer.lines[i].pop_double(0);
				atom.zcoord  = Buffer.lines[i].pop_double(0);
				atom.charge  = Buffer.lines[i].pop_double(0);
				molecule.add_atom(atom);
			}
			else if( i>noa && i<=(noa*2) ){
				zetasS.push_back( Buffer.lines[i].pop_double(0) );
				zetasP.push_back( Buffer.lines[i].pop_double(0) );
				zetasD.push_back( Buffer.lines[i].pop_double(0) );
			}
			else if( Buffer.lines[i].IF_word(_orbital,0,7) && inmat_fin == 0 ){
				orbN.push_back(i);
			}
			else if( Buffer.lines[i].IF_word(_inv_mat,0,15) ){
				inmat_in = i;				
			}
			else if( Buffer.lines[i].IF_word(_keywords2,0,9) ) {
				inmat_fin = i;				
			}
			else if( inmat_fin > 0 && Buffer.lines[i].IF_word(_orbital,0,7) ){
				molecule.betad = true;
				orbN_beta.push_back(i);
			}
		}

		string temp = "noname";
		for( unsigned i=0; i<orbN.size(); i++ ){
			int fin_ind = 0;
			if ( i==orbN.size()-1 )	fin_ind = inmat_in;
			else fin_ind = orbN[i+1];
			for(int j=orbN[i];j<fin_ind;j++){
				if ( j == orbN[i] ){
					molecule.occupied.push_back( Buffer.lines[j].pop_int(1) );
					molecule.orb_energies.push_back( Buffer.lines[j].pop_double(2) );
					molecule.MOnmb++;
				}else{ 
					for(int k=0;k<Buffer.lines[j].line_len;k++){
						if (k == 0 &&  Buffer.lines[j].words[k][0] == '-' ){
							if ( Buffer.lines[j].words[k].size() > 15 ){
								temp = Buffer.lines[j].words[k].substr(0,15);
								molecule.coeff_MO.push_back( D_E_conv(temp) );
								temp = Buffer.lines[j].words[k].substr( 15,Buffer.lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule.coeff_MO.push_back( D_E_conv(temp) );
									temp = Buffer.lines[j].words[k].substr( 30,Buffer.lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule.coeff_MO.push_back( D_E_conv(temp) );
										temp = Buffer.lines[j].words[k].substr( 45,Buffer.lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule.coeff_MO.push_back( D_E_conv(temp) );
											temp = Buffer.lines[j].words[k].substr( 60,Buffer.lines[j].words[k].size() );
											molecule.coeff_MO.push_back( D_E_conv(temp) );
										}else molecule.coeff_MO.push_back( D_E_conv(temp) );
									}else molecule.coeff_MO.push_back( D_E_conv(temp) );
								}else molecule.coeff_MO.push_back( D_E_conv(temp) );
							}else molecule.coeff_MO.push_back( D_E_conv( Buffer.lines[j].words[k]) );
						}else{
							if ( Buffer.lines[j].words[k].size() > 14 ){
								temp = Buffer.lines[j].words[k].substr(0,14);
								molecule.coeff_MO.push_back( D_E_conv(temp) );
								temp = Buffer.lines[j].words[k].substr( 14,Buffer.lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule.coeff_MO.push_back( D_E_conv(temp) );
									temp = Buffer.lines[j].words[k].substr( 29,Buffer.lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule.coeff_MO.push_back( D_E_conv(temp) );
										temp = Buffer.lines[j].words[k].substr( 44,Buffer.lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule.coeff_MO.push_back( D_E_conv(temp) );
											temp = Buffer.lines[j].words[k].substr(59,Buffer.lines[j].words[k].size() );
											molecule.coeff_MO.push_back( D_E_conv(temp) );
										}else molecule.coeff_MO.push_back( D_E_conv(temp) );
									}else molecule.coeff_MO.push_back( D_E_conv(temp) );
								}else molecule.coeff_MO.push_back( D_E_conv(temp) );
							}else molecule.coeff_MO.push_back( D_E_conv(Buffer.lines[j].words[k]) );
						}	
					}
				}
			}
		}
		for( int i=0; i<noa; i++ ){
			int sh = 0;
			Iaorbital orbS;
			orbS.alpha = zetasS[i];
			orbS.shell = 1;
			if 		( molecule.atoms[i].atomicN > 1  && molecule.atoms[i].atomicN < 11 ) sh		= orbS.shell = 2;
			else if ( molecule.atoms[i].atomicN > 10 && molecule.atoms[i].atomicN < 19 ) sh	= orbS.shell = 3;
			else if ( molecule.atoms[i].atomicN > 18 && molecule.atoms[i].atomicN < 37 ) sh	= orbS.shell = 4;
			molecule.atoms[i].add_orbital(orbS);
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
				molecule.atoms[i].add_orbital(orbPx);
				molecule.atoms[i].add_orbital(orbPy);
				molecule.atoms[i].add_orbital(orbPz);
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
				molecule.atoms[i].add_orbital(orbDxx);
				molecule.atoms[i].add_orbital(orbDyy);
				molecule.atoms[i].add_orbital(orbDzz);
				molecule.atoms[i].add_orbital(orbDxy);
				molecule.atoms[i].add_orbital(orbDyz);
				molecule.atoms[i].add_orbital(orbDxz);
			}
		}
	
		for(int i=inmat_in+1;i<inmat_fin;i++){
			for(int k=0;k<Buffer.lines[i].line_len;k++){
				temp = Buffer.lines[i].words[k];
				if ( k == 0 && temp[0] == '-' ){
					if ( temp.size() > 15 ){ 
						temp = temp.substr(0,15);
						inv_mat.push_back( D_E_conv(temp) );
						temp = Buffer.lines[i].words[k].substr( 15,Buffer.lines[i].words[k].size() );
						if ( temp.size() > 15 ){
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = Buffer.lines[i].words[k].substr( 30,Buffer.lines[i].words[k].size() );
							if ( temp.size() > 15) {
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = Buffer.lines[i].words[k].substr( 45,Buffer.lines[i].words[k].size() );
								if ( temp.size() > 15 ){
									temp = temp.substr(0,15);
									inv_mat.push_back( D_E_conv(temp) );
									temp = Buffer.lines[i].words[k].substr( 60,Buffer.lines[i].words[k].size() );
									inv_mat.push_back( D_E_conv(temp) );
								}else inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) );
				}else{
					if( temp.size() > 14 ){
						temp = temp.substr(0,14);
						inv_mat.push_back( D_E_conv(temp) );
						temp = Buffer.lines[i].words[k].substr( 14,Buffer.lines[i].words[k].size() );
						if ( temp.size() > 15 ){
							temp = temp.substr(0,15);
							inv_mat.push_back( D_E_conv(temp) );
							temp = Buffer.lines[i].words[k].substr( 29,Buffer.lines[i].words[k].size() );
							if ( temp.size() > 15 ) {
								temp = temp.substr(0,15);
								inv_mat.push_back( D_E_conv(temp) );
								temp = Buffer.lines[i].words[k].substr( 44,Buffer.lines[i].words[k].size() );
								if ( temp.size() > 15 ){
									temp = temp.substr(0,15);
									inv_mat.push_back( D_E_conv(temp) );
									temp = Buffer.lines[i].words[k].substr( 59,Buffer.lines[i].words[k].size() );
									inv_mat.push_back( D_E_conv(temp) );
								}else inv_mat.push_back( D_E_conv(temp) );
							}else inv_mat.push_back( D_E_conv(temp) );
						}else inv_mat.push_back( D_E_conv(temp) );
					}else inv_mat.push_back( D_E_conv(temp) ); 
				}
			}
		}	
		unsigned  nmo = molecule.MOnmb;
		int k   = 0;
	
	
		for( int i=0; i<orbN_beta.size(); i++ ){
			int fin_ind = 0;
			if ( i==orbN_beta.size()-1 ) fin_ind = Buffer.nLines;
			else fin_ind = orbN_beta[i+1];
			for(int j=orbN_beta[i];j<fin_ind;j++){
				if ( j == orbN_beta[i] ){
					molecule.occupied_beta.push_back( Buffer.lines[j].pop_int(1) );
					molecule.orb_energies_beta.push_back( Buffer.lines[j].pop_double(2) );
					molecule.MOnmb_beta++;
				}else{
					for(int k=0;k<Buffer.lines[j].line_len;k++){
						if (k == 0 &&  Buffer.lines[j].words[k][0] == '-' ){
							if ( Buffer.lines[j].words[k].size() > 15 ){
								temp = Buffer.lines[j].words[k].substr(0,15);
								molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = Buffer.lines[j].words[k].substr( 15,Buffer.lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = Buffer.lines[j].words[k].substr( 30,Buffer.lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = Buffer.lines[j].words[k].substr( 45,Buffer.lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
											temp = Buffer.lines[j].words[k].substr(60,Buffer.lines[j].words[k].size() );
											molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
										}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule.coeff_MO_beta.push_back( D_E_conv(Buffer.lines[j].words[k]) );
						}else{
							if ( Buffer.lines[j].words[k].size() > 14 ){
								temp = Buffer.lines[j].words[k].substr(0,14);
								molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
								temp = Buffer.lines[j].words[k].substr( 14,Buffer.lines[j].words[k].size() );
								if ( temp.size() > 15) {
									temp = temp.substr(0,15);
									molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
									temp = Buffer.lines[j].words[k].substr( 29,Buffer.lines[j].words[k].size() );
									if ( temp.size() > 15) {
										temp = temp.substr(0,15);
										molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
										temp = Buffer.lines[j].words[k].substr( 44,Buffer.lines[j].words[k].size() );
										if ( temp.size() > 15 ){
											temp = temp.substr(0,15);
											molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
											temp = Buffer.lines[j].words[k].substr(59,Buffer.lines[j].words[k].size() );
											molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
										}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
									}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
								}else molecule.coeff_MO_beta.push_back( D_E_conv(temp) );
							}else molecule.coeff_MO_beta.push_back( D_E_conv(Buffer.lines[j].words[k]) );
						}		
					}
				}
			}
		}
	
		Eigen::MatrixXd inv(nmo,nmo);
		Eigen::MatrixXd coeff_corrected(nmo,nmo);

		int bnmo = molecule.MOnmb_beta;
 
		Eigen::MatrixXd coeff_corrected_beta(bnmo,bnmo);
	
		for ( unsigned i=0;i<nmo;i++){
			for (unsigned j=0;j<nmo;j++){  coeff_corrected(i,j) = molecule.coeff_MO[i*nmo+j]; }
		}
	
		for (unsigned i=0;i<bnmo;i++){
			for (unsigned j=0;j<bnmo;j++) coeff_corrected_beta(i,j) = molecule.coeff_MO_beta[i*nmo+j];
		}
	
		for(unsigned i=0;i<nmo;i++){
			for(unsigned j=0;j<=i;j++){
				if (i==j)  inv(i,j) = inv_mat[k++];
				else       inv(j,i) = inv(i,j) = inv_mat[k++];
			}
		}
	
		coeff_corrected = coeff_corrected*inv;
		for (unsigned int i=0;i<nmo;i++){
			for (unsigned int j=0;j<nmo;j++)  molecule.coeff_MO[i*nmo+j] = coeff_corrected(i,j);
		}
		
		if ( molecule.betad ){
			coeff_corrected_beta = coeff_corrected_beta*inv; 
			for (unsigned i=0;i<nmo;i++){
				for (unsigned j=0;j<nmo;j++)	molecule.coeff_MO_beta[i*nmo+j] = coeff_corrected_beta(i,j);
			}
		}
		
		for(unsigned i=0;i<molecule.occupied.size();i++){
			molecule.num_of_electrons+= molecule.occupied[i];
		}
		
		molecule.update();
		
		m_log->input_message("HOMO energy: \n\t");
		m_log->input_message( std::to_string(molecule.homo_energy) );
		m_log->input_message("\nLUMO energy: \n\t");
		m_log->input_message( std::to_string(molecule.lumo_energy) );
		m_log->input_message("\nAtomic orbitals: \n\t");
		m_log->input_message(molecule.get_ao_number() );
		m_log->input_message("\n");
		m_log->inp_delim(2);
		
		for(unsigned i=0;i<molecule.atoms.size();i++){
			if (molecule.atoms[i].atomicN > 10 ){
				molecule.atoms[i].atomicN -= 10;
			}
			else if ( molecule.atoms[i].atomicN > 2 ){
				molecule.atoms[i].atomicN -= 2;
			}
		}
		
	}else{
		m_log->write_error("MGF file does not exist, skipping its parsing process!");
	}
}
/***************************************************************************************/
void mopac_files::get_overlap_m(){
	
	unsigned int nAO	= molecule.num_of_ao;
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
			if ( nLines == 205 ){
				int o = 0;
			}
			if ( in_indx == -1) {
				if ( Line.IF_word( _overlap,0,_overlap.size() ) ) {
					in_indx		= nLines;
				}
			}else if ( in_indx >0 && molecule.m_overlap.size() < fin_indx) {
				std::stringstream ssline(tmp_line);
				while ( ssline >> temp ){
					molecule.m_overlap.push_back(temp);
				}
			}
			nLines++;
		}
		buf.close();
		m_log->input_message("Size of Overlap matrix: \n\t");
		m_log->input_message( int( molecule.m_overlap.size() ) );
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
	unsigned int nAO	= molecule.num_of_ao;
	int nHomo			= molecule.num_of_electrons/2 - 1;
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
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule.coeff_MO) );
			m_log->input_message("Number of MO vectors: \n\t");
			m_log->input_message( int( molecule.coeff_MO.size() ) );
		}
		else {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule.coeff_MO_beta) );
			m_log->input_message("Number of beta MO vectors: \n\t");
			m_log->input_message( int( molecule.coeff_MO_beta.size() ) );
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
	if ( LMO )	{ keyword	= "LMO_ENERGY_LEVELS[";}
	
	unsigned int nLines		= 0;
	double temp 			= 0.0;
	unsigned int nMO		= 0;
	unsigned int nAO		= molecule.num_of_ao;
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
						m_log->write_warning("The molecular orbitals set are imcomplete!");
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
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule.orb_energies) );
			molecule.MOnmb = molecule.orb_energies.size();
			m_log->input_message("Number of MO energy levels: \n\t");
			m_log->input_message( int( molecule.orb_energies.size() ) );
			if ( molecule.MOnmb < molecule.num_of_atoms ){
				m_log->input_message("Problem in reading molecular energies! \n\t");
				m_log->input_message("Number of MO energy levels: \n\t");
				m_log->input_message( int( molecule.orb_energies.size() ) );
				if ( LMO ){
					m_log->input_message("mozyme run ");
					m_log->input_message(keyword);
					m_log->input_message("\n");
				}
			}
		}
		else {
			copy( mo_c.begin(),mo_c.end(),back_inserter(molecule.orb_energies_beta) );
			m_log->input_message("Number of beta MO energy levels: \n\t");
			m_log->input_message( int( molecule.orb_energies_beta.size() ) );
			molecule.MOnmb_beta = molecule.orb_energies_beta.size();
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