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
#include <algorithm> // review the need  

//Including PRIMoRDiA headers
//-------------------------------------------------------
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
#include "../include/gaussian_files.h"

using std::vector;
using std::string;
using std::stoi;
using std::stod;
using std::endl;
using std::cout;
using std::unique_ptr;

//================================================
gaussian_files::gaussian_files():
	name_f("noname")			,
	is_open(false)				{
}
/**********************************************************************/
gaussian_files::gaussian_files(const char* file_name):
	name_f(file_name)								,
	is_open(false)									{
	
	if ( IF_file(file_name) ) {
		is_open = true;
		if ( !check_file_ext(".fchk",name_f) ) {
			cout << "Warning! The file has wrong etension name!" << endl;
			m_log->write_warning("Warning! The file has wrong etension name!");
			is_open = false;
		}
	}else{
		m_log->write_error("Error opening GAUSSIAN file! Verify its presence in the current directory!");
		cout << "the file named " << file_name << " cannot be opened! " << endl;
	}
}
/**********************************************************************/
void gaussian_files::parse_fchk(){
	
	m_log->input_message("Starting to parse GAUSSIAN FCHK file!\n");
	molecule.name = get_file_name(name_f);
	molecule.name = remove_extension( molecule.name.c_str() );
	
	unsigned j,k,l;
	unsigned atomic_n_i = 0;
	unsigned atomic_n_f	= 0;
	unsigned coord_n_i	= 0;
	unsigned coord_n_f	= 0;
	unsigned shell_t_i	= 0;
	unsigned shell_t_f	= 0;
	unsigned prim_b_i	= 0;
	unsigned prim_b_f	= 0;
	unsigned shell_m_i	= 0;
	unsigned shell_m_f	= 0;
	unsigned prim_e_i	= 0;
	unsigned prim_e_f	= 0;
	unsigned cont_c_i	= 0;
	unsigned cont_c_f	= 0;
	unsigned conspt_c_i	= 0;
	unsigned conspt_c_f	= 0;
	unsigned alpha_e_i 	= 0;
	unsigned alpha_e_f	= 0;
	unsigned beta_e_i	= 0;
	unsigned beta_e_f	= 0;
	unsigned alpha_c_i	= 0;
	unsigned alpha_c_f	= 0;
	unsigned beta_c_i	= 0;
	unsigned beta_c_f	= 0;
	unsigned chgs_i		= 0;
	unsigned chgs_f		= 0;
	unsigned dens_i		= 0;
	unsigned dens_f		= 0;
	
	vector<double>	coords;
	vector<int>		shell_t;
	vector<int>		ngtos;
	vector<int>		shell_map;
	vector<double>	cont_c;
	vector<double>	cont_c_p;
	vector<double>	expos;
	
	int mult = 1;
	
	Ibuffer Buffer(name_f,true);
	for( unsigned i=0; i<Buffer.lines.size(); i++ ){
		if ( Buffer.lines[i].IF_line("Multiplicity",0,"I",1,3) ) { mult = Buffer.lines[i].pop_int(2); }
		if ( Buffer.lines[i].IF_line("Number",0,"electrons",2,5) ) { molecule.num_of_electrons = Buffer.lines[i].pop_int(4); } 
		else if ( Buffer.lines[i].IF_line("Number",0,"beta",2,6) ) {
			if ( mult%2 == 0 ) {
				molecule.betad = true;
			}
		}
		else if ( Buffer.lines[i].IF_line("Number",0,"beta",2,5) ) { molecule.num_of_electrons = Buffer.lines[i].pop_int(4); } 
		else if ( Buffer.lines[i].IF_line("Atomic",0,"numbers",1,5) ){ atomic_n_i = i;	}
		else if ( Buffer.lines[i].IF_line("Nuclear",0,"charges",1,5) ){ atomic_n_f = i; }
		else if ( Buffer.lines[i].IF_line("cartesian",1,"coordinates",2,6) ){ coord_n_i = i; }
		else if ( Buffer.lines[i].IF_line("Force",0,"Field",1,4) ){ if ( coord_n_f == 0 ) coord_n_f =i; }
		else if ( Buffer.lines[i].IF_line("Shell",0,"types",1,5) ){ shell_t_i = i; }
		else if ( Buffer.lines[i].IF_line("primitives",2,"shell",4,8) ){ shell_t_f = prim_b_i = i; }
		else if ( Buffer.lines[i].IF_line("Shell",0,"map",3,7) ){ shell_m_i = prim_b_f = i; }
		else if ( Buffer.lines[i].IF_line("Primitive",0,"exponents",1,5) ){ shell_m_f = prim_e_i = i; }
		else if ( Buffer.lines[i].IF_line("Contraction",0,"coefficients",1,5) ){ prim_e_f = cont_c_i = i; }
		else if ( Buffer.lines[i].IF_line("Contraction",1,"coefficients",2,6) ){ conspt_c_i = cont_c_f = i; }
		else if ( Buffer.lines[i].IF_line("Coordinates",0,"shell",3,7) ){ conspt_c_f = i; }		
		else if ( Buffer.lines[i].IF_line("Total",0,"Energy",1,4) ){ molecule.energy_tot = Buffer.lines[i].pop_double(3); }
		else if ( Buffer.lines[i].IF_line("Alpha",0,"Energies",2,6) ){ alpha_e_i = i; }
		else if ( Buffer.lines[i].IF_line("Beta",0,"Energies",2,6) ){ beta_e_i  = alpha_e_f = i; }
		else if ( Buffer.lines[i].IF_line("Alpha",0,"coefficients",2,6) ){
			if ( beta_e_i > 0 ) beta_e_f  = alpha_c_i = i;
			else alpha_e_f = alpha_c_i = i;
		}
		else if ( Buffer.lines[i].IF_line("Beta",0,"coefficients",2,6) ){  beta_c_i  = alpha_c_f = i; }
		else if ( Buffer.lines[i].IF_line("Total",0,"Density",2,6) ){
			dens_i = i;
			if ( beta_c_i > 0 )	beta_c_f = i;
			else alpha_c_f = i;
		}
		else if ( Buffer.lines[i].IF_line("QEq",0,"tensors",2,6) ) {
			if ( dens_f == 0 ){
				dens_f = i;			
			}
		}
		else if ( Buffer.lines[i].IF_line("Mulliken",0,"Charges",1,5) ){ 
			if ( dens_f == 0 ){
				dens_f = i;
			}			
			chgs_i = i;
		}
		else if ( Buffer.lines[i].IF_line("Optimization",0,"MaxStp",1,4) ){ chgs_f = i;}
		else if ( Buffer.lines[i].IF_line("ONIOM",0,"Charges",1,5) ){ if ( chgs_f == 0 ) chgs_f = i;}
	}
	k=0;
	unsigned noe = molecule.num_of_electrons;
	
	for( unsigned i=0; i<Buffer.lines.size(); i++){
		if ( i>atomic_n_i && i<atomic_n_f){
			for( j=0; j<Buffer.lines[i].line_len; j++ ){
				Iatom atom;
				atom.set_type( get_atomic_symbol( Buffer.lines[i].pop_int(0) ) );
				molecule.add_atom(atom);
			}		
		}
		else if ( i>coord_n_i && i<coord_n_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++ ){
				coords.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>shell_t_i && i< shell_t_f){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				shell_t.push_back(Buffer.lines[i].pop_int(0));
			}
		}
		else if( i>prim_b_i && i<prim_b_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				ngtos.push_back(Buffer.lines[i].pop_int(0));
			}
		}	
		else if( i>shell_m_i && i<shell_m_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				shell_map.push_back(Buffer.lines[i].pop_int(0));
			}
		}
		else if( i>prim_e_i && i<prim_e_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				expos.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>cont_c_i && i<cont_c_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				cont_c.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>conspt_c_i && i<conspt_c_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				cont_c_p.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>alpha_e_i && i<alpha_e_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.orb_energies.push_back(Buffer.lines[i].pop_double(0));
				molecule.MOnmb++;
			}
		}
		else if( i>beta_e_i && i<beta_e_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.orb_energies_beta.push_back(Buffer.lines[i].pop_double(0));
				molecule.MOnmb_beta++;
			}
		}
		else if( i>alpha_c_i && i<alpha_c_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.coeff_MO.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>beta_c_i && i<beta_c_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.coeff_MO_beta.push_back(Buffer.lines[i].pop_double(0));
			}
		}
		else if( i>dens_i && i<dens_f){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.m_dens.push_back( Buffer.lines[i].pop_double(0) );
			}
		}
		else if( i>chgs_i && i<chgs_f ){
			for(j=0;j<Buffer.lines[i].line_len;j++){
				molecule.atoms[k++].charge = Buffer.lines[i].pop_double(0);
			}
		}
	}

	j=0;
	for( unsigned i=0; i<molecule.atoms.size(); i++ ){
		molecule.atoms[i].xcoord = coords[j++];
		molecule.atoms[i].ycoord = coords[j++];
		molecule.atoms[i].zcoord = coords[j++];
	}
	k = 0;
	int r = 0;
	for( unsigned i=0; i<shell_map.size(); i++){
		if	( shell_t[i] == 0 ){
			Iaorbital orbS;
			orbS.symmetry = "S";
			for (j=0;j<ngtos[i];j++){ orbS.add_primitive(expos[k++],cont_c[r++]); }
			molecule.atoms[shell_map[i]-1].add_orbital(orbS);
		}
		else if ( shell_t[i] == -1 ){
			vector<double> ex_p;
			vector<double> conc_p; 
			Iaorbital orbS;
			orbS.symmetry = "S";
			for (j=0;j<ngtos[i];j++){
				double tmp_ex = expos[k++];
				double tmp_cc = cont_c[r++];
				orbS.add_primitive(tmp_ex,tmp_cc);
				ex_p.push_back(tmp_ex);
				conc_p.push_back(cont_c_p[k-1]);
			}
			molecule.atoms[shell_map[i]-1].add_orbital(orbS);
			
			Iaorbital orbPX,orbPY,orbPZ;
			for (j=0;j<ngtos[i];j++){ orbPX.add_primitive(ex_p[j],conc_p[j]); }
			orbPY = orbPZ = orbPX;
			orbPX.symmetry = "PX";
			orbPX.powx     = 1;
			orbPY.symmetry = "PY";
			orbPY.powy     = 1;
			orbPZ.symmetry = "PZ";
			orbPZ.powz     = 1;
			molecule.atoms[shell_map[i]-1].add_orbital(orbPX);
			molecule.atoms[shell_map[i]-1].add_orbital(orbPY);
			molecule.atoms[shell_map[i]-1].add_orbital(orbPZ);
		}
		else if ( shell_t[i] == 1 ){
			Iaorbital orbPX,orbPY,orbPZ;
			for (j=0;j<ngtos[i];j++){ orbPX.add_primitive(expos[k++],cont_c[r++]); }
			orbPY = orbPZ = orbPX;
			orbPX.symmetry = "PX";
			orbPX.powx     = 1;
			orbPY.symmetry = "PY";
			orbPY.powy     = 1;
			orbPZ.symmetry = "PZ";
			orbPZ.powz     = 1;
			molecule.atoms[shell_map[i]-1].add_orbital(orbPX);
			molecule.atoms[shell_map[i]-1].add_orbital(orbPY);
			molecule.atoms[shell_map[i]-1].add_orbital(orbPZ);
		}
		else if ( shell_t[i] == 2 ){
			Iaorbital orbDx2,orbDy2,orbDz2,orbDxy,orbDyz,orbDxz;
			for (j=0;j<ngtos[i];j++){ 
				orbDx2.add_primitive(expos[k++],cont_c[r++]);
			}
			orbDy2 = orbDz2 = orbDxy = orbDyz = orbDxz = orbDx2;
			orbDx2.symmetry = "XX";
			orbDx2.powx     = 2;
			orbDy2.symmetry = "YY";
			orbDy2.powy     = 2;
			orbDz2.symmetry = "ZZ";
			orbDz2.powz     = 2;
			orbDxy.symmetry = "XY";
			orbDxy.powx     = 1;
			orbDxy.powy     = 1;
			orbDyz.symmetry = "YZ";
			orbDyz.powy     = 1;
			orbDyz.powz     = 1;
			orbDxz.symmetry = "XZ";
			orbDxz.powx     = 1;
			orbDxz.powz     = 1;	
			molecule.atoms[shell_map[i]-1].add_orbital(orbDx2);
			molecule.atoms[shell_map[i]-1].add_orbital(orbDy2);
			molecule.atoms[shell_map[i]-1].add_orbital(orbDz2);
			molecule.atoms[shell_map[i]-1].add_orbital(orbDxy);
			molecule.atoms[shell_map[i]-1].add_orbital(orbDxz);
			molecule.atoms[shell_map[i]-1].add_orbital(orbDyz);
		}
		else if ( shell_t[i] == 3 ){
			Iaorbital fx3,fy3,fz3,fx2y,fx2z,fy2z,fy2x,fz2x,fz2y,fxyz;
			for (j=0;j<ngtos[i];j++){ 
				fx3.add_primitive(expos[k++],cont_c[r++]);
			}
			fxyz = fz2y = fz2x = fy2x = fy2z = fx2z = fx2y = fz3 = fy3 = fx3;
			fx3.symmetry = "XXX";
			fx3.powx = 3;
			fy3.symmetry = "YYY";
			fy3.powy = 3;
			fz3.symmetry = "ZZZ";
			fz3.powz = 3;
			fx2y.symmetry = "XXY";
			fx2y.powx =2;
			fx2y.powy =1;
			fy2x.symmetry = "YYX";
			fy2x.powx =1;
			fy2x.powy =2;
			fy2z.symmetry = "YYZ";
			fy2z.powz =1;
			fy2z.powy =2;
			fx2z.symmetry = "XXZ";
			fx2z.powx =2;
			fx2z.powz =1;
			fz2x.symmetry = "ZZX";
			fz2x.powx = 1;
			fz2x.powz = 2;
			fz2y.symmetry = "ZZY";
			fz2y.powy = 1;
			fz2y.powz = 2;
			fxyz.symmetry = "XYZ";
			fxyz.powx = 1;
			fxyz.powy = 1;
			fxyz.powz = 1;
			molecule.atoms[shell_map[i]-1].add_orbital(fx3);
			molecule.atoms[shell_map[i]-1].add_orbital(fy3);
			molecule.atoms[shell_map[i]-1].add_orbital(fz3);
			molecule.atoms[shell_map[i]-1].add_orbital(fx2y);
			molecule.atoms[shell_map[i]-1].add_orbital(fx2z);
			molecule.atoms[shell_map[i]-1].add_orbital(fy2z);
			molecule.atoms[shell_map[i]-1].add_orbital(fy2x);
			molecule.atoms[shell_map[i]-1].add_orbital(fz2y);
			molecule.atoms[shell_map[i]-1].add_orbital(fz2x);
			molecule.atoms[shell_map[i]-1].add_orbital(fxyz);
		}
	}
	//----------------------------------------------------
	this->get_overlap_m();
	for ( unsigned i=0; i<molecule.orb_energies.size(); i++){ molecule.orb_energies[i] *= 27.2114; }
	for ( unsigned i=0; i<molecule.orb_energies_beta.size(); i++){ molecule.orb_energies_beta[i] *= 27.2114; }
	molecule.update();
	m_log->input_message("HOMO energy: \n\t");
	m_log->input_message(molecule.homo_energy);
	m_log->input_message("\nLUMO energy: \n\t");
	m_log->input_message(molecule.lumo_energy);
	m_log->input_message("\nAtomic orbitals: \n\t");
	m_log->input_message( int(molecule.num_of_ao) );
	m_log->input_message("\n");
	//molecule.print();
	molecule.bohr_to_ang();
	//molecule.check_ed();
}
/**********************************************************************/
void gaussian_files::get_overlap_m(){
	m_log->input_message("Starting to store the overlap 1e intragrals from the GAUSSIAN LOG file!\n");
	string log_name = change_extension(name_f,".log");
	
	if ( IF_file(log_name.c_str() ) ){
		if ( !check_file_ext(".log",log_name.c_str() ) ) {
			cout << "Warning! The file has wrong etension name!" << endl;
			m_log->write_warning("Warning! The file has wrong etension name!");
		}
		molecule.m_overlap.resize( ( molecule.MOnmb*(molecule.MOnmb+1) )/2 );
		
		unsigned int over_in, over_fin;
		Ibuffer Buffer(log_name.c_str(),true);
		for(int i=0;i<Buffer.nLines;i++){
			if ( Buffer.lines[i].IF_line("***",0,"Overlap",1,3) ) 
				over_in  = i;
			if ( Buffer.lines[i].IF_line("***",0,"Kinetic",1,4) ) 
				over_fin = i;
		}
	
		int col_n = 0;
		int row_n = 0;
		int col_c = 0;
	
		for(int i=over_in+1;i<over_fin;i++){
			if ( Buffer.lines[i].words.size() == 1 ) col_n = stoi(Buffer.lines[i].words[0]);
			else if( Buffer.lines[i].words.size() > 1 && Buffer.lines[i].words[1].size() < 6) col_n = stoi(Buffer.lines[i].words[0]);
			else{
				row_n = stoi(Buffer.lines[i].words[0]) -1;
				col_c = col_n -1;
				for(int j=1;j<Buffer.lines[i].line_len;j++){
					molecule.m_overlap[col_c + (row_n*(row_n+1))/2] = Buffer.lines[i].pop_double_f(1);
					col_c++;
				}
			}
		}
	}else{
		m_log->write_error("Log file not opening! Error in parse overlap matrix for gaussian file.");
	}
}
/******************************************************************************/
gaussian_files::~gaussian_files(){}
//================================================================================
//END OF FILE
//================================================================================