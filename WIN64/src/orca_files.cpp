//orca_files.cpp

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

//Including c++ headers
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
#include "../include/orca_files.h"

using std::vector;
using std::string;
using std::stoi;
using std::stod;
using std::endl;
using std::cout;

//=======================================================
orca_files::orca_files()		:
	name_f("none")				,
	is_open(false)				{
}
/**********************************************************************************/
orca_files::orca_files(const char* file_name)	:
	is_open( false )							{
	
	if ( IF_file(file_name) ){
		is_open = true;
		if ( !check_file_ext(".out",file_name) ) {
			cout << "Warning! The file has wrong etension name!" << endl;
			m_log->write_warning("Warning! The file has wrong etension name!");
			is_open = false;
		}
	}else{
		m_log->write_error("Error opening ORCA file! Verify its presence in the current directory!");
		cout << "the file named " << file_name << " cannot be opened! " << endl;
	}
	name_f = file_name;
	molecule.name = get_file_name( name_f );
	molecule.name = remove_extension( molecule.name.c_str() );
}
/**********************************************************************************/
void orca_files::parse_out(){
	m_log->input_message("Starting to parse out file from ORCA.\n");
	
	int in_coords	= 0;
	int fin_coords	= 0;
	int orbs_in		= 0;
	int orbs_fin	= 0;
	int orbs_in_b	= 0;
	int orbs_fin_b	= 0;
	int chg_in		= 0;
	int chg_fin		= 0;
	int mo_in		= 0;
	int mo_fin		= 0;
	int ov_in		= 0;
	int ov_fin		= 0;
	
	vector<int> basis_in;
	vector<int> basis_fin;
	
	int col_n = 0;
	int row_n = 0;
	int col_c = 0;
	int line_indicator = 0;
	
	int aonum = 0;
	
	Ibuffer Buffer(name_f,true) ;
	for ( unsigned i=0; i<Buffer.nLines; i++){
		if		( Buffer.lines[i].IF_line("CARTESIAN",0,"(ANGSTROEM)",2,3) ){ in_coords = i; }
		else if ( Buffer.lines[i].IF_line("CARTESIAN",0,"(A.U.)",2,3) ){ fin_coords = i; }
		else if ( Buffer.lines[i].IF_line("Number",0,"Electrons",2,6) ){
			molecule.num_of_electrons = Buffer.lines[i].get_int(5);
		}
		else if ( Buffer.lines[i].IF_line("Total",0,"Energy",1,7) ){
			molecule.energy_tot = Buffer.lines[i].get_double(5);
		}
		else if ( Buffer.lines[i].IF_line("Basis",1,"set",2,7) ){
			basis_in.push_back(i);
		}
		else if( Buffer.lines[i].IF_line("end;",0,1) ){
			basis_fin.push_back(i);
		}
		else if( Buffer.lines[i].IF_line("contracted",2,"basis",3,7) ){
			aonum = Buffer.lines[i].get_int(6);
		}
		else if ( Buffer.lines[i].IF_line("ORBITAL",0,"ENERGIES",1,2) ) { 
			orbs_in = i; 
		}
		else if ( Buffer.lines[i].IF_line("UP",1,"ORBITALS",2,3) ) { 
			if ( orbs_in ==  0 ) {
				orbs_in = i;
			}
		}
		else if ( Buffer.lines[i].IF_line("DOWN",1,"ORBITALS",2,3) ) { 
			if ( orbs_in_b ==  0 ) {
				orbs_fin = i;
				orbs_in_b = i;
				molecule.betad = true;
			}
		}
		else if ( Buffer.lines[i].IF_line("OVERLAP",0,"MATRIX",1,2) ){
			ov_in = i;
		}
		else if ( Buffer.lines[i].IF_line("INITIAL",0,"MOREAD",2,3) ){
			ov_fin = i;
		}
		else if ( Buffer.lines[i].IF_line("DFT",0,"GENERATION",2,3) ){
			if ( ov_fin == 0 ) 	ov_fin = i;
		}
		else if ( Buffer.lines[i].IF_line("MOLECULAR",0,"ORBITALS",1,2) ){
			mo_in = i;
			if ( orbs_fin == 0 ) 
				orbs_fin = i;
			if ( orbs_in_b >  0 ) 
				orbs_fin_b = i;
		}
		else if ( Buffer.lines[i].IF_line("MULLIKEN",1,"ANALYSIS",3,5) ) { 
			if ( mo_fin == 0 ){
				mo_fin = i; 
			}
		}
		else if ( Buffer.lines[i].IF_line("MULLIKEN",0,"CHARGES",2,3) ) { chg_in = i; }
		else if ( Buffer.lines[i].IF_line("MULLIKEN",0,"CHARGES",2,6) ) { chg_in = i; }
		else if ( Buffer.lines[i].IF_line("MULLIKEN",0,"REDUCED",1,4) ) {  chg_fin = i; }
		else if ( Buffer.lines[i].IF_line("MULLIKEN",0,"REDUCED",1,7) ) {  chg_fin = i; }
	}
	for(unsigned j=in_coords; j<fin_coords; j++ ){
		if ( Buffer.lines[j].line_len == 4 ){
			double xx, yy, zz;
			xx = Buffer.lines[j].get_double(1);
			yy = Buffer.lines[j].get_double(2);
			zz = Buffer.lines[j].get_double(3);
			string type_= Buffer.lines[j].get_string(0);
			molecule.add_atom(xx,yy,zz,type_);
		}
	}
	vector<basis_orca> basisset;
	for( unsigned i=0; i<basis_in.size(); i++ ){
		basis_orca bo;
		for( unsigned j=basis_in[i]; j<basis_fin[i]; j++ ){
			if ( Buffer.lines[j].line_len == 7 ){
				bo.element_type = Buffer.lines[j].get_string(6);
			}
			else if ( Buffer.lines[j].line_len == 2 ){
					if ( Buffer.lines[j].words[0] != "NewGTO" ){
						bo.shell_sym.push_back(Buffer.lines[j].get_string(0) );
						bo.shell_size.push_back(Buffer.lines[j].get_int(1) );
					}
			}
			else if ( Buffer.lines[j].line_len == 3 ){
				bo.coefficients.push_back(Buffer.lines[j].get_double(1) );
				bo.exp.push_back(Buffer.lines[j].get_double(2) );
			}
		}
		basisset.push_back(bo);
	}

	for( unsigned i=(orbs_in +1); i<orbs_fin; i++ ){
		if ( Buffer.lines[i].line_len == 4  ) {
			if ( Buffer.lines[i].words[0] != "NO" ){
				molecule.orb_energies.push_back( Buffer.lines[i].get_double(3) ); 
				molecule.occupied.push_back( Buffer.lines[i].get_double(1) );
				molecule.MOnmb++;
			}
		}
	}
	
	for( unsigned i=(orbs_in_b +1); i<orbs_fin_b; i++ ){
		if ( Buffer.lines[i].line_len == 4 ) {
			if ( Buffer.lines[i].words[0] != "NO" ){
				molecule.orb_energies_beta.push_back( Buffer.lines[i].get_double(3) ); 
				molecule.occupied_beta.push_back( Buffer.lines[i].get_double(1) );
				molecule.MOnmb_beta++;
			}
		}
	}
	
	aonum = molecule.MOnmb;
	
	col_n = 0;
	row_n = 0;
	col_c = 0;
	line_indicator = 0;
	bool bet = false;
	
	molecule.coeff_MO.resize(aonum*aonum);
	if ( molecule.betad ) { 
		molecule.coeff_MO_beta.resize(aonum*aonum);
	}
	
	Buffer.clear();
	Ibuffer Buffer2(name_f,mo_in,mo_fin);
	unsigned l = 0;
	for( unsigned j=1; j<Buffer2.nLines; j++){
		if ( Buffer2.lines[j].line_len > 0 && line_indicator == 0 ) {
			col_n = Buffer2.lines[j].pop_int(0);
			line_indicator++;
			row_n = 0;
		}
		else if ( Buffer2.lines[j].line_len > 0 && line_indicator == 1  ) { line_indicator++; }
		else if ( Buffer2.lines[j].line_len > 0 && line_indicator == 2  ) { line_indicator++; }
		else if ( Buffer2.lines[j].line_len > 0 && line_indicator == 3  ) { line_indicator++; }
		else if ( Buffer2.lines[j].line_len >=3 && line_indicator == 4 && bet ){
			for( l=0; l<Buffer2.lines[j].line_len-2; l++ ){ 
				molecule.coeff_MO_beta[(col_n+l)*aonum+row_n] = Buffer2.lines[j].pop_double(2);
			}
			row_n++;
			if ( row_n == aonum ) line_indicator = 0;
			if ( col_n+l == aonum   ){
					line_indicator++;
			}
		}
		else if ( Buffer2.lines[j].line_len >= 3 && line_indicator == 4 ) {
			for( l=0; l<Buffer2.lines[j].line_len-2; l++ ){ 
				molecule.coeff_MO[(col_n+l)*aonum+row_n] = Buffer2.lines[j].pop_double(2);
			}
			row_n++;
			if ( row_n == aonum ) {
				line_indicator = 0;
				if ( col_n+l == aonum && molecule.betad ){
					bet = true;
				}else if ( col_n+l == aonum   ){
					line_indicator++;
				}
			}
		}
	}
	
	int counter = 0;
	int ls = 0; //line number of words 
	Buffer2.clear();
	Ibuffer Buffer3(name_f,chg_in,chg_fin);
	for( unsigned i=0; i<Buffer3.nLines; i++ ){
		if  ( molecule.atoms[counter].element.size() > 1 ){
			ls = 3;
		}else {
			ls = 4;
		}
		if ( !molecule.betad ) {
			if ( Buffer3.lines[i].line_len == ls && counter < molecule.num_of_atoms ) {
				molecule.atoms[counter++].charge = Buffer3.lines[i].get_double(ls-1);
			}
		}else{
			ls++;
			if ( Buffer3.lines[i].line_len == ls && counter < molecule.num_of_atoms ) {
				molecule.atoms[counter++].charge = Buffer3.lines[i].get_double(ls-2);
			}
		}
	}

		
	for ( unsigned i=0; i<molecule.num_of_atoms; i++ ){
		for ( unsigned j=0; j<basisset.size(); j++ ){
			if ( molecule.atoms[i].element == basisset[j].element_type ){
				int cnt =0;
				for(unsigned k=0; k<basisset[j].shell_size.size(); k++ ){
					if ( basisset[j].shell_sym[k] == "S" ){
						Iaorbital orb;
						orb.spherical = true;
						for ( unsigned l=0; l<basisset[j].shell_size[k]; l++ ){
							orb.add_primitive(basisset[j].coefficients[cnt],basisset[j].exp[cnt]);
							cnt++;
						}
						molecule.atoms[i].add_orbital(orb);
					}
					else if ( basisset[j].shell_sym[k] == "P" ){
						Iaorbital orbpx, orbpy, orbpz;
						orbpx.spherical = true;
						for ( int l=0;l<basisset[j].shell_size[k];l++){
							orbpx.add_primitive(basisset[j].coefficients[cnt],basisset[j].exp[cnt]);
							cnt++;
						}
						orbpy = orbpz 	= orbpx; 
						orbpx.symmetry 	= "PX";
						orbpy.symmetry 	= "PY";
						orbpz.symmetry 	= "PZ";
						orbpx.powx		= 1;
						orbpy.powy      = 1;
						orbpz.powz         = 1;
						molecule.atoms[i].add_orbital(orbpx);
						molecule.atoms[i].add_orbital(orbpy);
						molecule.atoms[i].add_orbital(orbpz);
					}
					else if ( basisset[j].shell_sym[k] == "D" ){
						Iaorbital orbdz2, orbdxz, orbdyz,orbdx2y2,orbdxy;
						orbdz2.spherical = true;
						for ( int l=0;l<basisset[j].shell_size[k];l++){
							orbdz2.add_primitive(basisset[j].coefficients[cnt],basisset[j].exp[cnt]);
							cnt++;
						}
						orbdxy = orbdx2y2  = orbdxz  = orbdxy = orbdz2 ; 
						orbdz2.symmetry = "D0";
						orbdxz.symmetry = "D1p";
						orbdyz.symmetry = "D1n";
						orbdx2y2.symmetry = "D2n";
						orbdxy.symmetry = "D2p";
						molecule.atoms[i].add_orbital(orbdz2);
						molecule.atoms[i].add_orbital(orbdxz);
						molecule.atoms[i].add_orbital(orbdyz);
						molecule.atoms[i].add_orbital(orbdx2y2);
						molecule.atoms[i].add_orbital(orbdxy);
					}
					else if ( basisset[j].shell_sym[k] == "F" ){
						Iaorbital orbf0, orbf1p, orbf1n,orbf2p,orbf2n,orbf3p,orbf3n;
						orbf0.spherical = true;
						for ( int l=0; l<basisset[j].shell_size[k]; l++ ){
							orbf0.add_primitive(basisset[j].coefficients[cnt],basisset[j].exp[cnt]);
							cnt++;
						}
						orbf3n = orbf3p = orbf2n = orbf2p = orbf1n = orbf1p = orbf0;
						orbf0.symmetry = "f0";
						orbf1p.symmetry = "f1p";
						orbf1n.symmetry = "f1n";
						orbf2p.symmetry = "f2p";
						orbf2n.symmetry = "f2n";
						orbf3p.symmetry = "f3p";
						orbf3n.symmetry = "f3n";
						molecule.atoms[i].add_orbital(orbf0);
						molecule.atoms[i].add_orbital(orbf1p);
						molecule.atoms[i].add_orbital(orbf1n);
						molecule.atoms[i].add_orbital(orbf2p);
						molecule.atoms[i].add_orbital(orbf2n);
						molecule.atoms[i].add_orbital(orbf3p);
						molecule.atoms[i].add_orbital(orbf3n);
					}
					else if ( basisset[j].shell_sym[k] == "G" ){
						Iaorbital G_orb;
						m_log->write_warning("Reading ortbital G, their contribution for volumetric data are not accounted by the software!");
						for( unsigned ll=0; ll<9; ll++ ){
							molecule.atoms[i].orbitals.push_back(G_orb);
						}
					}						
				}
			}
		}
	}
	if ( basisset.size() == 0 ) {
		molecule.num_of_ao = aonum; // ao num read from file rather than which one was counted
	}
	this->get_overlap(ov_in,ov_fin);
	molecule.update();
	//molecule.print();
	molecule.check();
	m_log->input_message("Orca file parsed!");
	
	m_log->input_message("Total Energy \n\t");
	m_log->input_message( std::to_string(molecule.energy_tot) );
	m_log->input_message( "\nHOMO number: \n\t");
	m_log->input_message( int (molecule.homoN) );
	m_log->input_message("\nHOMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.homo_energy) );
	m_log->input_message( "\nLUMO number: \n\t");
	m_log->input_message( int (molecule.lumoN) );
	m_log->input_message("\nLUMO energy: \n\t");
	m_log->input_message( std::to_string(molecule.lumo_energy) );
	m_log->input_message("\nAtomic orbitals: \n\t");
	m_log->input_message( aonum );
	m_log->input_message("\n");

}
/****************************************************************/
void orca_files::get_overlap(int ov_in,int ov_fin){
	unsigned int col_n,line_indicator,col_c,row_n;
	col_c = col_n = line_indicator = row_n = 0;
	unsigned int aonum = molecule.get_ao_number();
	
	Ibuffer Buffer(name_f,ov_in,ov_fin);
	
	vector<double> overlap_full(aonum*aonum);
	for( unsigned i=1; i<Buffer.nLines; i++ ){
		if ( Buffer.lines[i].line_len>1 && line_indicator == 0 ){
			col_n = Buffer.lines[i].get_int(0);
			line_indicator++;
			//cout << col_n << " " << line_indicator << endl;
		}
		else if ( Buffer.lines[i].line_len > 1 && line_indicator == 1 ){
			col_c = col_n;
			for ( unsigned j=0; j<Buffer.lines[i].line_len-1; j++){
				overlap_full[col_c*aonum + row_n] = Buffer.lines[i].pop_double(1);
				col_c++;
			}
			row_n++;
			if ( row_n == aonum ){
				line_indicator = 0;
				row_n = 0;
				if ( col_c == aonum-1 ){
					line_indicator++;
				}
			}
		}
	}
	
	for ( unsigned i=0; i<aonum; i++){
		for ( unsigned j=0; j<=i; j++){
			molecule.m_overlap.push_back( overlap_full[i*aonum+j] );
		}
	}
}
/*****************************************************************/
orca_files::~orca_files(){}
//================================================================================
//END OF FILE
//================================================================================