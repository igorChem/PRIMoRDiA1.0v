//interface.cpp

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
#include <memory>
#include <cstring>
#include <algorithm>
#include <string>
#include <fstream>
#include <ctime>
#include <experimental/filesystem>
#include <iomanip>

#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Imolecule.h"
#include "../include/QMparser.h"
#include "../include/Icube.h"
#include "../include/Ibuffer.h"
#include "../include/Iline.h"
#include "../include/gridgen.h"
#include "../include/primordia.h"
#include "../include/autoprimordia.h"
#include "../include/test_p.h"
#include "../include/interface.h"
#include "../include/pos_traj.h"
#include "../include/scripts.h"

/*********************************************************/
using std::unique_ptr;
using std::cout;
using std::endl;
using std::stoi; 
using std::move;
using std::string;
namespace fs = std::experimental::filesystem;
/***********************************************************************/
interface::interface()	:
	m_argc(0)			{
}
/***********************************************************************/
interface::interface(int argc, char** argv) :
	m_argc(argc)							{
	
	//---------------------------------
	// Filling argument vector
	for( int i =0; i<m_argc; i++ ){
		m_argv.emplace_back(argv[i]);
	}
	runtyp = m_argv[1];
	
	//---------------------------------
	// Setting time/date to print in the screen
	time_t now	= time(0);
	char* dt	= ctime(&now);
	cout << "Starting PRIMoRDiA software! "	<< endl;
	cout << "Calculations starting at: "	<< dt << endl;
	
	
	//---------------------------------
	// Setting general options
	for(int i=0;i<m_argc;i++){
		if      ( m_argv[i] == "-np")		NP			= stoi(m_argv[i+1]);
		else if ( m_argv[i] == "-verbose")	M_verbose 	= true;
	}
	
	//----------------------------------
	// Writting basic information in the log file
	m_log->initialize(M_verbose);
	m_log->input_message("Starting PRIMoRDiA software!\n ");
	m_log->input_message("Calculations Starting at: ");
	m_log->input_message(dt);
	m_log->input_message("\n");
	m_log->inp_delim(2);
	m_log->input_message("Provided Arguments: ");
	for ( unsigned i=0; i<m_argc; i++ ){
		m_log->input_message("\n\t");
		m_log->input_message(m_argv[i]);
	}
	m_log->input_message("\n");
	m_log->inp_delim(2);
}
/***********************************************************************/
void interface::run(){
	if ( runtyp == "-f" ){
		AutoPrimordia rds ( m_argv[2].c_str() );
		rds.init();
	}
	else if ( runtyp == "--help" || runtyp == "-h" ) this->write_help(); 
	else if ( runtyp == "-mo" )		this->MO_cube();
	else if ( runtyp == "-ed" )		this->ED_cube();
	else if ( runtyp == "-cp" )		this->Comp_cube();
	else if ( runtyp == "-input" )	this->write_input();
	else if ( runtyp == "-test" ) 	this->test_run();
	else if ( runtyp == "-cdiff" ){
		cube_diffs diffe(m_argv[2].c_str());
		diffe.write();
	}
	else if ( runtyp == "-cubed"){
		const char* file_name1 = m_argv[2].c_str();
		const char* file_name2 = m_argv[3].c_str();
		Icube cub1(file_name1);
		Icube cub2(file_name2);
		double diff = cub1.diff_integral(cub2);
		cout << "Absolute difference: "<< diff << endl;
		double similarity = cub1.similarity_index(cub2,"default");
		cout << "Similarity index: " << similarity << endl;
	}
	else if ( runtyp == "-int"){
		Icube cube(m_argv[2].c_str());
		cout << cube.calc_cube_integral() << endl;
	}	
	else if ( runtyp == "-p_traj" ){
		std::vector<int> rs_l;
		if ( m_argc > 2 ){
			for(int i=2;i<m_argc;i++){
				rs_l.push_back( stoi(m_argv[i] ) );
			}
			traj_rd traj(rs_l);
			traj.init_from_folder();
			traj.calculate_res_stats();
			traj.write_residues_reports();
			string from_f = "from_folder";
			scripts res_analy( from_f.c_str(), "residuos_analysis" );
			res_analy.write_r_residuos_barplot();
		}
	}
	else {
		cout << "No valid run option!" << endl;
		exit(-1);
	}	
}
/***********************************************************************/
void interface::MO_cube(){
	m_log->input_message("You are using PRIMoRDiA for generation of molecular orbital scalar field cube file!\n");
	
	unique_ptr<QMparser> qmfile ( new QMparser( m_argv[2].c_str(), m_argv[5]) );
	Imolecule molecule ( move(qmfile->get_molecule() ) );
	qmfile.reset(nullptr);
		
	unique_ptr<gridgen> main_grid( new gridgen( stoi(m_argv[3]), move(molecule) ) );
	if ( m_argv[5] == "orca"  ) 
		main_grid->calculate_orb_orca( stoi(m_argv[4]), false );
	else { 
		main_grid->calculate_orb( stoi(m_argv[4].c_str()), false );
	}
	main_grid->write_grid();
}
/***********************************************************************/
void interface::ED_cube(){
	m_log->input_message("You are using PRIMoRDiA for generation of total electron density scalar field cube file!\n");
	
	unique_ptr<QMparser> qmfile ( new QMparser( m_argv[2].c_str(),m_argv[4] ) );
	Imolecule molecule ( move(qmfile->get_molecule() ) );
	qmfile.reset(nullptr);
	
	unique_ptr<gridgen> main_grid( new gridgen( stoi(m_argv[3].c_str() ), move(molecule) ) );
	if ( m_argv[4].c_str() == "orca" )  
		main_grid->calculate_density_orca(); 
	else { 
		main_grid->calculate_density();
	}
	main_grid->write_grid();
}
/***********************************************************************/
void interface::Comp_cube(){		
	unique_ptr<QMparser> qmfile ( new QMparser ( m_argv[2].c_str(),m_argv[4]) );
	Imolecule molecule ( move(qmfile->get_molecule() ) );
	qmfile.reset(nullptr);
	unique_ptr<gridgen> dens ( new gridgen( stoi(m_argv[3].c_str() ),move(molecule) ) );
	dens->calculate_density();
	dens->write_grid();
	Icube comple  = dens->density.calculate_complement(false);
	string cpname = m_argv[2];
	cpname += "_complement.cube";
	comple.write_cube(cpname);	
}
/***********************************************************************/
void interface::write_help(){
		cout	<< "PRIMoRDiA help page\n"
				<< "PRIMoRDiA Macromolecular Reactivity Descriptor Acess\n"
				<< "The program must be run as follows:\n"
				<< "/path/to/executable [option run] [file_name] [other options] \n"
				<< "options run:\n"
				<< "--help/-h : Display this help message\n"
				<< "-f    : Reactivity descriptors run option\n"
				<< "-ed   : Electron density cube file generation run option\n"
				<< "-mo   : Molecular Orbital cube file generation run option\n"
				<< "-input: Produce input from the name list in the current folder\n"
				<< "-cp   : Electron density complement\n"
				<< "-cubed: cube file differences and similarity index calculation\n"
				<< "-cdiff: Calculates the similarity index from a list of cube files\n"
				<< "-int  : Calculates the integral of the cube file\n"
				<< "Generic options is the options must be placed after all the other arguments\n"
				<< "Generic options:\n"
				<< "-np [n] : program runs using n threads\n"
				<< "-log    : program produces a log file of its operations\n"
				<< "-verbose: program prints to the console messages about its operations\n"
				<< endl;
}
/***********************************************************************/
void interface::test_run(){
	test_p teste;
	//teste.init_general_test();
	//teste.test_reaction_analysis();
	//teste.test_traj_analysis();
}
/***********************************************************************/
void interface::write_input(){
	int option	= 0;
	int grid	= 0;
	double box[3];
	box[0] = 0;
	box[1] = 0;
	box[2] = 0;
	std::string lh			= "none";
	std::string program		= "none ";
	std::string mep			= "";
	std::string band_method	= "EW";
	std::string ext			= "";
	std::string RT			= "normal";
	
	int band				= 0;
	for( int i=0;i<m_argc;i++){
		if 		( m_argv[i] == "-op" )		option	= stoi(m_argv[i+1]);
		else if	( m_argv[i] == "-p" )		program	= m_argv[i+1];
		else if	( m_argv[i] == "-rt" )		RT		= m_argv[i+1];
		else if	( m_argv[i] == "-grid")		grid	= stoi(m_argv[i+1]); 
		else if	( m_argv[i] == "-lh")		lh		= m_argv[i+1]; 
		else if	( m_argv[i] == "-mep")		mep		= "mep";
		else if	( m_argv[i] == "-band")		band	= stoi(m_argv[i+1]);
		else if	( m_argv[i] == "-bandmethod") band_method = (m_argv[i+1]);
		else if	( m_argv[i] == "-box") {
			box[0] = stod(m_argv[i+1]); 
			box[1] = stod(m_argv[i+2]); 
			box[2] = stod(m_argv[i+3]); 
		}
	}
	
	std::vector<std::string> fnames; 
	fs::path c_path = fs::current_path();
	for (const auto & entry : fs::directory_iterator(c_path)){
		fnames.push_back( entry.path().filename() );
	}
	std::sort( fnames.begin(), fnames.end() );
		
	std::ofstream inp_file("primordia.input");
	inp_file << "#RT " << RT << " \n"
			 << "#PR eband 1\n";	
	
	for(int i=0;i<fnames.size();i++){
		if ( program == "mopac"){
			if ( check_file_ext(".aux",fnames[i].c_str()) ){
				if ( option == 1 )
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << program << " " << mep << endl;
				if ( option == 3 ){
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << band << " " << " "
								<< change_extension(fnames[i].c_str(),".pdb") << " " << program << " " << " "
								<< 0 << " " <<  0  << " " <<  0  << " " << 0 << " " << band_method << " " << mep << endl;
				}
			}else if( check_file_ext(".mgf",fnames[i].c_str()) ){
				if ( option == 2 ) {
					inp_file << option	<< " " << fnames[i] 
								<< " "	<< change_extension(fnames[i].c_str(),"_cat.mgf")	<< " "
								<< " "	<<	change_extension(fnames[i].c_str(),"_an.mgf")	<< " "
								<< lh	<<	" " << grid  << " " << 1 << " " << program << " " <<  mep << endl;
				}
			}		
			else if( check_file_ext(".out", fnames[i].c_str()) ){
				if ( option == 1 || option == 3 ){
					inp_file << option << " " << fnames[i] << " " << lh << " " << grid << " " << program << " " << mep << endl;
				}
			}			
		}
		else if ( program == "gamess"){
			if ( check_file_ext(".log",fnames[i].c_str()) ){
				if ( option == 1 )
					inp_file << option << " " << fnames[i] << " " << lh << " " << grid << " " << program << " " << mep << endl;
				if ( option == 2 ){
					inp_file << option	<< " " << fnames[i] 
								<< " "	<<	change_extension(fnames[i].c_str(),"_cat.log")	<< " "
								<< " "	<<	change_extension(fnames[i].c_str(),"_an.log")	<< " "
								<< lh	<< " " << grid 	<< " " << 1 << " " << program << " " << mep << endl;
				}
				if ( option == 3 ){
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << band << " " << " "
								<< change_extension(fnames[i].c_str(),".pdb") << " " << program << " " 
								<< 0 << " " <<  0  << " " <<  0  << " " << 0 << " " << band_method << " " << mep << endl;
				}
			}
		}
		else if ( program == "gaussian"){
			if ( check_file_ext(".fchk",fnames[i].c_str()) ){
				if ( option == 1 )
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << program << " " << mep << endl;
				if ( option == 2 ){
					inp_file << option	<< " " << fnames[i] 
								<< " "	<<	change_extension(fnames[i].c_str(),"_cat.fchk")	<< " "
								<< " "	<<	change_extension(fnames[i].c_str(),"_an.fchk")		<< " "
								<< lh	<<	" " << grid		<< " " << 1 << " " << program << " " << mep << endl;
				}
				if ( option == 3 ){
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << band << " " << " "
								<< change_extension(fnames[i].c_str(),".pdb") << " " << program << " " << " "
								<< 0 << " " <<  0  << " " <<  0  << " " << 0 << " " << band_method << " " << mep << endl;
				}
			}
		}
		else if ( program == "orca"){
			if ( check_file_ext(".out",fnames[i].c_str()) ){
				if ( option == 1 )
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << program << " " << mep << endl;
				if ( option == 2 ){
					inp_file << option	<< " " << fnames[i] 
								<< " "	<<	change_extension(fnames[i].c_str(),"_cat.out")	<< " "
								<< " "	<<	change_extension(fnames[i].c_str(),"_an.out")	<< " "
								<< lh	<<	" " << grid 	<< " " << 1 << " " << program << " " << mep << endl;
				}
				if ( option == 3 ){
					inp_file << option << " " << fnames[i] << " " <<	lh << " " << grid << " " << band << " " << " "
								<< change_extension(fnames[i].c_str(),".pdb") << " " << program << " " << " "
								<< 0 << " " <<  0  << " " <<  0  << " " << 0 << " " << band_method << " " <<mep << endl;
				}
			}
		}else{
				cout << "No program keyword recognized!" << endl;
		}
	}
	inp_file.close();
}
/***********************************************************************/
void interface::print_options(){
	for( int i=0;i<m_argc;i++){
		cout << m_argv[i] << endl;
	}
}
/***********************************************************************/
interface::~interface(){
	
}
//////////////////////////////////////////////////
//================================================================================
//END OF FILE
//================================================================================