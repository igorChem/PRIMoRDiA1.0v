//log.cpp
//source code for Ilog class.

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
 
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <experimental/filesystem>
//----------------------------------------------------
#include "../include/common.h"
#include "../include/log_class.h"

//-----------------------------------------------------
using std::cout;
using std::endl;
using std::string;
namespace fs = std::experimental::filesystem;

string delim_1 = "======================================================\n";
string delim_2 = "------------------------------------------------------\n";
/************************************************************/
Ilog::Ilog()			:
	screen_output(false){
}
/************************************************************/
void Ilog::initialize(bool sout){
	screen_output	= sout;
	
	time_t now	= time(0);
	char* dt	= ctime(&now);

	log_file.open("primordia.log",std::ios::out | std::ios_base::app);
	
}
/************************************************************/
void Ilog::input_message(std::string message){
	log_file << message;
	if ( screen_output ) cout << message << endl;
}
/************************************************************/
void Ilog::input_message(int message){
	if ( screen_output ) cout << message << endl; 
	log_file << message << endl;
}
/************************************************************/
void Ilog::input_message(double message){
	if ( screen_output ) cout << message << endl; 
	log_file << message << endl;
}
/************************************************************/
void Ilog::write_error( std::string _Error ){
	log_file << endl;
	log_file << "###************************ERROR************************###\n";
	log_file << "\t";
	log_file << _Error << endl;	
	log_file << "#*******************End of Error Message*******************\n\n";
}
/************************************************************/
void Ilog::write_warning( std::string _Warning ){
	log_file << endl;
	log_file << "##======================WARNING======================##\n";
	log_file << "\t";
	log_file << _Warning << endl;
	log_file << "#===============End of Warning Message===============##\n\n";
}
/************************************************************/
void Ilog::inp_delim( int i ){
	if ( i == 1 ){
		log_file << delim_1 << endl;
	}else if ( i == 2 ){
		log_file << delim_2 << endl;
	}
}
/************************************************************/
void Ilog::timer(){
	if ( screen_output ) {
		cout << chronometer.return_wall_time() << " seconds" << endl;
	}
	log_file << chronometer.return_wall_time() << " seconds" << endl;
}
/************************************************************/
void Ilog::abort(std::string message){
	cout << message << endl;
	exit(-1);
}
/************************************************************/
Ilog::~Ilog(){
	log_file << "Exiting PRIMoRDiA after " << chronometer.return_wall_time() << " seconds" << endl;
	log_file << "===========================================================================" << endl;
	log_file.close();
}
/************************************************************/
