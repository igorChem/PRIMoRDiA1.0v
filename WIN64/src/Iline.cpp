// source file for Iline class 
//--------------------------------------
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
#include <sstream>
#include <vector>
//--------------------------------------
#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/Iline.h"
//--------------------------------------
using std::string;
using std::move;
using std::stod;
using std::stoi;
//=====================================================
Iline::Iline()			:
	line_len(0)			{
	words.reserve(50);
}
/*****************************************************************************/
Iline::Iline(string Line):
	line_len(0)				{
		
	words.reserve(30);
	std::stringstream stream(Line);
	string word;
	word.reserve(20);
	
	while (stream >> word){
		words.emplace_back( word );
		line_len++;
	}
}
/*****************************************************************************/
Iline::Iline(char* Line)	:
	line_len(0)				{
		
	words.reserve(30);
	std::stringstream stream(Line);
	string word;
	word.reserve(20);
	
	while (stream >> word){
		words.emplace_back( word );
		line_len++;
	}
}
/*****************************************************************************/
Iline::Iline(const Iline& rhs_line)	:
	words(rhs_line.words)			,
	line_len(rhs_line.line_len)		{
}
/*****************************************************************************/
Iline& Iline::operator=(const Iline& rhs_line){
	if( this != &rhs_line ){
		words	= rhs_line.words;
		line_len= rhs_line.line_len;
	}
	return *this;
}
/*****************************************************************************/
Iline::Iline(Iline&& rhs_line) noexcept	:
	words( move(rhs_line.words) )		,
	line_len( rhs_line.line_len)		{
}
/*****************************************************************************/
Iline& Iline::operator=(Iline&& rhs_line) noexcept {
	if( this != &rhs_line ){
		words	= move(rhs_line.words);
		line_len	= rhs_line.line_len;
	}
	return *this;
}
/*****************************************************************************/
bool Iline::IF_word(string& s, unsigned  pos){
	if ( line_len == 0 || pos > line_len-1 ){
		return false;
	}
	return words[pos].compare(s) == 0;	
}
/*****************************************************************************/
bool Iline::IF_word(string& s, unsigned pos, unsigned fin){
	if ( line_len == 0 || pos > line_len-1 ){ 
		return false;
	}
	if ( words[pos].size()<fin ) {
		return false;
	}
	return words[pos].compare(0,fin,s) == 0;
}
/*****************************************************************************/
bool Iline::IF_line(std::string s){
	if ( line_len == 0 ) return false;
	if ( this->get_line() == s ) return true;
	else return false;
}
/*****************************************************************************/
bool Iline::IF_line(std::string s, int pos,int len){
	if ( line_len == 0 || line_len != len ) return false;
	if ( words[pos] == s ) return true;
	else return false;
}
/*****************************************************************************/
bool Iline::IF_line(string s1, int pos1,string  s2, int pos2, int len){
	if ( line_len == 0 || line_len != len ) return false;
	if ( words[pos1] == s1 && words[pos2] == s2 ) return true;
	else return false;
}
/*****************************************************************************/
string Iline::get_line(){
	string line = "";
	for(unsigned int i=0;i<words.size();i++){ line += words[i] + " "; }
	return line;
}
/*****************************************************************************/
int Iline::pop_int(int pos){
	int res = stoi(words[pos]);
	words.erase( words.begin()+pos );
	return res;
}
/*****************************************************************************/
int Iline::get_int(int pos){ 
	int res = 0;
	try{
		res = stoi(words[pos]);
	}catch( const std::invalid_argument& ){
		std::cout << "Impossible to convert to an integer!" << std::endl;
		m_log->write_error("In convert some string to int!\n verify you input file, error may be in the position of an argument!\n");
		m_log->input_message("The problematic line is: \n\t");
		m_log->input_message( this->get_line() );
	}
	return res;
}
/*****************************************************************************/
double Iline::pop_double(int pos){
	double res = stod(words[pos] );
	words.erase( words.begin()+pos );
	return res;
}
/*****************************************************************************/
double Iline::get_double(int pos){
	double res = 0.000;
	try{
		res = stod(words[pos]);
	}catch( const std::invalid_argument& ){
		std::cout << "Impossible to convert to double!" << std::endl;
		m_log->write_error("In convert some string to int!\n verify you input file, error may be in the position of an argument!\n");
		m_log->input_message("The problematic line is: \n\t");
		m_log->input_message( this->get_line() );
	}
	return res;
} 
/*****************************************************************************/
double Iline::pop_double_f(int pos){
	double res = D_E_conv( words[pos] );
	words.erase( words.begin()+pos );
	return res;
}
/*****************************************************************************/
double Iline::get_double_f(int pos){ return D_E_conv( words[pos] ); }
/*****************************************************************************/
string Iline::pop_string(int pos){
	string res = words[pos];
	words.erase( words.begin()+pos );
	return res;
}
/*****************************************************************************/
string& Iline::get_string(int pos){ return words[pos]; }
/*****************************************************************************/
void Iline::print(){
	for(unsigned int i=0; i<words.size();i++){
		std::cout << "The word # " << i << " is : " << words[i] << std::endl;
	}
}
/*****************************************************************************/
Iline::~Iline(){}
/*****************************************************************************/
//================================================================================
//END OF FILE
//================================================================================