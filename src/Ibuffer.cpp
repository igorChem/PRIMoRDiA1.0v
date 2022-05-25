//Ibuffer.cpp

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
 
//------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
//------------------------------------------
#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/Iline.h"
#include "../include/Ibuffer.h"
//------------------------------------------
using std::move;
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

/*==============================================*/
Ibuffer::Ibuffer()	:
	nLines(0)		,
	name("buffer")	{
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name	,
						bool parse	)	:
	nLines(0)							,
	name(file_name)						{
	
	char tmp_line[500];
	if ( parse ){
		if ( IF_file(file_name) ){
			ifstream buf(file_name);
			name  = file_name;
			while( !buf.eof() ){
				buf.getline(tmp_line,500);
				lines.emplace_back( tmp_line );
				nLines++;
			}
			buf.close();
			parsed = true;
		}else{
			string message = "Not possible to open the file: ";
			message += name;
			message += "\n";
			cout << message << endl;
			m_log->input_message(message);
			parsed = false;
		}
	}else{
		if ( IF_file(file_name) ){
			m_log->input_message("Checking if the file can be pased and how many lines it has.");
			ifstream buf(file_name);
			name = file_name;
			while( !buf.eof() ){
				buf.getline(tmp_line,200);
				nLines++;
			}
			buf.close();
		}else{
			string message = "Not possible to open the file: ";
			message += name;
			message += "\n";
			m_log->input_message(message);
			parsed = false;
		}
		parsed = false;
	}
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name	,
				int in					,
				int fin)				:
				nLines(0)				,
				name(file_name)			{
	
	int in_indx  = in;
	int fin_indx = fin;
	char tmp_line[500];
	
	if ( IF_file(file_name) ){
		ifstream buf(file_name);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			if ( nLines > in_indx && nLines < fin_indx  ){
				lines.emplace_back( move (tmp_line) );
			}
			if ( nLines == fin_indx ) { break; }
				nLines++;
		}
		nLines = lines.size();
		buf.close();
		parsed = true;
	}else{
		string message	= "Not possible to open the file: ";
		message			+= name;
		message += "\n";
		m_log->input_message(message);
		parsed = false;
	}
	
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name	,
				 string wrdin			,
				 string wrdfin			):
	nLines(0)							,
	name(file_name)						{
		
	int in_indx  = -1;
	int fin_indx = 0;
	char tmp_line[500];
		
	if ( IF_file(file_name) ){
		std::ifstream buf(file_name);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			Iline Line(tmp_line);
			if ( in_indx == -1 ){
				if ( Line.words[0].compare(0,wrdin.size(),wrdin) == 0  ){
					in_indx = nLines;
					lines.emplace_back( Line );
				}
			}else if ( in_indx >= 0 ){
				lines.emplace_back( Line );
				if ( Line.IF_word( wrdfin,0,wrdfin.size() ) ) {
					break;
				}
			}
			nLines++;
		}
		buf.close();
		nLines = lines.size();
	}else{
		string message = "Not possible to open the file: ";
		message+= name;
		message += "\n";
		m_log->input_message(message);
		parsed = false;
	}
}
/******************************************************************/
Ibuffer::Ibuffer(const char* file_name		, 
				vector<string>& wrds_in		,
				vector<string>& wrds_fin )	:
	nLines(0)								,
	name(file_name)							{
	
	int in_indx  = -1;
	int fin_indx = 0;
	char tmp_line[500];
	
	if ( IF_file(file_name) ){
		ifstream buf(file_name);
		while( !buf.eof() ){
			buf.getline(tmp_line,500);
			Iline Line (tmp_line);
			if ( in_indx == -1 ){
				for(unsigned int i=0;i<wrds_in.size(); i++){
					if ( Line.IF_word( wrds_in[i],0,wrds_in[i].size() ) ) {
						if ( in_indx == -1 ) 
							in_indx = nLines;
					}
				}
			}
			if ( in_indx >=0 ){
				for(unsigned int j=0;j<wrds_fin.size();j++){
					if ( Line.IF_word( wrds_fin[j],0,wrds_fin[j].size() ) ) {
						if ( fin_indx == 0 ) 
							fin_indx = nLines;
						break;
					}
				}
			}
			nLines++;
		}
		if ( in_indx >=0 ){
			buf.close();
			ifstream buf2(file_name);
			nLines = 0;
			int real_nlines = 0;
			while( !buf2.eof() ){
				buf2.getline(tmp_line,500);
				if ( nLines > in_indx && nLines < fin_indx  ){
					lines.emplace_back( tmp_line );
					real_nlines++;
				}
				if ( nLines == fin_indx ) {	break; }
				nLines++;
			}
			nLines = real_nlines;
			buf2.close();
		}else{
			nLines = 0;
		}
	}else{
		string message = "Not possible to open the file: ";
		message += name;
		message += "\n";
		m_log->input_message(message);
		parsed = false;
	}
}
/******************************************************************/
Ibuffer& Ibuffer::get_block(int in, int fin){
	if ( parsed ){
		vector<Iline> temp;
		for(int i=in;i<fin;i++) { temp.emplace_back( move(lines[i] ) ); }
		lines = move(temp);
		return *this;
	}else{
		if ( IF_file(name) ){
			if ( in !=0 ){
				ifstream buf(name);
				nLines = 0;
				string tmp_line;
				int real_nlines = 0;
				while( !buf.eof() ){
					getline(buf,tmp_line);
					if ( nLines > in && nLines < fin  ){
						lines.emplace_back( move (tmp_line) );
						real_nlines++;
					}
					if ( nLines == fin ) {	break; }
					nLines++;
				}
				nLines = real_nlines;
				buf.close();
				parsed = true;
			}else{
				nLines = 0;
				parsed = false;
			}
		}else{
			string message = "Not possible to open the file: ";
			message += name;
			message += "\n";
			m_log->input_message(message);
			parsed = false;
		}
		return *this;
	}
}
/******************************************************************/
void Ibuffer::print(){
	cout << name   << endl;
	cout << nLines << endl;
}
/*******************************************************************/
void Ibuffer::clear(){
	lines.clear();
	nLines = 0;
}
/******************************************************************/
Ibuffer::~Ibuffer(){}
//================================================================================
//END OF FILE
//================================================================================