// header file for declarations of Ilile class
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
#ifndef ILINE
#define ILINE
//================================
#include <iostream>
#include <string>
#include <vector>

//=====================================================================================================
/**
 * Class to hold information and manipulate text lines. Reunite tools to parse files by make tests with
 * the text lines content and features.
 * @class Iline
 * @author Igor Barden Grillo
 * @date 07/04/18
 * @file Iline.h
 * @brief class to represent and manipulate a line from text.
 */
class Iline {
	public:
		std::vector<std::string> words; //words of the line.
		unsigned int line_len; //line length
		Iline();  // default constructor
		Iline(std::string Line); // constructor overload from string
		Iline(char* Line); // constructor overload from char*
		Iline(const Iline& rhs_line); // copy constructor
		Iline& operator=(const Iline& rhs_line);  // assign operator overload
		Iline(Iline&& rhs_line) noexcept; // move constructor
		Iline& operator=(Iline&& rhs_line) noexcept; // move assign operator overload
		~Iline(); // destructor
		
		bool IF_word(std::string s,int pos,int fin);
		bool IF_line(std::string s);
		bool IF_line(std::string s, int pos,int len);
		bool IF_line(std::string s1, int pos1,std::string s2, int pos2, int len);
		std::string get_line(); 
		int pop_int(int pos);
		int get_int(int pos);
		double pop_double(int pos);
		double get_double(int pos);
		double pop_double_f(int pos);
		double get_double_f(int pos);
		std::string pop_string(int pos);
		std::string& get_string(int pos);
		void print();
};
//=====================================================
#endif
//================================================================================
//END OF FILE
//================================================================================