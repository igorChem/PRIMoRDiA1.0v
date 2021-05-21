//Ibuffer.h

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
 

#ifndef IBUFFER
#define IBUFFER
//------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
//------------------------------------------
#include "../include/common.h"
//------------------------------------------
class Iline; //foward declaration
//-------------------------------------------------------------------------------------------
/**
 * Class to hold and manipulate text files. Creates Ilines objects and sotre in a STL vector 
 * to hold each line of the file to be easily acessed and parsed.
 * @class Ibuffer
 * @author Igor Barden Grillo
 * @date 07/04/18
 * @file Ibuffer.h
 * @brief Class to open files and intatiate and hold Iline objects for each line in the file.
 */
class Ibuffer {
	public:
		const char* name; // file name to be open.
		unsigned int nLines; // number of lines in the file.
		bool parsed; // if the lines of the file were parsed and stored.
		std::vector<Iline> lines; //Iline objects to hold each file line.
		Ibuffer(); // default constructor 
		Ibuffer(const char* file_name,bool parse); // constructor from file path to be parsed or not.
		Ibuffer(const char* file_name,int in, int fin); // constructor to store info from file path from block of lines.
		Ibuffer(const char* file_name,std::string wrdin,std::string wrdfin); // constructor to store info from file path from block of lines.
		Ibuffer(const char* file_name, std::vector<std::string>& wrds_in,std::vector<std::string>& wrds_fin);// constructor to store info from file path from block of lines.
		Ibuffer(const Ibuffer& rhs_ibuf) = delete; // copy constructor deleted
		Ibuffer& operator=(const Ibuffer& rhs_ibuf) = delete;// assign operator overload deleted.
		Ibuffer& get_block(int in, int fin); // Retunrn a class object as a block of lines from the current one.
		void clear();
		~Ibuffer(); // Destructor.
		void print(); // print to the console informations about the parsed file. 
};

#endif

//================================================================================
//END OF FILE
//================================================================================