//Iaorbital.h

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

#ifndef IAORBITAL
#define IAORBITAL

//inclding c++ headers
#include <string>
#include <vector>
//#include <cstring>
//including PRIMoRDiA headers
#include "../include/common.h"

//-------------------------------------------------------------------------
/**
 * This class is meant to hold objects with gaussian primitive atomic orbitals function
 * information, and manipulate and change this information, like calculating the normalization constants.
 * @class Iprimitive
 * @author Igor Barden Grillo
 * @date 11/09/18
 * @file Iaorbital.h
 * @brief Gaussian type orbital primitive for the contracted GTO atomic orbital class.
 */
class Iprimitive {
	public:
		double n_fact; // Normalization coefficient.
		double exponent; // primitive gaussian exponent value .
		double c_coef; // primitive gaussian contraction coefficient.
		//----------------------------------------
		Iprimitive(); // default constructor
		Iprimitive(double exp,double contrac_coeff); // constructor from data.
		Iprimitive(const Iprimitive& prim_rhs); // copy constructor.
		Iprimitive(Iprimitive&& prim_rhs) noexcept; // move constructor.
		Iprimitive& operator=(const Iprimitive& prim_rhs); // assignment operator overload.
		Iprimitive& operator=(Iprimitive&& prim_rhs) noexcept; // move assignment operator overload.
		~Iprimitive(); // destructor
		//------------------------------------------
		void print(); // function to print information in the console.
	
};
//========================================================================
/** 
 * This class is meant to be hold atom object abstraction to hold, manipulate and modify
 * information extracted from quantum chemical output programs.
 * @class Iaorbital
 * @author Igor Barden Grillo.
 * @date 20/03/18
 * @file Iatom.h
 * @brief  A Iaorbital Class declaration to instatiate atomic orbital object to Quantum Mechanical 
 * properties calculations.  
 *   
 */ 
class Iaorbital {
	public:
		unsigned int shell; // shell orbital level.
		bool gto; // if the orbital is gaussian type.
		bool spherical; // if the orbital has spherical symmetry.
		std::string symmetry; //
		double n_factor; // normalization factor.
		double alpha; // atomic orbital exponent 
		unsigned int powx; //  
		unsigned int powy; //
		unsigned int powz; //
		std::vector<Iprimitive> gtos; // primitive array.
		//----------------------------------------------
		Iaorbital(); // default constructor. 
		Iaorbital(unsigned int level, std::string sym, double coef); // constructor from data.
		Iaorbital(const Iaorbital& rhs_orb); // copy constructor
		Iaorbital& operator=(const Iaorbital& rhs_orb); // assignment operator overload.
		Iaorbital(Iaorbital&& rhs_orb) noexcept; //move constructor.
		Iaorbital& operator=(Iaorbital&& rhs_orb) noexcept; //move assignment operator overload.
		~Iaorbital(); // destructor
		//-----------------------------------------------
		void add_primitive(double expo, double contrac); // create a Iprimitive object and store in gtos vector. 
		void add_primitive(Iprimitive gorb); // create a Iprimitive object and store in gtos vector. 
		bool normalize(); // calculate the normalization factor.
		void print(); // function to print information in the console.
		

};

#endif
//================================================================================
//END OF FILE
//================================================================================