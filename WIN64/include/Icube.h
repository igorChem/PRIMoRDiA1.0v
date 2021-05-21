// Icube.h 
// header file for the declarations ans documentation of the Icube class 

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

#ifndef CUBE
#define CUBE

#include <string>
#include <vector> 

#include "../include/common.h"
#include "../include/Imolecule.h"

//===================================================================================
/**
 * This class is meant to represent an cube file object in the gaussian format that are accepted by 
 * the most of quantum chemical GUI programs. The Icube object holds a Imolecule object, thus molecular information,
 * with name and cooridnates of the atoms, molecular orbital number the the Icube object may represent, and an array 
 * with the scalar field data.
 * @class Icube
 * @author Igor Barden Grillo
 * @date 20/03/18
 * @file Icube.h
 * @brief A class to represent volumetric chemical data that are placed in grid files.
 * The present class methods works with these type of files, parsing and rewriting them.
 * @see Imolecule
 */ 
class Icube {
	public:
	    // member variables
		std::string name;
		bool elec_dens;
		unsigned int MOn;
		unsigned int voxelN;
		std::string header;
		double origin[3];
		double gridsides[3];
		unsigned int grid[3];
		Imolecule molecule;
		std::vector<double> scalar;
		// constructors/destructor
		Icube();
		Icube(const char* file_nam);
		Icube(const Icube& rhs_cube); 
		Icube& operator=(const Icube& rhs_cube);  
		Icube(Icube&& rhs_cube) noexcept;
		Icube& operator=(Icube&& rhs_cube) noexcept; 
		~Icube();
		//friend functions
		friend bool operator==(const Icube& lhs_cube,const Icube& rhs_cube);
		friend Icube operator-(const Icube& lhs_cube,const Icube& rhs_cube);
		friend Icube operator+(const Icube& lhs_cube,const Icube& rhs_cube);
		friend Icube operator*(const Icube& lhs_cube,const Icube& rhs_cube);
		friend Icube operator/(const Icube& lhs_cube,const Icube& rhs_cube);
		friend Icube operator-(const Icube& lhs_cube,double value);
		friend Icube operator+(const Icube& lhs_cube,double value);
		friend Icube operator*(const Icube& lhs_cube,double value); 
		friend Icube operator/(const Icube& lhs_cube,double value);
		//member functions
		Icube scale_cube(double val);
		Icube log_cube();
		Icube SQ();
		double calc_cube_integral();
		Icube normalize();
		double diff_integral(const Icube& cube);
		double similarity_index(Icube& cube, std::string type);
		Icube calculate_complement(bool clog);
		double get_scalar(int x, int y, int z);
		double get_scalar(double x, double y, double z);	 
		void add_data(std::vector < std::vector < std::vector<double> > >& data);
		void write_cube(std::string cubeName);
		void get_cube_stats(double& mean, double& sum, double& min, double& max);
		void print();
};

//----------------------------------------------------------------------------------------------------------------------------------------------
/**
 * @class cube_diffs
 * @author Igor Barden Grillo
 * @date 22/05/19
 * @file Icube.h
 * @brief Class to calculate the similarity index between cube files given by a inpu files with the names of cubes.
 */
class cube_diffs {
	public:
		std::vector<double> diffs;
		std::vector<std::string> labels;
		cube_diffs(const char* name);
		cube_diffs(const cube_diffs& rhs) = delete;
		cube_diffs& operator=(const cube_diffs& rhs) = delete;
		~cube_diffs();
		void write();
};

#endif 
//================================================================================
//END OF FILE
//================================================================================