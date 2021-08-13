//claas header for grigen.cpp to generating electron density cubes from MO 
// grigen.h


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

#ifndef GRIDGEN_H
#define GRIDGEN_H
// C++ header files
#include <iostream>
#include <string>
//Our library includes
#include "../include/common.h"
#include "../include/Icube.h"

class Iaorbital;
class Imolecule;
//======================================================================================
/**
 * Class to automate the generation of cube file and Icube objects from the generation of
 * three dimensional scalar fields. 
 * @class gridgen
 * @author barden
 * @date 20/03/18
 * @file gridgen.h
 * @brief Class to generate grid scalar values from quantum computational calculations.
 */
class gridgen {
	public:
		//member variables
		std::string name;
		bool orbital;
		int Norb;
		unsigned int points;
		Imolecule molecule;
		double origin[3];
		double top_corner[3];
		double grid_sides[3];
		unsigned int grid_len[3];
		std::vector<double> AOxcoords;
		std::vector<double> AOycoords;
		std::vector<double> AOzcoords; 
		std::vector<Iaorbital> orbs;
		std::vector< std::vector < std::vector<double> > > psi;
		Icube density;
		// constructos/destructor
		gridgen();
		gridgen(int grd, Imolecule&& mol) noexcept;
		gridgen(const gridgen& rhs_grd) = delete;
		gridgen& operator=(const gridgen& rhs_grd) = delete;
		~gridgen();
		//member functions
		double calc_slater_orb(int i, int x, int y, int z);
		double calc_gauss_orb(int i, int x, int y, int z);
		double calc_orca_sphe(int i, int x, int y, int z);
		double calc_aorb(int i, int x, int y, int z);
		double calc_orb_voxel(int nm,int x,int y,int z, bool beta);
		double calc_orb_voxel_orca(int nm,int x,int y,int z,bool beta);
		void calculate_orb(int Nmo,bool beta);
		void calculate_orb_orca(int Nmo,bool beta);
		double electron_density_mo(int x, int y, int z);
		double electron_density_mo_orca(int x, int y, int z);
		double electron_density(int x, int y, int z);
		void calculate_density();
		void calculate_density_orca();
		void orbs_overlap();
		void calculate_mep_from_charges();
		Icube& get_cube();
		Icube& calc_HOMO();
		Icube& calc_LUMO();
		Icube& calc_HOMO_density();
		Icube& calc_LUMO_density();
		Icube calc_HOMO_band(int bandn);
		Icube calc_LUMO_band(int bandn);
		Icube calc_band_EAS(int bandn);
		Icube calc_band_NAS(int bandn);
		Icube calc_EBLC_EAS();
		Icube calc_EBLC_NAS();
		Icube& grid_from_atoms(std::vector<double> values);
		void write_grid(); 
		void set_lim(double* Min, double* gridSides, int* gridSize);
		void redefine_lim(int atom,double size);
		void redefine_lim(double xc,double yc, double zc,int size);

};

#endif

//================================================================================
//END OF FILE
//================================================================================
