// Imolecule.h
/* header file with the class declaration for chemical and geometrical information representation */

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

#ifndef IMOLECULE
#define IMOLECULE

//including c++ headers
#include <string>
#include <vector>

//including PRIMoRDiA headers
#include "../include/common.h"

class Iaorbital;
class Iatom; // foward declaration
//========================================================================
/**
 * @class Imolecule
 * @author Igor Barden Grillo
 * @date 20/03/18
 * @file Imolecule.h 
 * @brief Imolecule class to hold abstract representation for molecules and its quantum chemical informations
 * get from computational methods.
 * 
 * Main class of the program to hold quantum chemical molecular information, store atomic representations
 * objects.
 */
class Imolecule {
	public:
		//member variables
		std::string name;
		unsigned int num_of_atoms;
		unsigned int num_of_electrons;
		unsigned int num_of_ao; // number of atomic orbitals.
		int f_chg; // formal charge.
		float molar_mass;
		float mol_charge;
		double energy_tot; // total electronic energy. 
		double homo_energy;
		double lumo_energy;
		double total_dipmoment;
		double heat_of_formation;
		unsigned int MOnmb; //  number of alpha set molecular orbitals.
		unsigned int MOnmb_beta; // number of bet set molecular orbitals
		int homoN; // homo number
		int lumoN; // lumo number
		bool normalized; // if the atomic orbitals have the normalization factor calculated.
		bool bohr; // if the coordinates are in bohr.
		bool betad; // if the molecule has beta set of molecular orbitals.
		double ver_inf[3]; // coordinates of the inferior vertice.
		double ver_sup[3]; // coordinates of the superior vertice.
		double dipole_moment[3];
		std::vector<double> orb_energies;
		std::vector<double> orb_energies_beta;
		std::vector<double> coeff_MO;
		std::vector<double> coeff_MO_beta;
		std::vector<double> m_dens;
		std::vector<double> beta_dens;
		std::vector<double> m_overlap;
		std::vector <int> occupied;
		std::vector <int> occupied_beta;
		std::vector<Iatom> atoms;
		// constructors/destructor
		Imolecule(); 
		Imolecule(const Imolecule& rhs_molecule);
		Imolecule& operator=(const Imolecule& rhs_molecule);
		Imolecule(Imolecule&& rhs_molecule) noexcept;
		Imolecule& operator=(Imolecule&& rhs_molecule) noexcept;
		~Imolecule();
		//member funtions
		void add_atom(double x,double y,double z,std::string typ);
		void add_atom(Iatom atom);
		void print_coordinates();
		void write_xyz();
		void mol_vert_up();
		std::vector<double> extract_MO(int MO,bool beta);
		void print();
		void ang_to_bohr();
		void bohr_to_ang();
		void norm_orbs();
		int get_ao_number();
		double get_homo();
		double get_lumo();
		void print_basis();
		void update();
		void clear();
		bool check();
		double check_ed();
		void center_coord();
};

#endif
//================================================================================
//END OF FILE
//================================================================================