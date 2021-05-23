//scripts.h

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

#ifndef SCRIPTS
#define SCRIPTS
//--------------------------------------------------------------------------
#include <string>
//--------------------------------------------------------------------------
/**
 * @class scripts
 * @author Igor Barden Grillo
 * @date 15/03/20
 * @file scripts.h
 * @brief Class to handle and write scripts for pos-processment by other programs
 */
 
class local_rd;
class traj_rd;
class ReactionAnalysis;

class scripts{
	public:
		//member variables
		std::string file_name;
		std::ofstream script_file;
		std::string s_type;
		//constructors/destructor
		scripts();
		scripts( std::string Nm, std::string _type );
		scripts(const scripts& rhs);
		scripts& operator=(const scripts& rhs);
		~scripts();
		//member functions
		void write_r_dos( std::vector<double>& energies );		
		void write_pymol_cube(local_rd& lrdVol, bool fixed);
		void write_pymol_pdb();
		void write_r_heatmap(std::vector< std::vector<double> > rd_numerical,std::vector<std::string> rds,std::vector<std::string> residues);
		void write_r_residuos_barplot();
		void write_r_reaction_analysis(traj_rd& path_rd, std::vector<std::string>& pair_labels,ReactionAnalysis& r_info,std::string& nameb );
};
#endif
//================================================================================
//END OF FILE
//================================================================================