//pos_traj.h

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
#include <vector>
#include <string> 

class protein_lrd;
class primordia; 
//==================================================================
class traj_rd{
	public:	
		std::string user_path;
		std::vector<protein_lrd> frames;
		std::vector<int> res_list;
		std::vector< std::vector<double> > res_avg;
		std::vector< std::vector<double> > res_avg_all;
		std::vector< std::vector<double> > res_sd;
		
		std::vector< std::string > rds_labels;
		std::vector< std::string > atoms_labels;
		std::vector< std::vector< std::vector<double> > > atoms_rd; 
		//---------------------------------------------------
		traj_rd();
		traj_rd( std::vector<int> res );
		traj_rd( std::vector<primordia>& rds, std::vector<int>& ats, std::vector<int> res);
		~traj_rd();
		traj_rd(const traj_rd& rhs_traj) = delete;
		traj_rd& operator=(const traj_rd& rhs_traj) = delete;
		//-----------------------------------------------------
		void init_from_folder();
		void calculate_res_stats();
		void write_residues_reports();
		void residues_pattern_recognition();
		void gradient();
	 
};
