// source file for the trajectory pos-treatment class 
// pos_traj.cpp

// include statements from c++ library

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

#include <iostream>
#include <string> 
#include <cmath>
#include <cstring>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <fstream>
#include <experimental/filesystem>

#include "../include/common.h"
#include "../include/Iline.h"
#include "../include/pos_traj.h"
#include "../include/residue_lrd.h"
#include "../include/primordia.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/local_rd_cnd.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::to_string;
namespace fs = std::experimental::filesystem;
/****************************************************/
traj_rd::traj_rd()		:
	user_path("no_path"){	
}
/****************************************************/
traj_rd::traj_rd( vector<int> res ) :
	user_path("no_path")			,
	res_list(res)					{
	
	res_avg.resize(19);
	res_avg_all.resize(19);
	res_sd.resize(19);
	
	for( int i=0; i<res_list.size(); i++ ) res_list[i] -= 1;
}
/****************************************************/
traj_rd::traj_rd(	vector<primordia>& rds	,
					std::vector<int>& ats	, 
					vector<int> res )		:
					user_path("no_path")	,
					res_list(res)			{
	
	unsigned int i,j;
	
	string tmp_name  = "_";
	string tmp_name2 = "_";
	
	res_avg.resize(19);
	res_avg_all.resize(19);
	res_sd.resize(19);
	
	for( int i=0; i<res_list.size(); i++ ) {
		res_list[i] -= 1;
	}
		
	for( i=0; i<ats.size() ; i++){
		
		tmp_name  = rds[0].mol_info.atoms[ ats[i]-1 ].element;
		tmp_name += to_string( ats[i] );
		for ( int j=0; j<9; j++ ){
			atoms_labels.push_back(tmp_name);
		}		
		
		atoms_labels[0+i*9] +="\\nNucleophilicity";
		atoms_labels[1+i*9] +="\\nElectrophilicity";
		atoms_labels[2+i*9] +="\\nNetphilicity";
		atoms_labels[3+i*9] +="\\nSoftness";
		atoms_labels[4+i*9] +="\\nHardness\\n(LCP)";
		atoms_labels[5+i*9] +="\\nHardness\\n(Vee)";
		atoms_labels[6+i*9] +="\\nFukui\\nPotential";
		atoms_labels[7+i*9] +="\\nElectron\\nDensity";
		atoms_labels[8+i*9] +="\\nPartial\\nCharge";
	
		tmp_name2 = tmp_name;		
		tmp_name2+= "_Nphilicity";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_Ephilicity";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_NET";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_Soft";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_hard_lcp";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_hard_Vee";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_Fukui_pot";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_ED";
		rds_labels.push_back(tmp_name2);
		tmp_name2 = tmp_name;
		tmp_name2+="_chg";
		rds_labels.push_back(tmp_name2);
	}
	
	atoms_rd.resize( ats.size() );
	for( i=0; i<ats.size(); i++ ){
		atoms_rd[i].resize(9);
		for( j=0; j<9; j++ ){
			atoms_rd[i][j].resize( rds.size() );
		}
	}
	
	for ( i=0; i<atoms_rd.size(); i++ ){
		for( j=0; j<rds.size(); j++ ){
			atoms_rd[i][0][j] = rds[j].lrdCnd.lrds[0][ ats[i]-1 ];
			atoms_rd[i][1][j] = rds[j].lrdCnd.lrds[1][ ats[i]-1 ];
			atoms_rd[i][2][j] = rds[j].lrdCnd.lrds[3][ ats[i]-1 ];
			atoms_rd[i][3][j] = rds[j].lrdCnd.lrds[9][ ats[i]-1];
			atoms_rd[i][4][j] = rds[j].lrdCnd.lrds[5][ ats[i]-1 ];
			atoms_rd[i][5][j] = rds[j].lrdCnd.lrds[4][ ats[i]-1 ];
			atoms_rd[i][6][j] = rds[j].lrdCnd.lrds[6][ ats[i]-1 ];
			atoms_rd[i][7][j] = rds[j].lrdCnd.lrds[14][ ats[i]-1 ];
			atoms_rd[i][8][j] = rds[j].lrdCnd.lrds[13][ ats[i]-1 ];
		}
	}

	for ( i=0; i<rds.size(); i++){
		frames.push_back(rds[i].bio_rd);
	}
}
/****************************************************/
traj_rd::~traj_rd(){};
/****************************************************/
void traj_rd::init_from_folder(){
	
	fs::path c_path = fs::current_path();
	std::vector<string> fnames; 
	if ( user_path !="no_path" ){
		c_path =  user_path;
	}
	
	for ( const auto & entry : fs::directory_iterator(c_path) ){
		string tmp_name = entry.path();
		if ( check_file_ext( ".rslrd",tmp_name.c_str() ) ){
			fnames.push_back( tmp_name );
		}
	}
	//cout << c_path << endl;
	
	std::sort( fnames.begin(), fnames.end() );
	string Line;
	
	for( int i=0; i<fnames.size(); i++ ){
		protein_lrd Frame;
		if ( IF_file( fnames[i].c_str() ) ){
			std::ifstream fr_fle( fnames[i].c_str() );
			int line = 0;
			while( !fr_fle.eof() ){
				getline(fr_fle,Line);
				if ( line > 0 ){
					residue_lrd residue;
					Iline line_obj(Line);
					if ( line_obj.words[0].size() == 1 ){
						string lbl = line_obj.words[0];
						lbl +=  line_obj.words[1];
						Frame.labels.push_back(lbl);
						for( int j=2; j<line_obj.line_len; j++ ){
							residue.rd_sum[j-2] = line_obj.get_double(j);
						}
					}else{
						Frame.labels.push_back(line_obj.words[0]);
						for( int j=1; j<line_obj.line_len; j++ ){
							residue.rd_sum[j-1] = line_obj.get_double(j);
						}
					}
					Frame.residues_rd.push_back(residue);
				}
				line++;
			}
		if ( Frame.residues_rd.size() > 2 )
			frames.push_back(Frame);
		}
	}
}
/***************************************************************************/
void traj_rd::calculate_res_stats(){
		
	unsigned int i,j,k;
	
	for( i=0; i<15; i++ ){
		res_avg[i].resize( res_list.size() );
		res_sd[i].resize( res_list.size() );
		res_avg_all[i].resize( frames[0].residues_rd.size() );
	}
	
	int fsize = frames.size();
	
	for( k=0; k<res_list.size(); k++ ){
		for( j=0; j<15; j++ ){
			for( i=0; i<frames.size(); i++ ){
				res_avg[j][k] += frames[i].residues_rd[ res_list[k] ].rd_sum[j];
			}
			res_avg[j][k] /= fsize;
		}
	}
	
	for( i=0; i<fsize; i++ ){
		for( j=0; j<res_avg_all.size(); j++ ){
			for( k=0; k<frames[i].residues_rd.size(); k++ ){
				res_avg_all[j][k] += frames[i].residues_rd[k].rd_sum[j];
			}			
		}
	}
	
	for( i=0; i<res_avg_all.size(); i++){
		for( j=0; j<res_avg_all[i].size(); j++){
			res_avg_all[i][j] /= fsize;
		}
	}
	
	
	for( k=0; k<res_list.size(); k++ ){
		for( j=0; j<15; j++ ){
			for( i=0; i<frames.size(); i++ ){
				res_sd[j][k] +=((frames[i].residues_rd[res_list[k]].rd_sum[j] - res_avg[j][k])*
								(frames[i].residues_rd[res_list[k]].rd_sum[j] - res_avg[j][k]) );
			}
			res_sd[j][k] /= frames.size();
			res_sd[j][k] = sqrt(res_sd[j][k]);
		}
	}	 
}
/***************************************************************************/
void traj_rd::write_residues_reports(){
	string fname_res = "residues_data_frames";
	std::ofstream res_file_f( fname_res.c_str() );
	res_file_f	<< "frame Nucleophilicity Electrophilicity Radicality "
				<< "Netphilicity Hardness_Vee Hardness_LCP Fukui_pot_left Fukui_pot_right Fukui_pot_zero "
				<< "softness_dual hyper_softness Multiphilic Fukushima charge Electron_Density MEP "
				<< "hardness_TFD softness_avg hardness_int res\n";

	for( unsigned i=0;i<res_list.size();i++ ){
		for( unsigned j=0; j<frames.size(); j++ ){
			for( int k=0; k<frames[j].residues_rd[ res_list[i] ].rd_sum.size()+2; k++){
				if ( k == 0 ) { 
					res_file_f << j << " ";
				}
				res_file_f << frames[j].residues_rd[ res_list[i] ].rd_sum[k] << " ";
			}
			res_file_f << frames[j].labels[res_list[i]] << "\n";
		}
	}
	res_file_f.close();	
		
	string fname_res_avg = "residues_data_stat";
	std::ofstream res_avg_f( fname_res_avg.c_str() );	
	res_avg_f	<< "res Nucleophilicity Electrophilicity Radicality "
				<< "Netphilicity Hardness_Vee Hardness_LCP Fukui_pot_left Fukui_pot_right Fukui_pot_zero "
				<< "softness_dual hyper_softness Multiphilic Fukushima charge Electron_Density MEP "
				<< "hardness_TFD softness_avg hardness_int type\n";
		
	for( unsigned i=0; i<res_list.size(); i++){
		res_avg_f << frames[0].labels[ res_list[i] ] << " ";
		for( unsigned j=0; j<19; j++ ){
			res_avg_f << res_avg[j][i] << " ";
		}
		res_avg_f << "AVG \n";
	}
	for( unsigned i=0; i<res_list.size(); i++ ){
		res_avg_f << frames[0].labels[ res_list[i] ] << " ";
		for( unsigned j=0; j<19; j++ ){
			res_avg_f << res_sd[j][i] << " ";
		}
		res_avg_f << "SD \n";
	}
	
	res_avg_f.close();
	
	string fname_pro_avg = "protein_data_stat";
	std::ofstream pro_avg( fname_pro_avg.c_str() );
	pro_avg		<< "res Nucleophilicity Electrophilicity Radicality "
				<< "Netphilicity Hardness_Vee Hardness_LCP Fukui_pot_left Fukui_pot_right Fukui_pot_zero "
				<< "softness_dual hyper_softness Multiphilic Fukushima charge Electron_Density MEP "
				<< "hardness_TFD softness_avg hardness_int type\n";
	
	for( unsigned i=0; i<frames[0].residues_rd.size(); i++ ){
		pro_avg << frames[0].labels[i] << " ";
		for(int j=0; j<19; j++ ){
			pro_avg << res_avg_all[j][i] << " ";
		}
		pro_avg << "\n";
	}
	pro_avg.close();
	
}
/******************************************************************/
void traj_rd::gradient(){
	for (unsigned int i=0;i<frames.size();i++){
		frames[i].determine_hydrophilicity();
	}
	
	vector< vector<double> > grad_hidrophilic;
	vector< vector<double> > grad_hidrophobic;
	vector< vector<double> > grad_total;
	
	grad_hidrophilic.resize(15);
	grad_hidrophobic.resize(15);
	grad_total.resize(15);
	for( int i=0; i<15; i++ ){
		grad_hidrophilic[i].resize( frames.size() );
		grad_hidrophobic[i].resize( frames.size() );
		grad_total[i].resize( frames.size() );
	}
	
	for ( unsigned f=0;f<frames.size();f++){
		for ( int rd=0;rd<16;rd++ ){
			for ( unsigned rs=0;rs<frames[0].residues_rd.size()-1;rs++ ){
				if ( f != frames.size()-1){				
					grad_total[rd][f] += frames[f+1].residues_rd[rs].rd_sum[rd]
										-frames[f].residues_rd[rs].rd_sum[rd];
					if ( frames[f].hydrophobicity[rs] > 0 ){
						grad_hidrophobic[rd][f] += frames[f+1].residues_rd[rs].rd_sum[rd]
												-frames[f].residues_rd[rs].rd_sum[rd];
					}else{
						grad_hidrophilic[rd][f] += frames[f+1].residues_rd[rs].rd_sum[rd]
												-frames[f].residues_rd[rs].rd_sum[rd];
					}
				}
			}
		}
	}
	
	string fname_res = "gradient_total";
	std::ofstream res_file_tot( fname_res.c_str() );
	res_file_tot << "frame EAS NAS RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrophilic Fukushima Electron_Density Softness_dual mep comp_hardness\n";

	fname_res = "gradient_hidrophobic";
	std::ofstream res_file_phobic( fname_res.c_str() );
	res_file_phobic << "frame EAS NAS RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrophilic Fukushima Electron_Density Softness_dualmep comp_hardness\n";

	fname_res = "gradient_hidrophilic";
	std::ofstream res_file_philic( fname_res.c_str() );
	res_file_philic << "frame EAS NAS RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrophilic Fukushima Electron_Density Softness_dual mep comp_hardness\n";

	for( unsigned i=0; i<frames.size()-1; i++ ){
		for( int j=0; j<16; j++){
			if ( j == 0 ){
				int ii = i+1;
				res_file_tot 	<< ii <<" ";
				res_file_phobic << ii <<" ";
				res_file_philic << ii <<" ";
			}
			res_file_tot << grad_total[j][i] << " ";
			res_file_phobic << grad_hidrophobic[j][i] << " ";
			res_file_philic << grad_hidrophilic[j][i] << " ";
		}
		res_file_tot << "\n";
		res_file_phobic << "\n";
		res_file_philic << "\n";
	}
	
	res_file_tot.close();
	res_file_phobic.close();
	res_file_philic.close();

}
/******************************************************/

//////////////////////////////////////////////////////