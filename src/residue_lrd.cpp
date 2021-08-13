//residue_lrd.cpp

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
// include statements from c++ library


#include <vector>
#include <string>
#include <fstream> 

#include "../include/Iprotein.h"
#include "../include/residue_lrd.h" 
#include "../include/Iline.h"

using std::move;
using std::string;
using std::vector;
using std::cout;
using std::endl;

/****************************************************/
residue_lrd::residue_lrd(){
	rd_sum.resize(19);
	rd_avg.resize(19);
}
/****************************************************/
residue_lrd::residue_lrd(const residue_lrd& rhs):
	rd_sum(rhs.rd_sum),
	rd_avg(rhs.rd_avg){
}
/****************************************************/
residue_lrd::residue_lrd(residue_lrd&& rhs) noexcept:
	rd_sum( move(rhs.rd_sum) ),
	rd_avg( move(rhs.rd_avg) ){
	
}
/****************************************************/
residue_lrd& residue_lrd::operator=(const residue_lrd& rhs){
	if (this!=&rhs){
		rd_sum = rhs.rd_sum;
		rd_avg = rhs.rd_avg;
	}
	return *this;
}
/****************************************************/
residue_lrd& residue_lrd::operator=(residue_lrd&& rhs) noexcept{
	if (this!=&rhs){
		rd_sum = move(rhs.rd_sum);
		rd_avg = move(rhs.rd_avg);
	}
	return *this;
}
/****************************************************/
residue_lrd::~residue_lrd(){}
//===================================================
/****************************************************/
protein_lrd::protein_lrd(){}
/****************************************************/
protein_lrd::protein_lrd(std::vector<residue_lrd> res_rd):
	residues_rd(res_rd)									{
	
	if ( residues_rd.size() > 0 ) {	
		protein_avg_avg.resize( residues_rd[0].rd_sum.size() );
		protein_sum_avg.resize( residues_rd[0].rd_sum.size());
		protein_min.resize( residues_rd[0].rd_sum.size() );
		protein_max.resize( residues_rd[0].rd_sum.size() );
	
		unsigned int i,j;
	
		for ( j=0; j<protein_min.size(); j++){
		protein_min[j] = residues_rd[0].rd_sum[j];
		protein_max[j] = residues_rd[0].rd_sum[j];
		}
	
		for( i=0; i<residues_rd.size(); i++){
			for ( j=0; j<residues_rd[i].rd_sum.size(); j++){
				protein_avg_avg[j] += residues_rd[i].rd_avg[j];
				protein_sum_avg[j] += residues_rd[i].rd_sum[j];
				if ( i > 0 ){
					if ( protein_min[j] > residues_rd[i-1].rd_sum[j] ){
						protein_min[j] = residues_rd[i].rd_sum[j];
					}
					if ( protein_max[j] < residues_rd[i-1].rd_sum[j] ){
						protein_max[j] = residues_rd[i].rd_sum[j];
					}
				}
			}
		}
	
		for ( j=0; j<residues_rd[0].rd_sum.size(); j++){
			protein_avg_avg[j] /= residues_rd.size();
			protein_sum_avg[j] /= residues_rd.size();
		}	
	}
	else{
		m_log->input_message("PDB information not well retrieved! \nVerify your input file!\n");
	}
}

/****************************************************/
protein_lrd::protein_lrd(const protein_lrd& rhs):
	residues_rd(rhs.residues_rd)				,
	labels(rhs.labels)							,
	protein_avg_avg(rhs.protein_avg_avg)		,
	protein_sum_avg(rhs.protein_sum_avg)		,
	protein_max(rhs.protein_max)				,
	protein_min(rhs.protein_min)				,
	hydrophobicity(rhs.hydrophobicity)			{
	
}
/****************************************************/
protein_lrd::protein_lrd(protein_lrd&& rhs) noexcept:
	residues_rd( move(rhs.residues_rd) )			,
	labels( move(rhs.labels) )						,
	protein_avg_avg( move(rhs.protein_avg_avg) )	,
	protein_sum_avg( move(rhs.protein_sum_avg) )	,
	protein_max( move(rhs.protein_max) )			,
	protein_min( move(rhs.protein_min) )			,
	hydrophobicity( move(rhs.hydrophobicity) )		{
}
/****************************************************/
protein_lrd& protein_lrd::operator=(const protein_lrd& rhs){
	if (this!=&rhs){
		residues_rd		= rhs.residues_rd;
		labels			= rhs.labels;
		protein_avg_avg	= rhs.protein_avg_avg;
		protein_sum_avg	= rhs.protein_sum_avg;
		protein_max		= rhs.protein_max;
		protein_min		= rhs.protein_min;
		hydrophobicity	= rhs.hydrophobicity;
	}
	return *this;
}
/****************************************************/
protein_lrd& protein_lrd::operator=(protein_lrd&& rhs) noexcept{
	if (this!=&rhs){
		residues_rd		= move(rhs.residues_rd);
		labels			= move(rhs.labels);
		protein_avg_avg	= move(rhs.protein_avg_avg);
		protein_sum_avg	= move(rhs.protein_sum_avg);
		protein_max		= move(rhs.protein_max);
		protein_min		= move(rhs.protein_min);
		hydrophobicity	= move(rhs.hydrophobicity);
	}
	return *this;
}
/****************************************************/
protein_lrd::~protein_lrd(){}
/****************************************************/
void protein_lrd::write_protein_lrd(const Iprotein& prot){
	
	
	unsigned int i,j;
	string name_f_file = get_file_name( prot.name.c_str() );
	name_f_file += "_residues.rslrd";
	std::ofstream lrd_file( name_f_file.c_str() );
	
	lrd_file.precision(6);
	lrd_file  << std::fixed;
	
	lrd_file << "#res Nucleophilicity Electrophilicity RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic Electrophilic Fukushima Electron_Density Softness_dual MEP\n";
	for( i=0; i<prot.residues.size(); i++ ){
		lrd_file << (i+1)						
				 << prot.residues[i].type	<< " ";
		for( j=0; j<residues_rd[i].rd_sum.size(); j++){
			lrd_file << residues_rd[i].rd_sum[j]<< " ";
		}
		lrd_file << "\n";
	}
	
	lrd_file	<<  "Residues_Average ";
	if ( residues_rd.size() > 0 ){
		for( j=0; j<residues_rd[0].rd_sum.size(); j++){
			lrd_file <<  protein_sum_avg[j]	<< " ";
		}	
	}
	lrd_file << endl;
	lrd_file.close();
	
}
/*******************************************************/
void protein_lrd::recalculate_stats(){
	
	unsigned int i,j;
	
	if ( residues_rd.size() > 0 ) {	
		protein_avg_avg.resize( residues_rd[0].rd_sum.size() );
		protein_sum_avg.resize( residues_rd[0].rd_sum.size() );
		protein_min.resize( residues_rd[0].rd_sum.size() );
		protein_max.resize( residues_rd[0].rd_sum.size() );
		
		for ( j=0; j<residues_rd[0].rd_sum.size(); j++){
			protein_min[j] = residues_rd[0].rd_sum[j];
			protein_max[j] = residues_rd[0].rd_sum[j];
		}
		
		for( i=0;i<residues_rd.size();i++){
			for ( j=0; j<residues_rd[i].rd_sum.size(); j++){
				protein_avg_avg[j] += residues_rd[i].rd_avg[j];
				protein_sum_avg[j] += residues_rd[i].rd_sum[j];
				if ( i > 0 ){
					if ( protein_min[j] > residues_rd[i-1].rd_sum[j] ){
						protein_min[j] = residues_rd[i].rd_sum[j];
					}
					if ( protein_max[j] < residues_rd[i-1].rd_sum[j] ){
						protein_max[j] = residues_rd[i].rd_sum[j];
					}
				}
			}
		}
		for ( j=0; j<residues_rd[i].rd_sum.size(); j++ ){
			protein_avg_avg[j] /= residues_rd.size();
			protein_sum_avg[j] /= residues_rd.size();
		}
	}else{
		m_log->input_message("PDB information not well retrieved! \nVerify your input file!\n");
	}
	
}
//===============================================================================
void protein_lrd::determine_hydrophilicity(){
	
	unsigned int i,j;
	string tmp = "";
	vector<string> aliphatic={"ALA","GLY","ILE","LEU","PRO","MET","PHE","VAL"}; 
	hydrophobicity.resize( labels.size() );
	for( i=0; i<labels.size()-1; i++){
		tmp = labels[i];
		tmp = tmp.substr(tmp.size()-3,tmp.size());
		for ( j=0;j<aliphatic.size(); j++){
			if ( tmp == aliphatic[j] )
				hydrophobicity[i] = 1;
		}
	}
}
//================================================================================
//END OF FILE
//================================================================================