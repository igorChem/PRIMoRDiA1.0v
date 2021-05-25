//comp_hardness.h

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
#include <fstream>
#include "../include/comp_hardness.h"
#include "../include/local_rd_cnd.h"
#include "../include/global_rd.h"
#include "../include/Imolecule.h"
#include "../include/Iatom.h"
#include "../include/Iprotein.h"
#include "../include/residue_lrd.h"

using std::move;
/********************************************************************/
comp_hard::comp_hard():
	g_comp_hard(0.00) {
}
/********************************************************************/
comp_hard::comp_hard(const global_rd& grd,
					 double vm			){
	g_comp_hard= grd.hardness/vm;
	
}
/********************************************************************/
comp_hard::comp_hard(const global_rd& grd	, 
					 const local_rd_cnd& lrd,
					 double vm)				{
						 
	g_comp_hard= grd.hardness/vm;
	l_comp_hard.resize(4);
	for ( int i=0;i<4;i++){
		l_comp_hard[i].resize( lrd.EAS.size() );
	}
	
	for ( int i=0; i<lrd.EAS.size(); i++ ){
		l_comp_hard[0][i] = lrd.hardness_A[i]/lrd.molecule.atoms[i].wdw_volume; 	
		l_comp_hard[1][i] = lrd.hardness_B[i]/lrd.molecule.atoms[i].wdw_volume;	
		l_comp_hard[2][i] = lrd.hardness_C[i]/lrd.molecule.atoms[i].wdw_volume;	
		l_comp_hard[3][i] = lrd.hardness_D[i]/lrd.molecule.atoms[i].wdw_volume;				
	}	

}
/********************************************************************/
comp_hard::~comp_hard(){}
/********************************************************************/
comp_hard::comp_hard(const comp_hard& rhs)	:
	g_comp_hard( rhs.g_comp_hard) 			,
	l_comp_hard( rhs.l_comp_hard) 			,
	l_comp_hard_bio( rhs.l_comp_hard_bio)	{	
}
/********************************************************************/
comp_hard& comp_hard::operator=(const comp_hard& rhs){
	if ( this != &rhs ){
		g_comp_hard 	=rhs.g_comp_hard;
		l_comp_hard 	=rhs.l_comp_hard;
		l_comp_hard_bio =rhs.l_comp_hard_bio;
	}
	return *this;
}
/********************************************************************/
comp_hard::comp_hard( comp_hard&& rhs ) noexcept:
	g_comp_hard( move(rhs.g_comp_hard) )		,
	l_comp_hard( move(rhs.l_comp_hard) )		,
	l_comp_hard_bio( move(rhs.l_comp_hard_bio) ){
	
}
/********************************************************************/
comp_hard& comp_hard::operator=( comp_hard&& rhs ) noexcept{
	if ( this != &rhs ){
		g_comp_hard 	= move(rhs.g_comp_hard);
		l_comp_hard 	= move(rhs.l_comp_hard);
		l_comp_hard_bio = move(rhs.l_comp_hard_bio);

	}
	return *this;
}
/********************************************************************/
void comp_hard::calculate_protein(	const protein_lrd& lrd,
						const Iprotein& prot){

	double mol_vol = 0.000;
	l_comp_hard_bio.resize(4);
	unsigned int res_s = prot.residues.size();
	for( int i=0; i<l_comp_hard_bio.size(); i++ ){
		l_comp_hard_bio[i].resize( res_s );
	}
	for ( int j=0; j<res_s; j++){
		mol_vol = prot.residues[j].molar_vol;
		l_comp_hard_bio[0][j] = lrd.residues_rd[j].rd_sum[5]/mol_vol;
		l_comp_hard_bio[1][j] = lrd.residues_rd[j].rd_sum[6]/mol_vol;
		l_comp_hard_bio[2][j] = lrd.residues_rd[j].rd_sum[7]/mol_vol;
		l_comp_hard_bio[3][j] = lrd.residues_rd[j].rd_sum[8]/mol_vol;
	}
	std::string local_name_bio	= prot.name;
	local_name_bio +=".local_bio_CH";
	std::ofstream loc_nm_file( local_name_bio.c_str() );

	loc_nm_file << "#res CH_Hardness_A CH_Hardness_B CH_Hardness_C CH_Hardness_D\n";
	for( int i=0; i<res_s ;i++ ){
		loc_nm_file << (i+1)
					<< prot.residues[i].type << " "
					<< l_comp_hard_bio[0][i] << " "
					<< l_comp_hard_bio[1][i] << " "
					<< l_comp_hard_bio[2][i] << " "
					<< l_comp_hard_bio[3][i] << "\n";
	}
	loc_nm_file.close();
}
/********************************************************************/
void comp_hard::write_comp_hardness(const char* name){
	std::string glob_name = name;
	std::string local_name= name;
	glob_name+=  ".glob_CH";
	local_name+= ".local_CH";

	std::ofstream glob_file( glob_name.c_str() );
	std::ofstream local_file( local_name.c_str() );

	glob_file 	<< "Composite hardness\n"
				<< g_comp_hard;

	glob_file.close();

	local_file << "#atom CH_Hardness_A CH_Hardness_B CH_Hardness_C CH_Hardness_D\n";
	for( int i=0; i<l_comp_hard[0].size() ;i++ ){
		local_file << (i+1)
				<< l_comp_hard[0][i] << " "
				<< l_comp_hard[1][i] << " "
				<< l_comp_hard[2][i] << " "
				<< l_comp_hard[3][i] << "\n";
	}
	local_file.close();

}

//===================================================================