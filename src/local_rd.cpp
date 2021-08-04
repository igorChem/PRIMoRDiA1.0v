// source file for the local_rd class 
// local_rd.cpp

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
#include <memory>
#include <omp.h>
#include <algorithm>
#include <vector>
// include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/Icube.h"
#include "../include/global_rd.h"
#include "../include/local_rd.h"

using std::cout;
using std::endl;
using std::move;
using std::unique_ptr;
using std::string;
using std::vector;

//==========================================================================
//Class member functions definitions
/****************************************************************************/
local_rd::local_rd()	:
	name("nonamed")		,
	finite_diff(false)	,
	locHardness(false)	,
	charge(1)			{
		
	rd_names = {
		"HOMO",						// 0  bs
		"LUMO",						// 1  bs
		"Elec_Dens",				// 2  bs
		"ED_cation",				// 3  bs
		"ED_anion", 				// 4  bs
		"Nucleophilicity",			// 5  df
		"Electrophilicity",			// 6  df
		"Radical_sucseptibility",	// 7  df
		"Netphilicity",				// 8  df
		"L_Hardness_LCP",			// 9  if density
		"L_Hardness_Vee",			// 10 if density
		"Fukui_Potential_left",		// 11 df 
		"Fukui_Potential_right",	// 12 ex
		"Fukui_Potential_avg",		// 13 ex
		"L_Hardness_INT",			// 14 ex
		"Softness_Dual",			// 15 df
		"Softness_AVG",				// 16 ex
		"Hyper_Softess",			// 17 ex
		"Multiphilicity"			// 18 ex
		"MEP",						// 19 if density
		"MO_BAND"					// 20
	};
}
/***********************************************************************************/
local_rd::local_rd(const Icube& HOmo		,
				   const Icube& LUmo)		:
	name(HOmo.name)							,
	finite_diff(false)						,
	locHardness(false)						{
		
	
	lrds[0] = HOmo;
	lrds[1] = LUmo;
	lrds[5] = lrds[0]*lrds[0];
	lrds[6] = lrds[1]*lrds[1];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];
	
	
	rd_names = {
		"HOMO",						// 0  bs
		"LUMO",						// 1  bs
		"Elec_Dens",				// 2  bs
		"ED_cation",				// 3  bs
		"ED_anion", 				// 4  bs
		"Nucleophilicity",			// 5  df
		"Electrophilicity",			// 6  df
		"Radical_sucseptibility",	// 7  df
		"Netphilicity",				// 8  df
		"L_Hardness_LCP",			// 9  if density
		"L_Hardness_Vee",			// 10 if density
		"Fukui_Potential_left",		// 11 df 
		"Fukui_Potential_right",	// 12 ex
		"Fukui_Potential_avg",		// 13 ex
		"L_Hardness_INT",			// 14 ex
		"Softness_Dual",			// 15 df
		"Softness_AVG",				// 16 ex
		"Hyper_Softess",			// 17 ex
		"Multiphilicity"			// 18 ex
		"MEP",						// 19 if density
		"MO_BAND"					// 20 if band
	};

}
/***********************************************************************************/
local_rd::local_rd(const Icube& elec_dens	,
				   const Icube& HOmo		,
				   const Icube& LUmo		):
	name(HOmo.name)							,
	finite_diff(false)						,
	locHardness(true)						,
	charge(1)								{
		
	lrds[0] = HOmo;
	lrds[1] = LUmo;
	lrds[2] = elec_dens;
	lrds[5] = lrds[0];
	lrds[6] = lrds[1];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];
	
	
	rd_names = {
		"HOMO",						// 0  bs
		"LUMO",						// 1  bs
		"Elec_Dens",				// 2  bs
		"ED_cation",				// 3  bs
		"ED_anion", 				// 4  bs
		"Nucleophilicity",			// 5  df
		"Electrophilicity",			// 6  df
		"Radical_sucseptibility",	// 7  df
		"Netphilicity",				// 8  df
		"L_Hardness_LCP",			// 9  if density
		"L_Hardness_Vee",			// 10 if density
		"Fukui_Potential_left",		// 11 df 
		"Fukui_Potential_right",	// 12 ex
		"Fukui_Potential_avg",		// 13 ex
		"L_Hardness_INT",			// 14 ex
		"Softness_Dual",			// 15 df
		"Softness_AVG",				// 16 ex
		"Hyper_Softess",			// 17 ex
		"Multiphilicity"			// 18 ex
		"MEP",						// 19 if density
		"MO_BAND"					// 20 if band
	};
}
/***********************************************************************************/
local_rd::local_rd(const Icube& elecDens	,
					const Icube& cationDens	, 
					const Icube& anionDens	,
					int chg					):
	name(elecDens.name)						,
	finite_diff(true)						,
	locHardness(true)						,
	charge(chg)								{

	
	lrds[2] = elec_dens;
	lrds[3] = cationDens;
	lrds[4] = anionDens;
	lrds[5] = lrds[2] - lrds[3];
	lrds[6] = lrds[4] - lrds[2];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];	
	
	rd_names = {
		"HOMO",						// 0  bs
		"LUMO",						// 1  bs
		"Elec_Dens",				// 2  bs
		"ED_cation",				// 3  bs
		"ED_anion", 				// 4  bs
		"Nucleophilicity",			// 5  df
		"Electrophilicity",			// 6  df
		"Radical_sucseptibility",	// 7  df
		"Netphilicity",				// 8  df
		"L_Hardness_LCP",			// 9  if density
		"L_Hardness_Vee",			// 10 if density
		"Fukui_Potential_left",		// 11 df 
		"Fukui_Potential_right",	// 12 ex
		"Fukui_Potential_avg",		// 13 ex
		"L_Hardness_INT",			// 14 ex
		"Softness_Dual",			// 15 df
		"Softness_AVG",				// 16 ex
		"Hyper_Softess",			// 17 ex
		"Multiphilicity"			// 18 ex
		"MEP",						// 19 if density
		"MO_BAND",					// 20 if band
		"L_Hardness_TFD"			// 21 if density
	};
}
/***********************************************************************************/
local_rd::local_rd(const local_rd& lrd_rhs)	:
	name(lrd_rhs.name)						,
	finite_diff(lrd_rhs.finite_diff)		,
	locHardness(lrd_rhs.locHardness)		,
	charge(lrd_rhs.charge)					,
	rd_names(lrd_rhs.rd_names)				,
	lrds(lrd_rhs.lrds)						{
}
/***********************************************************************************/
local_rd& local_rd::operator=(const local_rd& lrd_rhs){
	if( this!=&lrd_rhs ){
		name			= lrd_rhs.name;
		finite_diff		= lrd_rhs.finite_diff;
		charge			= lrd_rhs.charge;
		locHardness		= lrd_rhs.locHardness;
		rd_names		= lrd_rhs.rd_names;
		lrds			= lrd_rhs.lrds;
	}
	return *this;
}
/***********************************************************************************/
local_rd::local_rd(local_rd&& lrd_rhs) noexcept	:
	name(lrd_rhs.name)							,
	finite_diff(lrd_rhs.finite_diff)			,
	locHardness(lrd_rhs.locHardness)			,
	charge(lrd_rhs.charge)						,
	rd_names( move(lrd_rhs.rd_names) )			,
	lrds( move(lrd_rhs.lrds) )					{
}
/***********************************************************************************/
local_rd& local_rd::operator=(local_rd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs ){
		name			= move(lrd_rhs.name);
		finite_diff		= lrd_rhs.finite_diff;
		charge			= lrd_rhs.charge;
		locHardness		= lrd_rhs.locHardness;
		rd_names		= move(rd_names);
		lrds			= move(lrds);
	}
	return *this;
}
/***********************************************************************************/
void local_rd::calculate_fukui_Band(const Icube& homo_b, const Icube& lumo_b){
	lrds[5] = homo_b;
	lrds[6] = lumo_b;
	lrds[20] = (homo_b+lumo_b)/2.0;
	lrds[5].normalize(5);
	lrds[6].normalize(5);
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];
}
/***********************************************************************************/
void local_rd::calculate_RD(const global_rd& grd){
	lrds[14] = lrds[5]*grd.grds[5] - lrds[5]*grd.grds[6];
	lrds[15] = lrds[8]*grd.grds[9];
	lrds[16] = lrds[7]*grd.grds[9];
	lrds[17] = lrds[8]*(grd.grds[9]*grd.grds[9]);
	lrds[18] = lrds[8]*grd.grds[10];
}
/***********************************************************************************/
void local_rd::calculate_Fukui_potential(const Imolecule& mol){
	
	vector<double> elec_H(lrds[5].voxelN);
	vector<double> nuc_H(lrds[5].voxelN);
	vector<double> rad_H(lrds[5].voxelN);
	unsigned int i,j,k,x,y,z;
	double ii	= 0;
	double jj	= 0;
	double kk	= 0;
	double xx	= 0;
	double yy	= 0;
	double zz	= 0;
	double r	= 0;
	double s1	= lrds[5].gridsides[0];
	double s2	= lrds[5].gridsides[1];
	double s3	= lrds[5].gridsides[2];
	double o1	= lrds[5].origin[0];
	double o2	= lrds[5].origin[1];
	double o3	= lrds[5].origin[2];
	unsigned int g1 = lrds[5].grid[0];
	unsigned int g2 = lrds[5].grid[1];
	unsigned int g3 = lrds[5].grid[2];
		
	#pragma omp declare reduction(vec_d_plus : std::vector<double> : \
				std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
				initializer(omp_priv = omp_orig)

	#pragma omp parallel default(shared) private(i,j,k,x,y,z,r,xx,yy,zz,ii,jj,kk)
	{
	#pragma omp parallel for reduction(vec_d_plus:elec_H,nuc_H,rad_H)
	for(i=0;i<g1;i++){
		for(j=0;j<g2;j++){
			for(k=0;k<g3;k++){
				for (x=0;x<g1;x++){
					for (y=0;y<g2;y++){
						for (z=0;z<g3;z++){
							xx	= x*s1 + o1;
							yy	= y*s2 + o2;
							zz	= z*s3 + o3;
							ii	= i*s1 + o1;
							jj	= j*s2 + o2;
							kk	= k*s3 + o3;
							xx 	-= ii;
							xx 	*= xx;
							yy 	-= jj;
							yy 	*= yy;
							zz 	-= kk;
							zz 	*= zz;
							r  	= sqrt(xx + yy + zz);
							if ( r == 0.000 ){
								elec_H[i*g1*g1+j*g2+k]	+= 0;
								nuc_H[i*g1*g1+j*g2+k] 	+= 0;
								rad_H[i*g1*g1+j*g2+k] 	+= 0;
							}else{
								elec_H[i*g1*g1+j*g2+k] += lrds[5].scalar[x*g1*g1+y*g2+z]/r;
								nuc_H[i*g1*g1+j*g2+k] += lrds[6].scalar[x*g1*g1+y*g2+z]/r;
								rad_H[i*g1*g1+j*g2+k] += lrds[7].scalar[x*g1*g1+y*g2+z]/r;
							}
						}
					}
				}
			}
		}
	}
	}
	
	double volume = std::abs(s1*s2*s3);
	lrds[11] = lrds[5];
	lrds[12] = lrds[5];
	lrds[13] = lrds[5];
	for(i=0;i<elec_H.size();i++) { lrds[11].scalar[i] = elec_H[i]; }
	for(i=0;i<nuc_H.size();i++) { lrds[12].scalar[i] = nuc_H[i]; }
	for(i=0;i<rad_H.size();i++) { lrds[13].scalar[i] = rad_H[i]; }
	
	lrds[11] = lrds[11]*volume;
	lrds[12] = lrds[12]*volume;
	lrds[13] = lrds[13]*volume;
	
}
/***********************************************************************************/
void local_rd::calculate_hardness(const global_rd& grd){
	locHardness = true;
	
	//Local Chemical Potential method
	double numofelec	= lrds[2].calc_cube_integral();
	double val1			= grd.grds[7]/numofelec;
	double val2			= grd.grds[8]*2;
		
	Icube elec_dens_norm= lrds[2]/ numofelec ;
	Icube elec_dens_hard= elec_dens_norm*val2;
	Icube homo_elec_dens= lrds[0].normalize() - elec_dens_norm;
	homo_elec_dens		= homo_elec_dens*val1;
	lrds[9]				= homo_elec_dens + elec_dens_hard;
	

	//local hardness com aproximação de potencial elétron-elétron
	vector<double> elec_H(lrds[5].voxelN);
	unsigned int i,j,k,x,y,z;
	double ii	= 0;
	double jj	= 0;
	double kk	= 0;
	double xx	= 0;
	double yy	= 0;
	double zz	= 0;
	double r	= 0;
	double s1	= lrds[5].gridsides[0];
	double s2	= lrds[5].gridsides[1];
	double s3	= lrds[5].gridsides[2];
	double o1	= lrds[5].origin[0];
	double o2	= lrds[5].origin[1];
	double o3	= lrds[5].origin[2];
	unsigned int g1 = lrds[5].grid[0];
	unsigned int g2 = lrds[5].grid[1];
	unsigned int g3 = lrds[5].grid[2];
	
	#pragma omp declare reduction(vec_d_plus : std::vector<double> : \
				std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                initializer(omp_priv = omp_orig)

	#pragma omp parallel default(shared) private(i,j,k,x,y,z,r,xx,yy,zz,ii,jj,kk)
	{
	#pragma omp parallel for reduction(vec_d_plus:elec_H)
	for(i=0;i<g1;i++){
		for(j=0;j<g2;j++){
			for(k=0;k<g3;k++){
				for (x=0;x<g1;x++){
					for (y=0;y<g2;y++){
						for (z=0;z<g3;z++){
							xx	= x*s1 + o1;
							yy	= y*s2 + o2;
							zz	= z*s3 + o3;
							ii	= i*s1 + o1;
							jj	= j*s2 + o2;
							kk	= k*s3 + o3;
							xx	-= ii;
							xx 	*= xx;
							yy 	-= jj;
							yy 	*= yy;
							zz 	-= kk;
							zz 	*= zz;
							r  	= sqrt(xx + yy + zz);
							if ( r == 0.000 ){
								elec_H[i*g1*g1+j*g2+k] += 0;
							}else{
								elec_H[i*g1*g1+j*g2+k] += lrds[2].scalar[x*g1*g1+y*g2+z]/r;
							}
						}
					}
				}
			}
		}
	}
	}
	
	double volume = std::abs(s1*s2*s3);
	lrds[10] = lrds[5];
	for(i=0;i<elec_H.size();i++) { lrds[10].scalar[i] = elec_H[i]; }
	lrds[10] = lrds[10]*volume;
	lrds[10] = lrds[10]*(1/numofelec);
	//----------------------------------------------------------------------------
	
	density_rc	= lrds[2].scale_cube(0.3333333);
	double Ck	= (3/10)*pow((3*M_PI*M_PI),2/3);
	double Cx	= (3/4*M_PI)*pow((3*M_PI*M_PI),1/3);
	lrds[21]	= (2/9*numofelec)*density_rc;
	Icube temp1	= density_rc;
	Icube temp2	= density_rc;
	Icube temp3	= density_rc;
	temp1 = temp1*5*Ck - 2*Cx;
	temp2 = 0.458*density_rc;
	temp3 = temp2 + 1; 
	temp3 = temp3.scale_cube(3.0);
	
	temp2 = temp2+2;
	temp2 = temp2/temp3;
	temp2 = -0.0466*temp2;
	lrds[21] = lrds[21]*(temp1-temp2)
	lrds[21] = lrds[21]+lrds[10];
}
/***********************************************************************************/
local_rd operator-(const local_rd& lrd_lhs,const local_rd& lrd_rhs){
	local_rd Result(lrd_lhs);
	for(unsigned int i=0; i<lrd_lhs.lrds.size(); i++){
		Result.lrds[i] = lrd_lhs.lrds[i] - lrd_rhs.lrds[i];
	}
	return Result;
}
/***********************************************************************************/
void local_rd::write_LRD(){
	
	std::string typestr;
	std::string typestr2;
	
	if ( finite_diff )  {
		typestr = "Descriptor calculated with Finite Differences approximation \n";
		typestr2 = "FD";
	}else{
		typestr = " Descriptor calculated with Frozen Orbital Approximation \n";
		typestr2 = "FOA";
	}
	
	string name_type = name + typestr2;
	std::vector<string>	cube_names;
	rd_names = {
		"HOMO",						// 0  bs
		"LUMO",						// 1  bs
		"Elec_Dens",				// 2  bs
		"ED_cation",				// 3  bs
		"ED_anion", 				// 4  bs
		"Nucleophilicity",			// 5  df
		"Electrophilicity",			// 6  df
		"Radical_sucseptibility",	// 7  df
		"Netphilicity",				// 8  df
		"L_Hardness_LCP",			// 9  if density
		"L_Hardness_Vee",			// 10 if density
		"Fukui_Potential_left",		// 11 df 
		"Fukui_Potential_right",	// 12 ex
		"Fukui_Potential_avg",		// 13 ex
		"L_Hardness_INT",			// 14 ex
		"Softness_Dual",			// 15 df
		"Softness_AVG",				// 16 ex
		"Hyper_Softess",			// 17 ex
		"Multiphilicity"			// 18 ex
		"MEP",						// 19 if density
		"MO_BAND",					// 20 if band
		"L_Hardness_TFD"			// 21 if density
	};
	
	cube_names.push_back(name_type+"_left_Fukui");
	cube_names.push_back(name_type+"_right_Fukui");
	cube_names.push_back(name_type+"_zero_Fukui");
	cube_names.push_back(name_type+"_net_Fukui_ph1");
	cube_names.push_back(name_type+"_net_Fukui_ph2");
	cube_names.push_back(name_type+"_left_Fukui");
	cube_names.push_back(name_type+"_left_Fukui");
	cube_names.push_back(name_type+"_left_Fukui");
	cube_names.push_back(name_type+"_left_Fukui");
	string fukui_suc_elec	= name + typestr2 + "_left_Fukui"; 
	string fukui_suc_nuc	= name + typestr2 + "_right_Fukui"; 
	string fukui_suc_rad	= name + typestr2 + "_zero_Fukui"; 
	string deltaFukui		= name + typestr2 + "_net_Fukui_ph1";
	string deltaFukui2		= name + typestr2 + "_net_Fukui_ph2";
	string local_hardnessA	= name + typestr2 + "_hardness_Vee";
	string local_hardnessB	= name + typestr2 + "_hardness_LCP";
	
	lrds[5].header	= "Left Fukui function, electrophilic attack succescitibily or local electofilicity\n"	+ typestr;
	lrds[6].header	= "Right Fukui Function, nucleophilic attack succescitibily or local nucleofilicity\n"	+ typestr;
	lrds[7].header	= "Average Fukui Function, radical attack succescitibily\n"								+ typestr;
	lrds[8].header	= "Net Fukui Function, dual descriptor of attack succescitibility\n"					+ typestr;
	
	string loc_softdual		= name + typestr2 + "_softdual"; 
	string local_mult		= name + typestr2 + "_multifilic";
	Softness_Dual.header	= "Local sofntess dual \n" + typestr;
	multifilic.header		= "Multifilic Descriptor \n" + typestr;
	
	
	EAS.write_cube(fukui_suc_elec + ".cube");
	EAS.name = fukui_suc_elec;
	NAS.write_cube(fukui_suc_nuc + ".cube");
	NAS.name = fukui_suc_nuc;
	RAS.write_cube(fukui_suc_rad + ".cube");
	RAS.name = fukui_suc_rad;
	Dual.write_cube(deltaFukui + ".cube");
	Dual.name = deltaFukui;
	Dual.write_cube(deltaFukui2 + ".cube");
	
	if ( extra_RD ){
		multifilic.write_cube(local_mult + ".cube");
		multifilic.write_cube(local_mult + "2_.cube");
		multifilic.name = local_mult;
		Softness_Dual.write_cube(loc_softdual + ".cube");
		Softness_Dual.write_cube(loc_softdual + "_2.cube");
		Softness_Dual.name = loc_softdual;
	}
	
	if( locHardness == true ) {
		Hardness.write_cube(local_hardness + ".cube");
		Hardness.name = local_hardness;
	}
	m_log->input_message("Finishing writting the local reactivity descriptos grids.\n");
}
/***********************************************************************************/
local_rd::~local_rd(){}
//================================================================================
//END OF FILE
//================================================================================

