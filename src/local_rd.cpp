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
#include "../include/Itimer.h"

using std::cout;
using std::endl;
using std::move;
using std::unique_ptr;
using std::string;
using std::vector;

std::vector<std::string> descriptor_names = {
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
	"Multiphilicity",			// 18 ex
	"MEP",						// 19 if density
	"MO_BAND",					// 20 if Band
	"L_Hardness_TFD"			// 21 if TFD
};


//==========================================================================
//Class member functions definitions
/****************************************************************************/
local_rd::local_rd()			:
	name("nonamed")				,
	FD(false)					,
	LH(false)					,
	band(false)					,
	charge(1)					,
	TFD(false)					,
	rd_names(descriptor_names)	{
		
	lrds.resize( rd_names.size() );
}
/***********************************************************************************/
local_rd::local_rd(Icube HOmo		,
				   Icube LUmo)		:
	name(HOmo.name)							,
	FD(false)								,
	LH(false)								,
	band(false)								,
	charge(1)								,
	TFD(false)								,
	rd_names(descriptor_names)				{
	
	lrds.resize( rd_names.size() );
	
	lrds[0] = HOmo;
	lrds[1] = LUmo;
	lrds[5] = lrds[0]*lrds[0];
	lrds[6] = lrds[1]*lrds[1];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];
}
/***********************************************************************************/
local_rd::local_rd(Icube elec_dens	,
					Icube HOmo		,
					Icube LUmo		):
	name(HOmo.name)							,
	FD(false)								,
	LH(false)								,
	band(false)								,
	charge(1)								,
	TFD(false)								,
	rd_names(descriptor_names)				{
	
	lrds.resize( rd_names.size() );

	lrds[0] = HOmo;
	lrds[1] = LUmo;
	lrds[2] = elec_dens;
	lrds[5] = lrds[0]*lrds[0];
	lrds[6] = lrds[1]*lrds[1];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];
	
}
/***********************************************************************************/
local_rd::local_rd(Icube elecDens			,
					Icube cationDens		, 
					Icube anionDens			,
					int chg					):
	name(elecDens.name)						,
	FD(true)								,
	LH(true)								,
	charge(chg)								,
	TFD(false)								,
	band(false)								,
	rd_names(descriptor_names)				{
	
	lrds.resize( rd_names.size() );

	lrds[2] = elecDens;
	lrds[3] = cationDens;
	lrds[4] = anionDens;
	lrds[5] = lrds[2] - lrds[3];
	lrds[6] = lrds[4] - lrds[2];
	lrds[7] = (lrds[5]+lrds[6])/2.0;
	lrds[8] = lrds[6]-lrds[5];	

}
/***********************************************************************************/
local_rd::local_rd(const local_rd& lrd_rhs)	:
	name(lrd_rhs.name)						,
	FD(lrd_rhs.FD)							,
	LH(lrd_rhs.LH)							,
	charge(lrd_rhs.charge)					,
	TFD(lrd_rhs.TFD)						,
	rd_names(lrd_rhs.rd_names)				,
	lrds(lrd_rhs.lrds)						{
}
/***********************************************************************************/
local_rd& local_rd::operator=(const local_rd& lrd_rhs){
	if( this!=&lrd_rhs ){
		name			= lrd_rhs.name;
		FD				= lrd_rhs.FD;
		charge			= lrd_rhs.charge;
		LH				= lrd_rhs.LH;
		TFD				= lrd_rhs.TFD;
		rd_names		= lrd_rhs.rd_names;
		lrds			= lrd_rhs.lrds;
	}
	return *this;
}
/***********************************************************************************/
local_rd::local_rd(local_rd&& lrd_rhs) noexcept	:
	name(lrd_rhs.name)							,
	FD(lrd_rhs.FD)								,
	LH(lrd_rhs.LH)								,
	charge(lrd_rhs.charge)						,
	TFD(lrd_rhs.TFD)							,
	rd_names( move(lrd_rhs.rd_names) )			,
	lrds( move(lrd_rhs.lrds) )					{
}
/***********************************************************************************/
local_rd& local_rd::operator=(local_rd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs ){
		name			= move(lrd_rhs.name);
		FD				= lrd_rhs.FD;
		charge			= lrd_rhs.charge;
		LH				= lrd_rhs.LH;
		rd_names		= move(lrd_rhs.rd_names);
		lrds			= move(lrd_rhs.lrds);
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
	lrds[14] = lrds[5]*grd.grds[5];
	lrds[14] = lrds[5]*grd.grds[5] - (lrds[6]*grd.grds[6]);
	lrds[15] = lrds[8]*grd.grds[9];
	lrds[16] = lrds[7]*grd.grds[9];
	lrds[17] = lrds[8]*(grd.grds[9]*grd.grds[9]);
	lrds[18] = lrds[8]*grd.grds[10];
}
/***********************************************************************************/
void local_rd::calculate_Fukui_potential(){
	
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
	//--------------------------------
	double initi_time = omp_get_wtime();
	m_log->input_message("Time for Calculate Fukui Potential: ");
	//--------------------------------
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
								elec_H[i*g1*g1+j*g2+k]	+= lrds[5].scalar[x*g1*g1+y*g2+z]/r;
								nuc_H[i*g1*g1+j*g2+k]	+= lrds[6].scalar[x*g1*g1+y*g2+z]/r;
								rad_H[i*g1*g1+j*g2+k]	+= lrds[7].scalar[x*g1*g1+y*g2+z]/r;
							}
						}
					}
				}
			}
		}
	}
	}
	double fin_time = omp_get_wtime() - initi_time;
	m_log->input_message(fin_time);
	//---------------------------------------------
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
	LH = true;
	//Local Chemical Potential method
	double numofelec	= lrds[2].calc_cube_integral();
	double val1			= grd.grds[7]/numofelec;
	double val2			= grd.grds[8]*2;
		
	Icube elec_dens_norm= lrds[2]/ numofelec ;
	Icube elec_dens_hard= elec_dens_norm*val2;
	Icube homo_elec_dens= lrds[5] - elec_dens_norm;
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
	//--------------------------------
	double initi_time = omp_get_wtime();
	m_log->input_message("Time for Calculate Local Hardness: ");
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
	double fin_time = omp_get_wtime() - initi_time;
	m_log->input_message(fin_time);
	double volume = std::abs(s1*s2*s3);
	lrds[10] = lrds[2];
	for(i=0;i<elec_H.size();i++) { lrds[10].scalar[i] = elec_H[i]; }
	lrds[10] = lrds[10]*volume;
	lrds[10] = lrds[10]*(1/numofelec);
	//----------------------------------------------------------------------------
	
	if ( TFD ){
		Icube density_rc1	= lrds[2].scale_cube(0.6666667);
		Icube density_rc2	= lrds[2].scale_cube(0.3333333);
		double Ck	= 2.8172;
		double Cx	= 0.7386;
		double con1	= 10.0/(9.0*numofelec);
		double con2	= 4/(9.0*numofelec);
		lrds[21] = density_rc1*(Ck*con1);
		lrds[21] = lrds[21]+lrds[10];
		density_rc2 = density_rc2*(Cx*con2);
		lrds[21] = lrds[21] - density_rc2;
	}
}
/***********************************************************************************/
void local_rd::calculate_MEP(const Imolecule& mol){
	
	vector<double> MEP(lrds[5].voxelN);
	unsigned int x,y,z;
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

	#pragma omp parallel default(shared) private(x,y,z,r,xx,yy,zz,ii,jj,kk)
	{
	#pragma omp parallel for reduction(vec_d_plus:MEP)
	for (x=0;x<g1;x++){
		for (y=0;y<g2;y++){
			for (z=0;z<g3;z++){
				xx	= x*s1 + o1;
				yy	= y*s2 + o2;
				zz	= z*s3 + o3;
				for (int na=0; na<mol.atoms.size(); na++ ){
					ii	= mol.atoms[na].xcoord - xx;
					jj	= mol.atoms[na].ycoord - yy;
					kk	= mol.atoms[na].zcoord - zz;
					ii *= ii;
					jj *= jj;
					kk *= kk;
					r = sqrt(ii+jj+kk);
					if ( r == 0.000 ){
						MEP[x*g1*g1+y*g2+z] += 0;
					}else{
						MEP[x*g1*g1+y*g2+z] += mol.atoms[na].charge /r;
						
					}
				}
			}
		}
	}
	}
	
	lrds[19] = lrds[2];
	for(unsigned i=0;i<MEP.size();i++) { lrds[19].scalar[i] = MEP[i]; }
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
	
	if ( FD )  {
		typestr = "Descriptor calculated with Finite Differences approximation \n";
		typestr2 = "FD_";
	}else{
		typestr = "Descriptor calculated with Frozen Orbital Approximation \n";
		typestr2 = "FOA_";
	}
	
	string name_type = name + typestr2;
	std::vector<string>	cube_names;

	cube_names.push_back(name_type+rd_names[0]+"_ph1");	// 0
	cube_names.push_back(name_type+rd_names[0]+"_ph2");	// 1
	cube_names.push_back(name_type+rd_names[1]+"_ph1");	// 2
	cube_names.push_back(name_type+rd_names[1]+"_ph2");	// 3
	cube_names.push_back(name_type+rd_names[2]);		// 4
	cube_names.push_back(name_type+rd_names[3]);		// 5
	cube_names.push_back(name_type+rd_names[4]);		// 6
	cube_names.push_back(name_type+rd_names[5]);		// 7
	cube_names.push_back(name_type+rd_names[6]);		// 8
	cube_names.push_back(name_type+rd_names[7]);		// 9
	cube_names.push_back(name_type+rd_names[8]+"_ph1");	// 10
	cube_names.push_back(name_type+rd_names[8]+"_ph2");	// 11
	cube_names.push_back(name_type+rd_names[9]);		// 12
	cube_names.push_back(name_type+rd_names[10]);		// 13
	cube_names.push_back(name_type+rd_names[11]);		// 14
	cube_names.push_back(name_type+rd_names[12]);		// 15
	cube_names.push_back(name_type+rd_names[13]);		// 16
	cube_names.push_back(name_type+rd_names[14]);		// 17
	cube_names.push_back(name_type+rd_names[15]+"_ph1");// 18
	cube_names.push_back(name_type+rd_names[15]+"_ph2");// 19
	cube_names.push_back(name_type+rd_names[16]);		// 20
	cube_names.push_back(name_type+rd_names[17]);		// 21
	cube_names.push_back(name_type+rd_names[18]+"_ph1");// 22
	cube_names.push_back(name_type+rd_names[18]+"_ph2");// 23
	cube_names.push_back(name_type+rd_names[19]+"_ph1");// 24
	cube_names.push_back(name_type+rd_names[19]+"_ph2");// 25
	cube_names.push_back(name_type+rd_names[20]);		// 26
	cube_names.push_back(name_type+rd_names[21]);		// 27
	
	lrds[0].header	= "Highest energy Occupied Molecular Orbital\n" + typestr;
	lrds[1].header	= "Lowest energy Unnocupied Molecular Orbital\n"+ typestr;
	lrds[2].header	= "Total electron density calculated with PRIMoRDiA\n"+ typestr;
	lrds[3].header	= "Total electron density calculated with PRIMoRDiA\n"+ typestr;
	lrds[4].header	= "Total electron density calculated with PRIMoRDiA\n"+ typestr;
	lrds[5].header	= "Left Fukui function, electrophilic attack succescitibily or local electofilicity\n"	+ typestr;
	lrds[6].header	= "Right Fukui Function, nucleophilic attack succescitibily or local nucleofilicity\n"	+ typestr;
	lrds[7].header	= "Average Fukui Function, radical attack succescitibily\n"								+ typestr;
	lrds[8].header	= "Dual deacriptor definition from Fukui function, or netfilicity\n"					+ typestr;
	lrds[9].header	= "Local hardness, calculated with working equation based on local chemical potential definition\n" + typestr;
	lrds[10].header	= "Local hardness, calculated with working equation based on eletron-eletron potential \n" + typestr;
	lrds[11].header	= "Fukui potential, using the left Fukui function to approximate the electron density\n" + typestr;
	lrds[12].header	= "Fukui potential, using the right Fukui function to approximate the electron density\n"					+ typestr;
	lrds[13].header	= "Fukui potential, using the zero Fukui function to approximate the electron density\n"					+ typestr;
	lrds[14].header	= "Local hardness, using the Fukui function to distribute the global hardness\n" + typestr;
	lrds[15].header	= "Local softness, using the the net Fukui function to distribute the global sofntess\n"+ typestr;
	lrds[16].header	= "Local softness, using the the average Fukui function to distribute the global sofntess\n"+ typestr;
	lrds[17].header	= "Hyper local softness\n" + typestr;
	lrds[18].header	= "Multifilicity\n" + typestr;	
	lrds[19].header	= "Molecular Electrostatic Potential\n" + typestr;
	lrds[20].header	= "Molecular Orbitals band localization\n" + typestr;
	lrds[21].header	= "Local hardness based on the Thomas-Fermi-Dirac functionals definitions\n" + typestr;
	
	for(unsigned i=0;i<lrds.size();i++){
		lrds[i].name = rd_names[i];
	}
		
	if ( FD ){
		lrds[2].write_cube(cube_names[2]+".cube"); //electron density
		lrds[3].write_cube(cube_names[5]+".cube"); //electron density cation
		lrds[4].write_cube(cube_names[6]+".cube"); //electron density anion
	}else{
		lrds[0].write_cube(cube_names[0]+".cube"); // HOMO ph1
		lrds[0].write_cube(cube_names[1]+".cube"); // HOMO ph2
		lrds[1].write_cube(cube_names[2]+".cube"); // LUMO ph1
		lrds[1].write_cube(cube_names[3]+".cube"); // LUMO ph2
	}
	
	//Default ouutput
	lrds[5].write_cube(cube_names[7]+".cube"); // left Fukui
	lrds[6].write_cube(cube_names[8]+".cube"); // right Fukui
	lrds[7].write_cube(cube_names[9]+".cube"); // zero Fukui
	lrds[8].write_cube(cube_names[10]+".cube"); // dual Fukui ph1
	lrds[8].write_cube(cube_names[11]+".cube"); // dual Fukui ph2
	lrds[11].write_cube(cube_names[14]+".cube"); // left Fukui potential 
	lrds[14].write_cube(cube_names[17]+".cube"); // local hardnes int
	lrds[15].write_cube(cube_names[18]+".cube"); // local softness dual ph1
	lrds[15].write_cube(cube_names[19]+".cube"); // local softness dual ph2


	if ( LH ){
		lrds[9].write_cube(cube_names[12]+".cube"); //local hardness lcp
		lrds[10].write_cube(cube_names[13]+".cube"); // local hardness Vee
		lrds[19].write_cube(cube_names[24]+".cube"); // MEP ph1
		lrds[19].write_cube(cube_names[25]+".cube"); // MEP ph2
		if ( TFD ){
			lrds[21].write_cube(cube_names[27]+".cube"); // local hardness TFD complete functional
		}
		if ( !FD )  lrds[2].write_cube(cube_names[5]+".cube"); // total electron density
	}
	
	if ( extra_RD ){
		lrds[12].write_cube(cube_names[15]+".cube"); // right Fukui potential
		lrds[13].write_cube(cube_names[16]+".cube"); // zero Fukui potential
		lrds[16].write_cube(cube_names[20]+".cube"); // local sofntess average
		lrds[17].write_cube(cube_names[21]+".cube"); // hyper local softness
		lrds[18].write_cube(cube_names[22]+".cube"); // multiphilicity ph1
		lrds[18].write_cube(cube_names[23]+".cube"); // multiphilicity ph2
	}
	if ( band ){
		lrds[20].write_cube(cube_names[26]+".cube"); // MO band localization
	}

	m_log->input_message("Finishing writting the local reactivity descriptos grids.\n");
}
/***********************************************************************************/
local_rd::~local_rd(){

}
//================================================================================
//END OF FILE
//================================================================================

