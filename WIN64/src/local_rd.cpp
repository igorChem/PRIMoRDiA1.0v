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
}
/***********************************************************************************/
local_rd::local_rd(const Icube& HOmo		,
					const Icube& LUmo)		:
	name(HOmo.name)							,
	finite_diff(false)						,
	locHardness(false)						,
	homo(HOmo)								,
	lumo(LUmo)								{
}
/***********************************************************************************/
local_rd::local_rd(const Icube& elec_dens	,
				   const Icube& HOmo		,
				   const Icube& LUmo		):
	name(HOmo.name)							,
	finite_diff(false)						,
	locHardness(true)						,
	charge(1)								,
	homo(HOmo)								,
	lumo(LUmo)								,
	elec_dens(elec_dens)					{
}
/***********************************************************************************/
local_rd::local_rd(const Icube& elecDens		,
						const Icube& cationDens	, 
						const Icube& anionDens	,
						int chg					):
	name(elecDens.name)							,
	finite_diff(true)							,
	locHardness(true)							,
	charge(chg)									,
	elec_dens(elecDens)							,
	cation(cationDens)							,
	anion(anionDens)							{

	homo	= elec_dens	- cation;
	lumo	= anion		- elec_dens;
}
/***********************************************************************************/
local_rd::local_rd(const local_rd& lrd_rhs)	:
	name(lrd_rhs.name)						,
	finite_diff(lrd_rhs.finite_diff)		,
	locHardness(lrd_rhs.locHardness)		,
	charge(lrd_rhs.charge)					,
	elec_dens(lrd_rhs.elec_dens)			,
	cation(lrd_rhs.cation)					,
	anion(lrd_rhs.anion)					,
	homo(lrd_rhs.homo)						,
	lumo(lrd_rhs.lumo)						,
	EAS(lrd_rhs.EAS)						,
	NAS(lrd_rhs.NAS)						,
	RAS(lrd_rhs.RAS)						,
	Dual(lrd_rhs.Dual)						,
	Hardness(lrd_rhs.Hardness)				,
	Softness_Dual(lrd_rhs.Softness_Dual)	,
	Hyper_Softness(lrd_rhs.Hyper_Softness)	,
	multifilic(lrd_rhs.multifilic)			{
}
/***********************************************************************************/
local_rd& local_rd::operator=(const local_rd& lrd_rhs){
	if( this!=&lrd_rhs ){
		name			= lrd_rhs.name;
		finite_diff		= lrd_rhs.finite_diff;
		charge			= lrd_rhs.charge;
		locHardness		= lrd_rhs.locHardness;
		elec_dens		= lrd_rhs.elec_dens;
		cation			= lrd_rhs.cation;
		anion			= lrd_rhs.anion;
		homo			= lrd_rhs.homo;
		lumo			= lrd_rhs.lumo;	
		EAS				= lrd_rhs.EAS;
		NAS				= lrd_rhs.NAS;
		RAS				= lrd_rhs.RAS;
		Dual			= lrd_rhs.Dual;
		Hardness		= lrd_rhs.Hardness;
		Softness_Dual	= lrd_rhs.Softness_Dual;
		Hyper_Softness	= lrd_rhs.Hyper_Softness;
		multifilic		= lrd_rhs.multifilic;
	}
	return *this;
}
/***********************************************************************************/
local_rd::local_rd(local_rd&& lrd_rhs) noexcept	:
	name(lrd_rhs.name)							,
	finite_diff(lrd_rhs.finite_diff)			,
	locHardness(lrd_rhs.locHardness)			,
	charge(lrd_rhs.charge)						,
	elec_dens( move(lrd_rhs.elec_dens) )		,
	cation( move(lrd_rhs.cation) )				,
	anion( move(lrd_rhs.anion) )				,
	homo( move(lrd_rhs.homo) )					,
	lumo( move(lrd_rhs.lumo) )					,
	EAS( move(lrd_rhs.EAS) )					,
	NAS( move(lrd_rhs.NAS) )					,
	RAS( move(lrd_rhs.RAS) )					, 
	Dual( move(lrd_rhs.Dual) )					,
	Hardness( move(lrd_rhs.Hardness) )			,
	Softness_Dual( move( lrd_rhs.Softness_Dual) ),
	Hyper_Softness( move( lrd_rhs.Hyper_Softness) ),
	multifilic( move(lrd_rhs.multifilic) )		{
}
/***********************************************************************************/
local_rd& local_rd::operator=(local_rd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs ){
		name			= move(lrd_rhs.name);
		finite_diff		= lrd_rhs.finite_diff;
		charge			= lrd_rhs.charge;
		locHardness		= lrd_rhs.locHardness;
		elec_dens		= move(lrd_rhs.elec_dens);
		cation			= move(lrd_rhs.cation);
		anion			= move(lrd_rhs.anion);
		homo			= move(lrd_rhs.homo);
		lumo			= move(lrd_rhs.lumo);
		EAS				= move(lrd_rhs.EAS);
		NAS				= move(lrd_rhs.NAS);
		RAS				= move(lrd_rhs.RAS);
		Dual			= move(lrd_rhs.Dual);
		Hardness		= move(lrd_rhs.Hardness);
		Softness_Dual	= move(lrd_rhs.Softness_Dual);
		Hyper_Softness	= move(lrd_rhs.Hyper_Softness);
		multifilic		= move(lrd_rhs.multifilic);
	}
	return *this;
}
/***********************************************************************************/
void local_rd::calculate_fukui(){
	if ( !finite_diff ) {
		EAS	= homo;
		NAS	= lumo;
		RAS	= (homo	+	lumo)/2.0;
		Dual= lumo - homo;
	}else{
		EAS	= homo/charge;
		NAS	= lumo/charge;
		RAS	= (EAS	+ NAS)/2.0;
		Dual = NAS	- EAS; 
	}
}
/***********************************************************************************/
void local_rd::calculate_RD(const global_rd& grd){
	Softness_Dual	= Dual*grd.softness;
	Hyper_Softness	= RAS*grd.softness;
	multifilic		= Dual*grd.Electrophilicity;
}
/***********************************************************************************/
void local_rd::calculate_hardness(const global_rd& grd, string method){
	locHardness = true;
	if ( method == "LCP" ){
		double numofelec	= elec_dens.calc_cube_integral();
		double val1			= grd.chemical_pot/numofelec;
		double val2			= grd.hardness*2; 	
	
		Icube elec_dens_norm= elec_dens / numofelec ;
		Icube elec_dens_hard= elec_dens_norm*val2;
		Icube homo_elec_dens= homo.normalize() - elec_dens_norm;
		homo_elec_dens		= homo_elec_dens*val1;
		Hardness			= homo_elec_dens + elec_dens_hard;
	}
	else if ( method == "mepFukui" ){
		vector<double> elec_H(homo.voxelN);
		vector<double> nuc_H(homo.voxelN);
		unsigned int i,j,k,x,y,z;
		double ii	= 0;
		double jj	= 0;
		double kk	= 0;
		double xx	= 0;
		double yy	= 0;
		double zz	= 0;
		double r	= 0;
		double s1	= homo.gridsides[0];
		double s2	= homo.gridsides[1];
		double s3	= homo.gridsides[2];
		double o1	= homo.origin[0];
		double o2	= homo.origin[1];
		double o3	= homo.origin[2];
		unsigned int g1 = homo.grid[0];
		unsigned int g2 = homo.grid[1];
		unsigned int g3 = homo.grid[2];
		
		#pragma omp declare reduction(vec_d_plus : std::vector<double> : \
					std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)

		#pragma omp parallel default(shared) private(i,j,k,x,y,z,r,xx,yy,zz,ii,jj,kk)
		{
		#pragma omp parallel for reduction(vec_d_plus:elec_H,nuc_H)
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
								r  		= sqrt(xx + yy + zz);
								if ( r == 0.000 ){
									elec_H[i*g1*g1+j*g2+k]	+= 0;
									nuc_H[i*g1*g1+j*g2+k] 	+= 0;
								}else{
									elec_H[i*g1*g1+j*g2+k] += homo.scalar[x*g1*g1+y*g2+z]/r;
								}
							}
						}
					}
				}
			}
		}
		
		}
		double volume = std::abs(s1*s2*s3);
		Hardness = homo;
		for(i=0;i<elec_H.size();i++) { Hardness.scalar[i] = elec_H[i]; }
		Hardness = Hardness * volume;
	}
	else if ( method == "mepEE" ){
		vector<double> elec_H(homo.voxelN);
		unsigned int i,j,k,x,y,z;
		double ii	= 0;
		double jj	= 0;
		double kk	= 0;
		double xx	= 0;
		double yy	= 0;
		double zz	= 0;
		double r	= 0;
		double s1	= homo.gridsides[0];
		double s2	= homo.gridsides[1];
		double s3	= homo.gridsides[2];
		double o1	= homo.origin[0];
		double o2	= homo.origin[1];
		double o3	= homo.origin[2];
		unsigned int g1 = homo.grid[0];
		unsigned int g2 = homo.grid[1];
		unsigned int g3 = homo.grid[2];
		
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
									elec_H[i*g1*g1+j*g2+k] += elec_dens.scalar[x*g1*g1+y*g2+z]/r;
								}
							}
						}
					}
				}
			}
		}
		}
		double volume = std::abs(s1*s2*s3);
		Hardness = homo;
		for(i=0;i<elec_H.size();i++) { Hardness.scalar[i] = elec_H[i]; }
		Hardness	= Hardness * volume;
		int norm	= 2*elec_dens.calc_cube_integral();
		Hardness	= Hardness/norm;
	}	
	else if ( method == "fukui" ){
		Hardness = lumo*grd.lumo_en - homo*grd.homo_en;
	}
	else{
		locHardness = false;
	}
}
/***********************************************************************************/
local_rd operator-(const local_rd& lrd_lhs,const local_rd& lrd_rhs){
	local_rd Result(lrd_lhs);
	if( lrd_lhs.homo == lrd_rhs.homo ) {
		Result.EAS				= lrd_lhs.EAS			-	lrd_rhs.EAS;
		Result.NAS				= lrd_lhs.NAS			-	lrd_rhs.NAS;
		Result.RAS				= lrd_lhs.RAS			-	lrd_rhs.RAS;
		Result.Dual				= lrd_lhs.Dual			-	lrd_rhs.Dual;
		Result.Hyper_Softness 	= lrd_lhs.Hyper_Softness-	lrd_rhs.Hyper_Softness;
		Result.Softness_Dual	= lrd_lhs.Softness_Dual	-	lrd_rhs.Softness_Dual;
		Result.multifilic		= lrd_lhs.multifilic	-	lrd_rhs.multifilic;
		if ( lrd_lhs.locHardness ) Result.Hardness = lrd_lhs.Hardness - lrd_rhs.Hardness;
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
	
	string fukui_suc_elec	=	name + typestr2 + "_EAS"; 
	string fukui_suc_nuc	=	name + typestr2 + "_NAS"; 
	string fukui_suc_rad	=	name + typestr2 + "_RAS"; 
	string deltaFukui		=	name + typestr2 + "_dual";
	string deltaFukui2		=	name + typestr2 + "_dual2";
	string local_hardness	=	name + typestr2 + "_Hardness";
	EAS.header				= "Left Fukui function, electrophilic attack succescitibily\n"		+ typestr;
	NAS.header				= "Right Fukui Function, nucleophilic attack succescitibily\n"		+ typestr;
	RAS.header				= "Average Fukui Function, radical attack succescitibily\n"			+ typestr;
	Dual.header				= "Net Fukui Function, dual descriptor of attack succescitibility\n"+ typestr;
	Hardness.header			= "Local Hardness \n" + typestr;
	
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

