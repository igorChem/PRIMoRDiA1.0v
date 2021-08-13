// source file for the local_rd class 
// local_rd.cpp
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
#include <iostream>
#include <string> 
#include <vector>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <experimental/filesystem>
// include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/QMdriver.h"
#include "../include/global_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/Iprotein.h"
#include "../include/residue_lrd.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::move;
using std::abs;
using std::unique_ptr;
namespace fs = std::experimental::filesystem;

const double precision = 1e-08;


std::vector<std::string> rd_names = {
	"nucleophilicity"				, //0
	"electrophilicity"				, //1
	"radicality"					, //2
	"netphilicity"					, //3
	"hardness_Vee"					, //4
	"hardness_lcp"					, //5
	"Fukui_pot_left"				, //6
	"Fukui_pot_right"				, //7
	"Fukui_pot_zero"				, //8
	"softness_dual"					, //9
	"hyper_softness"				, //10
	"multiphilicity"				, //11
	"fukushima"						, //12
	"charge"						, //13
	"electron_density"				, //14
	"mep"							, //15
	"hardness_TFD"					, //16
	"softness_avg"					, //17
	"hardness_int"					  //18
};

/***********************************************************************************/
local_rd_cnd::local_rd_cnd()	:
	name("nonamed")				,
	FD(false)					,
	charge(0)					{
	
	names = rd_names;
	lrds.resize( rd_names.size() );
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(unsigned int nof):
	FD(false)								,
	names(rd_names)							,
	charge(1)								{
	
	lrds.resize( rd_names.size() );
	
	for(unsigned int i=0;i<lrds.size(); i++){
		lrds.resize( nof );
	}
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const Imolecule& mol_neut	,
							const Imolecule& mol_cation	, 
							const Imolecule& mol_anion)	:
	name(mol_neut.name)									, 
	FD(true)											,
	names(rd_names)										,
	charge(1)											{
	
	lrds.resize( rd_names.size() );
	
	for(unsigned int i=0;i<lrds.size(); i++){
		lrds[i].resize( mol_neut.atoms.size() );
	}
	
	for (unsigned int i=0;i< mol_neut.atoms.size();i++){
		lrds[0][i]	= ( (-mol_neut.atoms[i].charge)		- (-mol_cation.atoms[i].charge) )/charge;
		lrds[1][i]	= ( (-mol_anion.atoms[i].charge)	- (-mol_neut.atoms[i].charge) ) /charge;
		lrds[2][i]	= ( ( (-mol_anion.atoms[i].charge)	- (-mol_cation.atoms[i].charge) )/2 ) /charge;
		lrds[3][i]	= lrds[1][i] - lrds[0][i];
		lrds[13][i]	= mol_neut.atoms[i].charge;
	}
}         
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const local_rd_cnd& lrd_rhs)	:
	name(lrd_rhs.name)									,
	FD(lrd_rhs.FD)										,
	names(lrd_rhs.names)								,
	charge(lrd_rhs.charge)								,
	lrds(lrd_rhs.lrds)									{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(const local_rd_cnd& lrd_rhs){
	if( this!=&lrd_rhs){
		name	= lrd_rhs.name;
		FD		= lrd_rhs.FD;
		names	= lrd_rhs.names;
		charge	= lrd_rhs.charge;
		lrds	= lrd_rhs.lrds;
	}
	return *this;
}                
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(local_rd_cnd&& lrd_rhs) noexcept	:
	name( move(lrd_rhs.name) )								,
	FD( move(lrd_rhs.FD) ) 									,
	names( move(lrd_rhs.names) )							,
	charge( move(lrd_rhs.charge) )							,
	lrds( move(lrd_rhs.lrds) )								{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(local_rd_cnd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs){
		name	= move(lrd_rhs.name);
		FD		= move(lrd_rhs.FD);
		names	= move(lrd_rhs.names);
		charge	= move(lrd_rhs.charge);
		lrds	= move(lrd_rhs.lrds);
	}       	
	return *this;
}
/***********************************************************************************/
local_rd_cnd operator-(const local_rd_cnd& lrd_lhs,const local_rd_cnd& lrd_rhs){
	local_rd_cnd Result = lrd_lhs;
	for (unsigned int i=0;i<lrd_lhs.lrds.size();i++){
		for (unsigned int j=0;j<lrd_lhs.lrds[i].size();j++){
			Result.lrds[i][j] = lrd_lhs.lrds[i][j] - lrd_rhs.lrds[i][j];
		}
	}
	return Result;
}
/***********************************************************************************/
void local_rd_cnd::calculate_frontier_orbitals( const Imolecule& molecule, unsigned band){
	name					= molecule.name;
	double value_h			= 0.0;
	double value_l			= 0.0;
	unsigned init_orb		= 0;
	unsigned n_aorbs		= 0;
	unsigned ao 			= molecule.num_of_ao;
	unsigned homon 			= abs(molecule.homoN);
	unsigned lumon 			= abs(molecule.lumoN);
	unsigned cnt 			= 0;
	unsigned init			= homon-band;
	unsigned fin 			= lumon+band;
	
	for( unsigned atom=0; atom<molecule.atoms.size(); atom++ ){
		//defining the indices of the atomic orbitals
		if ( atom == 0 ) {
			init_orb = 0;
			n_aorbs  = molecule.atoms[0].orbitals.size();
		}else{
			for( unsigned j=0; j<atom; j++ ){
				init_orb += molecule.atoms[j].norb;
			}
			for( unsigned j=0; j<=atom; j++ ){
				n_aorbs  += molecule.atoms[j].norb;
			}
		}
		//-------------------------------------------------
		
		//calculating the occupied molecular orbitals
		unsigned init = homon-band;
		for( unsigned i=init; i<=homon; i++ ){
			for( unsigned mu=init_orb; mu<n_aorbs; mu++ ){
				for ( unsigned nu=init_orb; nu<n_aorbs; nu++ ) {
					value_h+=molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + nu]*
							molecule.m_overlap[nu+(mu*(mu+1) )/2];
				}
			}
		}
		
		if ( atom == 0 && band > 0 ){
			m_log->input_message("Number of occupied MO used to calculate condensed to atom descriptors: ");
			m_log->input_message( int(band) );
			m_log->input_message("\n");
		}		
		//-------------------------------------------------
		
		//calculating the virtual molecular orbitals
		
		for( unsigned i=lumon; i<=fin; i++ ){
			for( unsigned mu=init_orb; mu<n_aorbs; mu++ ){
				for ( unsigned nu=init_orb; nu<n_aorbs; nu++ ) {
					value_l +=molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + nu]*
							molecule.m_overlap[nu+(mu*(mu+1))/2];
				}
			}
		}
		if ( atom == 0  && band > 0 ){
			m_log->input_message("Number of  virtual MO used to calculate condensed to atom descriptors: ");
			m_log->input_message( int(band) ) ;
			m_log->input_message("\n");
		}
		
		//--------------------------------------------------
		
		// calculting the beta orbitals
		if ( molecule.betad ){
			//calculating the occupied molecular orbitals
			for( unsigned i=init; i<=homon; i++ ){
				for( unsigned mu=init_orb; mu<n_aorbs; mu++){
					for ( unsigned nu=init_orb; nu<n_aorbs; nu++ ) {
						value_h+=molecule.coeff_MO_beta[ao*i + mu]*
								molecule.coeff_MO_beta[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1) )/2];
					}
				}
			}
			if ( atom == 0 && band > 0 ){
				m_log->input_message("Number of occcupied virtual MO used to calculate condensed to atom descriptors: ");
				m_log->input_message( int(band) );
				m_log->input_message("\n");
			}
			value_h /= 2;
			//------------------------------------------------
			
			//calculating the virtual molecular orbitals
			for( unsigned i=lumon; i<=fin; i++ ){
				if ( molecule.orb_energies_beta[i] <= molecule.lumo_energy+energy_crit ){
					for(unsigned mu=init_orb; mu<n_aorbs; mu++ ){
						for (unsigned nu=init_orb; nu<n_aorbs; nu++ ) {
							value_l +=molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
						}
					}
				}
			}
			if ( atom == 0  && band > 0 ){
				m_log->input_message("Number of  virtual beta MO used to calculate condensed to atom descriptors: ");
				m_log->input_message( int(band) );
				m_log->input_message("\n");
			}
			value_l /= 2;
		}
		lrds[0][atom] = value_h;
		lrds[1][atom] = value_l;
		lrds[13][atom]= molecule.atoms[atom].charge;
		lrds[12][atom]= value_h + value_l;
	}
	if ( band > 1 ){
		lrds[0] = norm_dvec(lrds[0],5);
		lrds[1] = norm_dvec(lrds[1],5);
	}else{
		lrds[0] = norm_dvec(lrds[0],1);
		lrds[1] = norm_dvec(lrds[1],1);
	}
	for( unsigned i=0; i<lrds[0].size(); i++ ){
		lrds[3][i] = lrds[1][i] - lrds[0][i];
		lrds[2][i] = (lrds[1][i] + lrds[0][i])/2;
	}	

}
/***********************************************************************************/
void local_rd_cnd::energy_weighted_fukui_functions( const Imolecule& molecule ){
	double pre_coef		= exp( -abs(energy_crit) );
	double coefficient	= 0.0;
	int mo_count 		= 0;
	double value_h		= 0.0;
	double value_l		= 0.0;
	unsigned init_orb	= 0;
	unsigned n_aorbs	= 0;
	unsigned ao 		= molecule.num_of_ao;
	unsigned homon 		= abs(molecule.homoN);
	unsigned lumon 		= abs(molecule.lumoN);

	for( unsigned atom=0; atom<molecule.atoms.size(); atom++ ){
		// defining the atomic orbitals indices
		if ( atom == 0 ) {
			init_orb = 0;
			n_aorbs  = molecule.atoms[0].orbitals.size();
		}else{
			for( unsigned j=0; j<atom; j++ ){
				init_orb += molecule.atoms[j].norb;
			}
			for( unsigned j=0; j<=atom; j++ ){
				n_aorbs  += molecule.atoms[j].norb;
			}
		}		
		//------------------------------------------------
		
		//calculating the occupied molecular orbitals
		for( unsigned i=0; i<=homon; i++ ){
			coefficient = exp(-abs( molecule.orb_energies[i]-molecule.homo_energy ) );
			if ( coefficient > pre_coef ){ 
				mo_count++;
				for( unsigned mu=init_orb; mu<n_aorbs; mu++){
					for ( unsigned nu=init_orb; nu<n_aorbs; nu++ ) {
						value_h +=coefficient*
								molecule.coeff_MO[ao*i + mu]*
								molecule.coeff_MO[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
			}
		}
		
		if ( atom == 0 ){
			m_log->input_message("Number of occupied MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(mo_count);m_log->input_message("\n");	
		}
		
		//------------------------------------------------
			
		//calculating the virtual molecular orbitals
		for( unsigned i=lumon;i<molecule.orb_energies.size();i++){
			coefficient = exp(-abs( molecule.orb_energies[i]-molecule.lumo_energy ) );
			if ( coefficient > pre_coef ){ 
				mo_count++;
				for(unsigned mu=init_orb;mu<n_aorbs;mu++){
					for (unsigned nu=init_orb;nu<n_aorbs;nu++) {
						value_l +=coefficient*
								molecule.coeff_MO[ao*i + mu]*
								molecule.coeff_MO[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
			}
		}
		if ( atom == 0 ){
			m_log->input_message("Number of  virtual  MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(mo_count);
			m_log->input_message("\n");	
		}		
		
		//------------------------------------------------
			
		//calculating the beta molecular orbitals
		if ( molecule.betad ){
			mo_count = 0;
			for( int i=0;i<=homon;i++){
				coefficient = exp(-abs(molecule.orb_energies_beta[i]-molecule.homo_energy ) );
				if ( coefficient > pre_coef ){ 
					mo_count++;
					for(int mu=init_orb;mu<n_aorbs;mu++){
						for (int nu=init_orb;nu<n_aorbs;nu++) {
							value_h +=coefficient*
									molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
						}
					}
				}
			}
			value_h /=2;
			
			if ( atom == 0 ){
				m_log->input_message("Number of occupied beta MO used to calculate condensed to atom descriptors: ");
				m_log->input_message(mo_count);
				m_log->input_message("\n");			
			}
			
			for( int i=lumon;i<molecule.orb_energies.size();i++){
				coefficient = exp(-abs(molecule.orb_energies_beta[i]-molecule.lumo_energy ) );
				if ( coefficient >  pre_coef  ){ 
					mo_count++;
					for(int mu=init_orb;mu<n_aorbs;mu++){
						for (int nu=init_orb;nu<n_aorbs;nu++) {
							value_l +=coefficient*
									molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
						}
					}
				}
			}
			value_l /=2;
			
			if ( atom == 0 ){
				m_log->input_message("Number of  virtual  MO used to calculate condensed to atom descriptors: ");
				m_log->input_message(mo_count);
				m_log->input_message("\n");	
			}
		}
		lrds[0][atom] = value_h;
		lrds[1][atom] = value_l;
		lrds[12][atom]= value_h + value_l;
		lrds[13][atom]	= molecule.atoms[atom].charge;
	}
	lrds[0] = norm_dvec(lrds[0],5);
	lrds[1] = norm_dvec(lrds[1],5);
	for( unsigned i=0; i<lrds[0].size(); i++ ){
		lrds[3][i] = lrds[1][i] - lrds[0][i];
		lrds[2][i] = (lrds[1][i] + lrds[0][i])/2;
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_fukui_potential( const Imolecule& molecule ){
	double r, xi, yi, zi = 0;
	unsigned nof = molecule.atoms.size();
	for (unsigned i=0; i<nof; i++ ){
		for (unsigned j=0; j<nof; j++ ){
			if ( i != j ){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				lrds[6][i] += lrds[0][j]/r;
				lrds[7][i] += lrds[1][j]/r;
				lrds[8][i] += lrds[2][j]/r;
			}
		}
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_RD(const global_rd& grd){
	for( unsigned j=0; j<lrds[0].size(); j++ ){
		lrds[9][j]	= lrds[3][j]*grd.grds[9];
		lrds[10][j]	= lrds[2][j]*grd.grds[9]*grd.grds[9];
		lrds[11][j]	= lrds[3][j]*grd.grds[10];
		lrds[17][j]	= lrds[2][j]*grd.grds[9];
		lrds[18][j]	= lrds[0][j]*grd.grds[5] - lrds[1][j]*grd.grds[6];
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_hardness(const global_rd& grd, const Imolecule& molecule){
	unsigned nof		= molecule.atoms.size();
	unsigned init_orb	= 0;
	unsigned n_aorbs	= 0;
	double value		= 0.0;
	double occ 			= 2.0;
	unsigned ao 		= molecule.num_of_ao; 

	//-----------------------------------------------------
	//Estimating electron density
	for( unsigned atom=0; atom<nof; atom++ ){
		if ( atom == 0 ) {
			init_orb = 0;
			n_aorbs  = molecule.atoms[0].orbitals.size();
		}else{
			for(unsigned j=0; j<atom; j++ ){
				init_orb += molecule.atoms[j].norb;
			}
			for(unsigned j=0; j<=atom; j++ ){
				n_aorbs  += molecule.atoms[j].norb;
			}
		}
		if ( molecule.betad ) { occ = 1.0 ;}
		for( unsigned i=0; i<=molecule.homoN; i++ ){
			for(unsigned mu=init_orb; mu<n_aorbs; mu++ ){
				for (unsigned nu=init_orb; nu<n_aorbs; nu++ ){
					value+=occ*
							molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + nu];
				}
			}
		}
		if ( molecule.occupied_beta.size() > 0 ){
			for( unsigned i=0; i<molecule.homoN; i++ ){
				if ( molecule.occupied_beta[i] > 0 ){
					for(unsigned mu=init_orb; mu<n_aorbs; mu++ ){
						for (unsigned nu=init_orb; nu<n_aorbs; nu++ ){
							value +=molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu];
						}
					}
				}
			}
		}
		lrds[14][atom] = value;
	}
	
	double xi, yi, zi, r = 0.000;
	//-----------------------------------------------------
	//calculating local hardness with method (electron-electron interaction)
	for (unsigned i=0; i<nof; i++ ){
		for (unsigned j=0; j<nof; j++ ){
			if ( i != j ){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				lrds[4][i] += lrds[14][j]/r;
			}
		}
		lrds[4][i] /= molecule.num_of_electrons;
	}
	
	//----------------------------------------------------
	//calculating local hardness with method (LCP)
	for(unsigned i=0; i<nof; i++ ){
		lrds[5][i] = lrds[0][i]*( grd.grds[7] /nof )
						+ (lrds[14][i]/nof)*grd.grds[8]*2;
	}
	
	//----------------------------------------------------
	//calculating local hardness with method (tomas-fermi-dirac statistical treatment)
	if (TFD){
		double Ck	= (3/10)*pow((3*M_PI*M_PI),2/3);
		double Cx	= (3/4*M_PI)*pow((3*M_PI*M_PI),1/3);
		std::vector<double> temp1 = lrds[14];
		std::vector<double> temp2 = lrds[14];
		std::vector<double> temp3 = lrds[14];
		std::vector<double> temp4 = lrds[14];
		
		for ( unsigned i=0; i<nof; i++ ){
			temp1[i]	= pow(temp1[i],0.33333333);
			lrds[16][i]	= (2/(9*molecule.num_of_electrons))*temp1[i];
			temp2[i]	= temp2[i]*5*Ck - 2*Cx;
			temp3[i]	= 0.458*temp1[i];
			temp4[i]	= pow((temp3[i]+1),3);
			temp4[i]	= -00466*( (temp3[i]+2) / temp4[i] );
			lrds[16][i]	= lrds[16][i]*(temp2[i] - temp4[i]);
			lrds[16][i] += lrds[4][i];
		}
	}
}
/*************************************************************************************/
void local_rd_cnd::calculate_mep(const Imolecule& molecule){
	double xi, yi, zi, r = 0.000;
	for( unsigned i=0; i<molecule.atoms.size(); i++ ){
		for( unsigned j=0; j<molecule.atoms.size(); j++ ){
			if (i != j){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				lrds[15][i] += molecule.atoms[j].atomicN/r;
				lrds[15][i] -= lrds[4][i];
			}
		}
	}
}
/*************************************************************************************/
protein_lrd local_rd_cnd::rd_protein(const Iprotein& prot){

	vector<residue_lrd> res_rd( prot.residues.size() );
	unsigned nof = lrds[0].size();
	unsigned cnt = 0;
	for( unsigned i=0; i<prot.residues.size(); i++ ){
		for( unsigned j=0; j<prot.residues[i].atom_s; j++ ){
			for( unsigned k=0; k<lrds.size(); k++ ){
				res_rd[i].rd_sum[k]	+= lrds[k][cnt++];
			}
		}
	}

	for( unsigned i=0; i<res_rd.size(); i++ ){
		for( unsigned j=0; j<res_rd[i].rd_sum.size(); j++ ){
			res_rd[i].rd_avg[j]	+= res_rd[i].rd_sum[j]/nof;
		}
	}
	
	protein_lrd protein_react_descriptors(res_rd);
	protein_react_descriptors.write_protein_lrd(prot);
	
	for( int i=0; i<prot.residues.size(); i++){
		protein_react_descriptors.labels.push_back( std::to_string(i+1) );
		protein_react_descriptors.labels[i] += prot.residues[i].type;
	}	
	return protein_react_descriptors;
}
/*******************************************************************************************/
void local_rd_cnd::write_rd_protein_pdb(const Iprotein& protein){
	
	Iprotein prot = protein;
	pdb rd_results;
	rd_results.name = get_file_name( protein.name.c_str() );
	prot.load_b_column(lrds[0]);
	prot.title	= rd_results.name + "_nucleophilicity";
	prot.remark	= "Electrophilic Attack Suscptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[1]);
	prot.title	= rd_results.name +"_electrophilicity";
	prot.remark	= "Nucleophilic Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[2]);
	prot.title	= rd_results.name +"_radicality";
	prot.remark	= "Radical Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[3]);
	prot.title	= rd_results.name +"_netphilicity";
	prot.remark	= "Local netphilicity Descriptor Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[4]);
	prot.title	= rd_results.name + "_hardness_lcp";
	prot.remark	= "Local Hardness (local chemical potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[5]);
	prot.title 	= rd_results.name +"_hardness_Vee";
	prot.remark	= "Local Hardness (Electron-electron potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[6]);
	prot.title	= rd_results.name +"_fukui_pot_left";
	prot.remark	= "Local Hardness (left-Fukui potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[7]);
	prot.title	= rd_results.name +"_fukui_pot_right";
	prot.remark	= "Local Hardness (right-Fukui potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[8]);
	prot.title	= rd_results.name +"_fukui_pot_zero";
	prot.remark	= "Local Hardness (zero-Fukui potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);	
	prot.load_b_column(lrds[9]);
	prot.title	= rd_results.name +"_softness_dual";
	prot.remark	= "Local Softness (dual descriptor distribution based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[10]);
	prot.title	= rd_results.name +"_hyper_softness";
	prot.remark	= "Local Hyper Softness at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[11]);
	prot.title	= rd_results.name +"_multiphilicity";
	prot.remark	= "Local Multiphilicity at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[12]);
	prot.title	= rd_results.name +"_fukushima";
	prot.remark	= "Localization of frontier band orbitals at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[13]);
	prot.title	= rd_results.name +"_charge";
	prot.remark	= "Local Multiphilic descriptor  at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[14]);
	prot.title	= rd_results.name +"_electron_density";
	prot.remark	= "Electron Density at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);	
	prot.load_b_column(lrds[15]);
	prot.title	= rd_results.name +"_mep";
	prot.remark	= "Molecular electrostatic potential per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[16]);
	prot.title	= rd_results.name +"_hardness_TFD";
	prot.remark	= "Local hardness with complete DFT functional per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[17]);
	prot.title	= rd_results.name +"_softness_avg";
	prot.remark	= "Local softness using zero Fukui function per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(lrds[18]);
	prot.title	= rd_results.name +"_hardness_int";
	prot.remark	= "Local hardness distributed by the Fukui function per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors in PDBs.\n");
	m_log->inp_delim(1);
	
	fs::create_directory(name+"_PDB_RD");
	rd_results.write_models(name+"_PDB_RD");
}
/*************************************************************************************/
void local_rd_cnd::write_LRD(){
	std::string temps;
	if ( FD ) { 
		temps = name+"FD.lrd";
	}else{
		temps = name+"FOA.lrd";
	}
	const char* Name_s = temps.c_str();
	std::ofstream lrd_file;
	lrd_file.open(Name_s);
	
	lrd_file.precision(6);
	lrd_file << std::fixed;
	lrd_file << name << " " << "\n" <<  std::left;
	lrd_file << "n atom ";
	
	for( unsigned i=0; i<names.size(); i++){
		lrd_file << names[i] << " ";
	}
	lrd_file << endl;
	
	for( unsigned i=0; i<lrds[0].size(); i++){
		lrd_file	<< (i+1)
					<< " ";
		for(unsigned j=0; j<lrds.size(); j++ ){
			lrd_file	<< lrds[j][i]
						<< " ";
		}
		lrd_file << endl;
	}
	lrd_file.close();
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors.\n");
	m_log->inp_delim(2);
}
//================================================================================
//END OF FILE
//================================================================================