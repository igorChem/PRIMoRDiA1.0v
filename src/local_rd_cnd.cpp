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

/***********************************************************************************/
local_rd_cnd::local_rd_cnd()	:
	name("nonamed")				,
	finite_diff(true)			,
	charge(0)					,
	ed(false)					,
	mep(false)					{
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(Imolecule&& mol) noexcept:
	name(mol.name)									,
	finite_diff(false)								,
	charge(1)										,
	ed(false)										,
	mep(false)										,
	molecule( move(mol) )							{
	
	EAS.resize(molecule.num_of_atoms);
	NAS.resize(molecule.num_of_atoms);
	RAS.resize(molecule.num_of_atoms);
	dual.resize(molecule.num_of_atoms);
	Softness_Dual.resize(molecule.num_of_atoms);
	Hyper_softness.resize(molecule.num_of_atoms);
	hardness_A.resize(molecule.num_of_atoms);
	hardness_B.resize(molecule.num_of_atoms);
	hardness_C.resize(molecule.num_of_atoms);
	hardness_D.resize(molecule.num_of_atoms);
	multifilic.resize(molecule.num_of_atoms);
	electrophilicity.resize(molecule.num_of_atoms);
	electron_density.resize(molecule.num_of_atoms);
	fukushima.resize(molecule.num_of_atoms);
	
}
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const Imolecule& mol_neut	,
							const Imolecule& mol_cation	, 
							const Imolecule& mol_anion)	:
	name(mol_neut.name)									, 
	finite_diff(true)									,
	charge(1)											,
	ed(false)											,
	mep(false)											,
	molecule( Imolecule() )								{
	
	molecule.atoms				= mol_neut.atoms;
	molecule.num_of_atoms		= mol_neut.num_of_atoms;
	molecule.num_of_electrons	= mol_neut.num_of_electrons;

	//molecule.update();
	
	EAS.resize(molecule.num_of_atoms);
	NAS.resize(molecule.num_of_atoms);
	RAS.resize(molecule.num_of_atoms);
	dual.resize(molecule.num_of_atoms);
	Softness_Dual.resize(molecule.num_of_atoms);
	Hyper_softness.resize(molecule.num_of_atoms);
	hardness_A.resize(molecule.num_of_atoms);
	hardness_B.resize(molecule.num_of_atoms);
	hardness_C.resize(molecule.num_of_atoms);
	hardness_D.resize(molecule.num_of_atoms);
	multifilic.resize(molecule.num_of_atoms);
	electrophilicity.resize(molecule.num_of_atoms);
	electron_density.resize(molecule.num_of_atoms);
	fukushima.resize(molecule.num_of_atoms);	
	
	for (unsigned int i=0;i<molecule.num_of_atoms;i++){
		EAS[i]				= ( (-molecule.atoms[i].charge)		- (-mol_cation.atoms[i].charge) )/charge;
		NAS[i]				= ( (-mol_anion.atoms[i].charge)	- (-molecule.atoms[i].charge) ) /charge;
		RAS[i]				= ( ( (-mol_anion.atoms[i].charge)	- (-mol_cation.atoms[i].charge) )/2 ) /charge;
		dual[i]				= NAS[i] - EAS[i];
	}
}         
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(const local_rd_cnd& lrd_rhs)	:
	name(lrd_rhs.name)									,
	ed(lrd_rhs.ed)										,
	mep(lrd_rhs.mep)									,
	finite_diff(lrd_rhs.finite_diff)					,
	charge(lrd_rhs.charge)								,
	molecule( lrd_rhs.molecule )						,
	EAS(lrd_rhs.EAS)									, 
	NAS(lrd_rhs.NAS)									,
	RAS(lrd_rhs.RAS)									,
	dual(lrd_rhs.dual)									,
	Softness_Dual(lrd_rhs.Softness_Dual)				,
	Hyper_softness(lrd_rhs.Hyper_softness)				,
	hardness_A(lrd_rhs.hardness_A)						,
	hardness_B(lrd_rhs.hardness_B)						,
	hardness_C(lrd_rhs.hardness_C)						,
	hardness_D(lrd_rhs.hardness_D)						,
	multifilic(lrd_rhs.multifilic)						,
	electrophilicity(lrd_rhs.electrophilicity)			,
	electron_density(lrd_rhs.electron_density)			,
	fukushima(lrd_rhs.fukushima)						{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(const local_rd_cnd& lrd_rhs){
	if( this!=&lrd_rhs){
		name			= lrd_rhs.name;
		ed 				= lrd_rhs.ed;
		mep				= lrd_rhs.mep;
		finite_diff		= lrd_rhs.finite_diff;
		charge			= lrd_rhs.charge;
		molecule		= lrd_rhs.molecule;
		EAS				= lrd_rhs.EAS;
		NAS				= lrd_rhs.NAS;
		RAS				= lrd_rhs.RAS;
		dual			= lrd_rhs.dual;
		Softness_Dual	= lrd_rhs.Softness_Dual;
		electron_density= lrd_rhs.electron_density;
		Hyper_softness	= lrd_rhs.Hyper_softness;
		hardness_A		= lrd_rhs.hardness_A;
		hardness_B		= lrd_rhs.hardness_B;
		hardness_C		= lrd_rhs.hardness_C;
		hardness_D		= lrd_rhs.hardness_D;
		multifilic		= lrd_rhs.multifilic;
		electrophilicity= lrd_rhs.electrophilicity;
		fukushima		= lrd_rhs.fukushima;
	}       	
	return *this;
}                
/***********************************************************************************/
local_rd_cnd::local_rd_cnd(local_rd_cnd&& lrd_rhs) noexcept	:
	name( move(lrd_rhs.name) )								,
	ed( move(lrd_rhs.ed) )									,
	mep( move(lrd_rhs.mep) )								,
	finite_diff( move(lrd_rhs.finite_diff) )				,
	charge( move(lrd_rhs.charge) )							,
	molecule( move(lrd_rhs.molecule) )						,
	EAS( move(lrd_rhs.EAS) )								, 
	NAS( move(lrd_rhs.NAS) )								,
	RAS( move(lrd_rhs.RAS) )								,
	dual( move(lrd_rhs.dual) )								,
	Softness_Dual( move(lrd_rhs.Softness_Dual) )			,
	electron_density( move(lrd_rhs.electron_density) )		,
	Hyper_softness( move(lrd_rhs.Hyper_softness) )			,
	hardness_A( move(lrd_rhs.hardness_A) )					,
	hardness_B( move(lrd_rhs.hardness_B) )					,
	hardness_C( move(lrd_rhs.hardness_C) )					,
	hardness_D( move(lrd_rhs.hardness_D) )					,
	multifilic( move(lrd_rhs.multifilic) )					,
	electrophilicity( move(lrd_rhs.electrophilicity) )		,
	fukushima( move(lrd_rhs.fukushima) )					{
}
/***********************************************************************************/
local_rd_cnd& local_rd_cnd::operator=(local_rd_cnd&& lrd_rhs) noexcept {
	if( this!=&lrd_rhs){
		name				= move(lrd_rhs.name);  
		ed 					= move(lrd_rhs.ed);
		mep 				= move(lrd_rhs.mep);
		finite_diff			= move(lrd_rhs.finite_diff);
		charge				= move(lrd_rhs.charge);
		molecule			= move(lrd_rhs.molecule);
		EAS					= move(lrd_rhs.EAS);
		NAS					= move(lrd_rhs.NAS);
		RAS					= move(lrd_rhs.RAS);
		dual				= move(lrd_rhs.dual);
		Softness_Dual		= move(lrd_rhs.Softness_Dual);
		Hyper_softness		= move(lrd_rhs.Hyper_softness);
		electron_density	= move(lrd_rhs.electron_density);
		hardness_A			= move(lrd_rhs.hardness_A);
		hardness_B			= move(lrd_rhs.hardness_B);
		hardness_C			= move(lrd_rhs.hardness_C);
		hardness_D			= move(lrd_rhs.hardness_D);
		multifilic			= move(lrd_rhs.multifilic);
		electrophilicity	= move(lrd_rhs.electrophilicity);
		fukushima			= move(lrd_rhs.fukushima);
	}       	
	return *this;
}
/***********************************************************************************/
local_rd_cnd operator-(const local_rd_cnd& lrd_lhs,const local_rd_cnd& lrd_rhs){
	local_rd_cnd Result = lrd_lhs;
	if ( lrd_lhs.molecule.num_of_atoms == lrd_rhs.molecule.num_of_atoms) {
		for (unsigned int i=0;i<lrd_lhs.molecule.num_of_atoms;i++){
			Result.EAS[i]				= lrd_lhs.EAS[i]				- lrd_rhs.EAS[i];
			Result.NAS[i]				= lrd_lhs.NAS[i]				- lrd_rhs.NAS[i];
			Result.RAS[i]				= lrd_lhs.RAS[i]				- lrd_rhs.RAS[i];
			Result.dual[i]				= lrd_lhs.dual[i]				- lrd_rhs.dual[i];
			Result.Softness_Dual[i]		= lrd_lhs.Softness_Dual[i] 		- lrd_rhs.Softness_Dual[i];
			Result.Hyper_softness[i]	= lrd_lhs.Hyper_softness[i] 	- lrd_rhs.Hyper_softness[i];
			Result.hardness_A[i]		= lrd_lhs.hardness_A[i]			- lrd_rhs.hardness_A[i];
			Result.hardness_B[i]		= lrd_lhs.hardness_B[i]			- lrd_rhs.hardness_B[i];
			Result.hardness_C[i]		= lrd_lhs.hardness_C[i]			- lrd_rhs.hardness_C[i];
			Result.hardness_D[i]		= lrd_lhs.hardness_D[i]			- lrd_rhs.hardness_D[i];
			Result.electron_density[i] 	= lrd_lhs.electron_density[i] 	- lrd_rhs.electron_density[i];
			Result.multifilic[i]		= lrd_lhs.multifilic[i]			- lrd_rhs.multifilic[i];				
			Result.electrophilicity[i]	= lrd_lhs.electrophilicity[i]	- lrd_rhs.electrophilicity[i];
			Result.fukushima[i]			= lrd_lhs.fukushima[i]			- lrd_rhs.fukushima[i];
		}
	}
	return Result;
}		
/***********************************************************************************/
void local_rd_cnd::calculate_Fukui(){
	int i =0;
	if ( !finite_diff ) {
		unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
		if ( qm_calc->data_ok ){
			#pragma omp parallel for default(shared) private(i) 
			for (i=0;i<molecule.num_of_atoms;i++){
				EAS[i]			= qm_calc->EAS_in_atoms(i,0);
				NAS[i]			= qm_calc->NAS_in_atoms(i,0);
				fukushima[i]	= qm_calc->fukushima(i,0);
				RAS[i]			= (EAS[i] + NAS[i])/2;
				dual[i]			= NAS[i] - EAS[i];
			}
		}else{
			m_log->write_warning("Skipping Local RD calculations!");
		}
		molecule = move( qm_calc->molecule );
	}
}
/***********************************************************************************/
void local_rd_cnd::calculate_Fukui_band(int band){
	unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
	unsigned int nof,i;
	nof  = molecule.num_of_atoms;
	if ( qm_calc->data_ok ){
		omp_set_num_threads(NP);
		#pragma omp parallel for default(shared) private(i)	
		for (i=0;i<nof;i++){
			EAS[i]		= qm_calc->EAS_in_atoms(i,band);
			NAS[i]		= qm_calc->NAS_in_atoms(i,band);
			fukushima[i]= qm_calc->fukushima(i,band);
			RAS[i]		= (EAS[i] + NAS[i])/2;
			dual[i]		= NAS[i] - EAS[i];
		}
	}else{
		m_log->write_warning("Skipping Local RD calculations!");
	}
	molecule = move( qm_calc->molecule );
	
}
/***********************************************************************************/
void local_rd_cnd::calculate_Fukui_EW(int band){
	unsigned int nof ,i;
	nof = molecule.num_of_atoms;
	unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
	if ( qm_calc->data_ok ){
		omp_set_num_threads(NP);
		#pragma omp parallel for default(shared) private(i) 
		for (i=0;i<nof;i++){
			EAS[i]			= qm_calc->EAS_in_atoms_EW(i);
			NAS[i] 			= qm_calc->NAS_in_atoms_EW(i);
			fukushima[i]	= qm_calc->fukushima(i,band);
			RAS[i] 			= (EAS[i] + NAS[i])/2;
			dual[i]			= NAS[i] - EAS[i];
		}
	}else{
		m_log->write_warning("Skipping Local RD calculations!");
	}
	molecule = move(qm_calc->molecule);
}
/***********************************************************************************/
void local_rd_cnd::calculate_RD(const global_rd& grd){
	for(unsigned int i=0;i<molecule.num_of_atoms;i++){
			Softness_Dual[i]	= dual[i]*grd.softness;
			Hyper_softness[i]	= RAS[i]*grd.softness;
			multifilic[i]		= dual[i]*grd.Electrophilicity;
			electrophilicity[i]	= NAS[i]*grd.Electrophilicity;
		}
}
/***********************************************************************************/
void local_rd_cnd::calculate_Hardness(const global_rd& grd){
	unsigned int nof = molecule.num_of_atoms;
	unsigned int i,j;
	
	omp_set_num_threads(NP);
	//-----------------------------------------------------
	//Estimating electron density
	if ( !finite_diff ){
		unique_ptr<QMdriver> qm_calc ( new QMdriver( move(molecule) ) );
		vector<Iatom> atoms = qm_calc->molecule.atoms;
		#pragma omp parallel for default(shared) private(i) 
		for(i=0;i<nof;i++){
			electron_density[i] = qm_calc->density_in_atoms(i);
		}
		molecule = move(qm_calc->molecule);
	}else{
		#pragma omp parallel for default(shared) private(i) 
		for(i=0;i<nof;i++){
			electron_density[i] = molecule.atoms[i].atomicN - molecule.atoms[i].charge;
		}
	}
	//----------------------------------------------------
	//calculating local hardness with method A (LCP)
	#pragma omp parallel for default(shared) private(i) 
	for(i=0;i<nof;i++){
		hardness_A[i] = EAS[i]*(grd.chemical_pot/nof)
						+ (electron_density[i]/nof)*grd.hardness*2;
	}
	double r  = 0;
	double xi, yi,zi;
		
	//-----------------------------------------------------
	//calculating local hardness with method B (electron-electron interaction)
	#pragma omp parallel for default(shared) private(i)
	for (i=0;i<nof;i++){
		for (j=0;j<nof;j++){
			if ( i != j ){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				hardness_B[i] += electron_density[j]/r;
			}
		}
		hardness_B[i] /= molecule.num_of_electrons;
	}
	//-----------------------------------------------------
	//calculating local hardness with method C (Fukui potential)
	r  = 0;
	for (i=0;i<nof;i++){
		for (j=0;j<nof;j++){
			if ( i != j){
				xi = molecule.atoms[i].xcoord - molecule.atoms[j].xcoord;
				xi *= xi; 
				yi = molecule.atoms[i].ycoord - molecule.atoms[j].ycoord;
				yi *= yi;
				zi = molecule.atoms[i].zcoord - molecule.atoms[j].zcoord;
				zi *= zi;
				r  = sqrt(xi+yi+zi);
				hardness_C[i] += EAS[j]/r;
			}
		}
	}
	//--------------------------------------------------------
	//calculating local hardness with method D (Fukui distribution)
	for(i=0;i<nof;i++){	hardness_D[i] = NAS[i]*grd.lumo_en - EAS[i]*grd.homo_en; }
	
}
/*************************************************************************************/
protein_lrd local_rd_cnd::rd_protein(const Iprotein& prot){
	
	unsigned int i,j;
	
	vector<residue_lrd> res_rd( prot.residues.size() );
	unsigned int nof = molecule.num_of_atoms;
	int cnt = 0;
	for( i=0; i<prot.residues.size(); i++ ){
		for( j=0; j<prot.residues[i].atom_s; j++ ){
			res_rd[i].rd_sum[0]	+= EAS[cnt];
			res_rd[i].rd_sum[1] += NAS[cnt];
			res_rd[i].rd_sum[2] += RAS[cnt];
			res_rd[i].rd_sum[3] += dual[cnt];
			res_rd[i].rd_sum[4] += Hyper_softness[cnt];
			res_rd[i].rd_sum[5] += hardness_A[cnt];
			res_rd[i].rd_sum[6] += hardness_B[cnt];
			res_rd[i].rd_sum[7] += hardness_C[cnt];
			res_rd[i].rd_sum[8] += hardness_D[cnt];
			res_rd[i].rd_sum[9] += multifilic[cnt];
			res_rd[i].rd_sum[10] += electrophilicity[cnt];
			res_rd[i].rd_sum[11] += fukushima[cnt];
			res_rd[i].rd_sum[12] += electron_density[cnt];
			res_rd[i].rd_sum[13] += Softness_Dual[cnt];
			res_rd[i].rd_sum[14] += molecule.atoms[cnt].charge;
			cnt++;
		}
	}

	for( i=0; i<res_rd.size(); i++ ){
		for( j=0; j<res_rd[i].rd_sum.size(); j++ ){
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
	prot.load_b_column(EAS);
	prot.title	= rd_results.name + "_EAS";
	prot.remark	= "Electrophilic Attack Suscptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(NAS);
	prot.title	= rd_results.name +"_NAS";
	prot.remark	= "Nucleophilic Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(RAS);
	prot.title	= rd_results.name +"_RAS";
	prot.remark	= "Radical Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(dual);
	prot.title	= rd_results.name +"_Netphilicity";
	prot.remark	= "Local netphilicity Descriptor Attack Susceptibility at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_A);
	prot.title	= rd_results.name + "_hardness_A";
	prot.remark	= "Local Hardness (local chemical potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_B);
	prot.title 	= rd_results.name +"_hardness_B";
	prot.remark	= "Local Hardness (Electron-electron potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_C);
	prot.title	= rd_results.name +"_hardness_C";
	prot.remark	= "Local Hardness (Fukui potential based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(hardness_D);
	prot.title	= rd_results.name +"_hardness_D";
	prot.remark	= "Local Hardness (Fukui distribution based) at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(Softness_Dual);
	prot.title	= rd_results.name +"_softness";
	prot.remark	= "Local Softness Dual at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(Hyper_softness);
	prot.title	= rd_results.name +"_hypersoftness";
	prot.remark	= "Local Hyper Softness at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(fukushima);
	prot.title	= rd_results.name +"_fukushima";
	prot.remark	= "Localization of frontier band orbitals at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(multifilic);
	prot.title	= rd_results.name +"_multiphilic";
	prot.remark	= "Local Multiphilic descriptor  at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	prot.load_b_column(electrophilicity);
	prot.title	= rd_results.name +"_electrophilicity";
	prot.remark	= "Local Electrophilicity descriptor  at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);	
	prot.load_b_column(electron_density);
	prot.title	= rd_results.name +"_electron_density";
	prot.remark	= "Electronic density per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors in PDBs.\n");
	m_log->inp_delim(1);

	vector<double> pcharges( molecule.atoms.size() );
	for(unsigned int i=0; i<molecule.atoms.size();i++) {
		pcharges[i] = molecule.atoms[i].charge;
	}
	prot.load_b_column(pcharges);
	prot.title 	= rd_results.name +"_mep";
	prot.remark	= "Partial charge per atom at b-factor column calculated by PRIMoRDiA";
	rd_results.models.emplace_back(prot);
	//rd_results.write_pdb(rd_results.name +"_RD.pdb");
	fs::create_directory(name+"_PDB_RD");
	rd_results.write_models(name+"_PDB_RD");
}
/*************************************************************************************/
void local_rd_cnd::write_LRD(){
	std::string temps;
	if ( finite_diff ) { temps = name+"FD.lrd";
	}else{	temps = name+"FOA.lrd";	}
	const char* names = temps.c_str();
	std::ofstream lrd_file;
	lrd_file.open(names);
	
	lrd_file.precision(8);
	lrd_file << std::fixed;
	lrd_file << name << " " << "\n" <<  std::left;
	lrd_file << "n atom charge Nucleophilicity Electrophilicity RAS Netphilicity Softness Hardness_A Hardness_B Hardness_C Hardness_D Multiphilic T_Electrohilicity Fukushima Electron_density Softness_dual\n";
	
	for(int i=0;i<molecule.num_of_atoms;i++){
		lrd_file	<< (i+1) 
					<< " "
					<< molecule.atoms[i].element	<< " "
					<< molecule.atoms[i].charge		<< " "
					<< EAS[i]						<< " "
					<< NAS[i]						<< " "
					<< RAS[i]						<< " "
					<< dual[i]						<< " "
					<< Hyper_softness[i]			<< " "
					<< hardness_A[i]				<< " "
					<< hardness_B[i]				<< " "
					<< hardness_C[i]				<< " "
					<< hardness_D[i]				<< " "
					<< multifilic[i]				<< " "
					<< electrophilicity[i]			<< " "
					<< fukushima[i] 				<< " "
					<< electron_density[i]			<< " "
					<< Softness_Dual[i]				<< "\n";
	}
	lrd_file.close();
	m_log->input_message("Finishing the writting of the condensed local reactivity descriptors.\n");
	m_log->inp_delim(2);
}
//================================================================================
//END OF FILE
//================================================================================