//QMdriver.cpp
//------------------------------------------

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
//C++ headers
#include <iostream>
#include <cmath>
//------------------------------------------
//PRIMoRDiA headers
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h" 
#include "../include/QMdriver.h"
//------------------------------------------
using std::move;
using std::cout;
using std::endl;
using std::vector;
using std::abs;
using std::string;
/***************************************************************/
QMdriver::QMdriver(){}
/***************************************************************/
QMdriver::QMdriver(Imolecule&& mol) noexcept:
	molecule( move(mol) )					,
	data_ok(true)							{
	
	m_log->input_message("Checking some molecular information to performa the local condensed to atoms reactivity descriptors");
	
	if ( molecule.m_overlap.size() <= 0		||
		 molecule.coeff_MO.size() <= 0		|| 
		 molecule.orb_energies.size() <=0  ) {
		cout << "There are required molecular information that were not properly stored!" << endl;
		m_log->input_message("Molecular information checking for local consended reactivity descriptors ended with problems.");
		m_log->input_message("Current molecular storage containers size: ");
		m_log->input_message("Overlap matrix: ");
		m_log->input_message(int ( molecule.m_overlap.size() ) );
		m_log->input_message("Orbital energies: ");
		m_log->input_message(int ( molecule.orb_energies.size() ) );
		m_log->input_message("Molecular orbitals coefficient: ");
		m_log->input_message(int ( molecule.coeff_MO.size() ) );
		data_ok = false;
	}


}
/***************************************************************/
double QMdriver::MO_in_atoms(int atom,int MO,bool beta){
	double value	= 0.0;
	int init_orb	= 0;
	int n_aorbs 	= 0;
	unsigned int j;
	int ao = molecule.num_of_ao;
	if ( ao == 0 ) { ao = molecule.MOnmb; }
	
	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	if ( !beta ){
		for(int mu=init_orb;mu<n_aorbs;mu++){
				value += molecule.coeff_MO[ao*MO + mu]*
						 molecule.coeff_MO[ao*MO + mu];
		}
	}else if ( beta ){
		for(int mu=init_orb;mu<n_aorbs;mu++){
				value += molecule.coeff_MO_beta[ao*MO + mu]*
						 molecule.coeff_MO_beta[ao*MO + mu];
		}
	}
	return value;
}
/***************************************************************/
double QMdriver::EAS_in_atoms(int atom,int band){
	double value	= 0.0;
	int init_orb	= 0;
	int n_aorbs	= 0;
	unsigned int j;
	int ao = molecule.num_of_ao;
	int homon = abs(molecule.homoN);
	int cnt = 0;
	
	if ( ao == 0 ) { ao = molecule.MOnmb; }
	
	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	
	
	int init = homon-band;
	for( int i=init;i<=homon;i++){
		if ( molecule.orb_energies[i] >= molecule.homo_energy-energy_crit ){
			for(int mu=init_orb;mu<n_aorbs;mu++){
				for (int nu=init_orb;nu<n_aorbs;nu++) {
					value +=molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + nu]*
							molecule.m_overlap[nu+(mu*(mu+1) )/2];
				}
			}
			cnt++;
		}
	}
	if ( atom == 0 && band > 0 ){
		m_log->input_message("Number of occupied MO used to calculate condensed to atom descriptors: ");
		m_log->input_message(cnt);
	}
	if ( molecule.betad ){
		cnt = 0;
		for( int i=init;i<=homon;i++){
			if ( molecule.orb_energies_beta[i] >= molecule.homo_energy-energy_crit ){
				for(int mu=init_orb;mu<n_aorbs;mu++){
					for (int nu=init_orb;nu<n_aorbs;nu++) {
						value +=molecule.coeff_MO_beta[ao*i + mu]*
								molecule.coeff_MO_beta[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1) )/2];
					}
				}
				cnt++;
			}
		}
		if ( atom == 0  && band > 0 ){
			m_log->input_message("Number of occcupied virtual MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(cnt);
		}
		value /= 2;
	}
	return value;
}
/***************************************************************/
double QMdriver::NAS_in_atoms(int atom,int band){
	double value	= 0.0;
	int init_orb	= 0;
	int n_aorbs		= 0;
	unsigned int j;
	int ao			= molecule.num_of_ao;
	int lumon		= abs(molecule.lumoN);
	int fin			= lumon+band;
	int cnt			= 0;
	
	if ( ao == 0 ) { ao = molecule.MOnmb; }
	
	if ( atom == 0 ) {
		init_orb		= 0;
		n_aorbs	= molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	
	for( int i=lumon;i<=fin;i++){
		if ( molecule.orb_energies[i] <= molecule.lumo_energy+energy_crit ){
			for(int mu=init_orb;mu<n_aorbs;mu++){
				for (int nu=init_orb;nu<n_aorbs;nu++) {
					value +=molecule.coeff_MO[ao*i + mu]*
								molecule.coeff_MO[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1))/2];
				}
			}
			cnt++;
		}
	}
	if ( atom == 0  && band > 0 ){
			m_log->input_message("Number of  virtual MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(cnt);
	}
	if ( molecule.betad ) {
		cnt = 0;
		for( int i=lumon;i<=fin;i++){
			if ( molecule.orb_energies_beta[i] <= molecule.lumo_energy+energy_crit ){
				for(int mu=init_orb;mu<n_aorbs;mu++){
					for (int nu=init_orb;nu<n_aorbs;nu++) {
						value +=molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
				cnt++;
			}
		}
		if ( atom == 0  && band > 0 ){
			m_log->input_message("Number of  virtual beta MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(cnt);
		}
		value /=2;
	}
	return value;
}
/***************************************************************/
double QMdriver::EAS_in_atoms_EW(int atom){
	double value		= 0.0;	
	int init_orb		= 0;
	int n_aorbs		= 0;
	unsigned int j;
	int ao 				= molecule.num_of_ao;
	int homon 		= molecule.homoN;
	int mo_count 	= 0;
	if ( ao == 0 ) { ao = molecule.MOnmb; }

	double pre_coef = exp(-abs( energy_crit ) );
	
	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	
	double coefficient = 0.0;
	
	for( int i=0;i<=homon;i++){
		coefficient = exp(-abs(molecule.orb_energies[i]-molecule.homo_energy ) );
		if ( coefficient > pre_coef ){ 
			mo_count++;
			for(int mu=init_orb;mu<n_aorbs;mu++){
				for (int nu=init_orb;nu<n_aorbs;nu++) {
					value +=coefficient*
								molecule.coeff_MO[ao*i + mu]*
								molecule.coeff_MO[ao*i + nu]*
								molecule.m_overlap[nu+(mu*(mu+1))/2];	
				}
			}
		}
	}
	if ( atom == 0  ){
		m_log->input_message("Number of  occupied  MO used to calculate condensed to atom descriptors: ");
		m_log->input_message(mo_count);
	}
	if ( molecule.betad ){
		mo_count = 0;
		for( int i=0;i<=homon;i++){
			coefficient = exp(-abs(molecule.orb_energies_beta[i]-molecule.homo_energy ) );
			if ( coefficient > pre_coef ){ 
				mo_count++;
				for(int mu=init_orb;mu<n_aorbs;mu++){
					for (int nu=init_orb;nu<n_aorbs;nu++) {
						value +=coefficient*
									molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
			}
		}
		if ( atom == 0 ){
			m_log->input_message("Number of  occupied  beta MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(mo_count);
		}
		value /=2;
	}
	return value;
}
/***************************************************************/
double QMdriver::NAS_in_atoms_EW(int atom){
	double value		= 0.0;
	int init_orb		= 0;
	int n_aorbs		= 0;
	unsigned int j;
	int ao				= molecule.num_of_ao;
	int lumon			= abs(molecule.homoN);
	int mo_count	= 0;
	double pre_coef = exp(-abs( energy_crit ) );
	if ( ao == 0 ) { ao = molecule.MOnmb; }

	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	
	int fin = molecule.orb_energies.size();
	
	double coefficient = 0.0;
	for( int i=lumon;i<fin;i++){
		coefficient = exp(-abs(molecule.orb_energies[i]-molecule.lumo_energy ) );
		if ( coefficient > pre_coef ){ 
			mo_count++;
			for(int mu=init_orb;mu<n_aorbs;mu++){
				for (int nu=init_orb;nu<n_aorbs;nu++) {
					value +=coefficient*
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
	}
	if ( molecule.betad ){
		mo_count = 0;
		for( int i=lumon;i<=fin;i++){
			coefficient = exp(-abs(molecule.orb_energies_beta[i]-molecule.lumo_energy ) );
			if ( coefficient >  pre_coef  ){ 
				mo_count++;
				for(int mu=init_orb;mu<n_aorbs;mu++){
					for (int nu=init_orb;nu<n_aorbs;nu++) {
						value +=coefficient*
									molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu]*
									molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
			}
		}
		if ( atom == 0 ){
			m_log->input_message("Number of  virtual  beta MO used to calculate condensed to atom descriptors: ");
			m_log->input_message(mo_count);
		}
		value /= 2;
	}
	return value;
}
/***************************************************************/
double QMdriver::density_in_atoms(int atom){
	unsigned i,j	= 0;
	double value	= 0.0;
	int size_alf		= 0;
	int size_bet	= -1;
	int init_orb	= 0;
	int n_aorbs 	= 0;
	int ao			= molecule.num_of_ao;
	if ( ao == 0 ) { ao = molecule.MOnmb; }

	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	double occ = 2.0;
	if ( molecule.betad ) { occ = 1.0 ;}
	for(i=0;i<=molecule.homoN;i++){
		for(int mu=init_orb;mu<n_aorbs;mu++){
			for (int nu=init_orb;nu<n_aorbs;nu++) {
				value +=occ*
							molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + nu];
							//molecule.m_overlap[ nu+( mu*(mu+1) ) / 2 ];
			}
		}
	}
	if ( molecule.occupied_beta.size() > 0 ){
		for(i=0;i<molecule.homoN;i++){
			if ( molecule.occupied_beta[i] > 0 ){
				for(int mu=init_orb;mu<n_aorbs;mu++){
					for (int nu=init_orb;nu<n_aorbs;nu++) {
						value +=molecule.coeff_MO_beta[ao*i + mu]*
									molecule.coeff_MO_beta[ao*i + nu];
									///molecule.m_overlap[nu+(mu*(mu+1))/2];
					}
				}
			}
		}
	}
	return value;
}
/***************************************************************/
double QMdriver::fukushima(int atom, int band){
	double value				= 0.0;
	unsigned int init_orb	= 0;
	unsigned int n_aorbs	= 0;
	unsigned int j;
	unsigned int ao			= molecule.num_of_ao;
	
	if ( atom == 0 ) {
		init_orb = 0;
		n_aorbs  = molecule.atoms[0].orbitals.size();
	}else{
		for(j=0;j<atom;j++){
			init_orb += molecule.atoms[j].norb;
		}
		for(j=0;j<=atom;j++){
			n_aorbs  += molecule.atoms[j].norb;
		}
	}
	for( int i=(molecule.homoN-band);i<=molecule.homoN;i++){
		if ( molecule.orb_energies[i] >= molecule.homo_energy-energy_crit ){
			for(int mu=init_orb;mu<n_aorbs;mu++){
				value +=molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + mu];
			}
		}
	}
	if ( molecule.betad ){
		for( int i=(molecule.homoN-band);i<=molecule.homoN;i++){
			if ( molecule.orb_energies_beta[i] >= molecule.homo_energy-energy_crit ){
				for(int mu=init_orb;mu<n_aorbs;mu++){
					value +=molecule.coeff_MO_beta[ao*i + mu]*
								molecule.coeff_MO_beta[ao*i + mu];
				}
			}
		}
	}
	for( int i=molecule.lumoN;i<=molecule.lumoN+band;i++){
		if ( molecule.orb_energies[i] <= molecule.lumo_energy+energy_crit ){
			for(int mu=init_orb;mu<n_aorbs;mu++){
				value +=molecule.coeff_MO[ao*i + mu]*
							molecule.coeff_MO[ao*i + mu];
			}
		}
	}
	if ( molecule.betad ){
		for( int i=molecule.lumoN;i<=molecule.lumoN+band;i++){
			if ( molecule.orb_energies_beta[i] <= molecule.lumo_energy+energy_crit ){
				for(int mu=init_orb;mu<n_aorbs;mu++){
					value +=molecule.coeff_MO_beta[ao*i + mu]*
								molecule.coeff_MO_beta[ao*i + mu];
				}
			}
		}
	}
	return value;
}
/***************************************************************/
QMdriver::~QMdriver(){}
//===============================================================
//END OF FILE 
//===============================================================

