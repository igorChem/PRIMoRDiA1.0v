// source file for the Imolecule class to represent and deal with moleculae information
// imolecule.cpp 

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
//c++ headers
#include <iostream> 
#include <algorithm>
#include <vector>
#include <string>  
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
//--------------------------------------------------------
//PRIMoRDiA headers
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
//--------------------------------------------------------
using std::cout; 
using std::endl;
using std::move;
using std::string;
using std::vector;
//=======================================================================================
/***************************************************************************************/
Imolecule::Imolecule()		:
	name("nonamed")			,
	num_of_atoms(0)			,
	num_of_ao(0)			,
	f_chg(0)				,
	num_of_electrons(0)		,
	molar_mass(0.0)			, 
	mol_charge(0.0)			,
	energy_tot(0.0)			,
	homo_energy(0.0)		,
	lumo_energy(1000.0)		,
	total_dipmoment(0.0)	,
	heat_of_formation(0.0)	,
	MOnmb(0)				,
	MOnmb_beta(0)			,
	homoN(0.0)				,
	lumoN(0.0)				,
	normalized(false)		,
	bohr(false)				,
	betad(false)			{
	
	for (unsigned int i=0; i<3;i++){
		ver_inf[i]				= 0.0;
		ver_sup[i]				= 0.0;
		dipole_moment[i]	= 0.0;
	}
} 
/***************************************************************************************/
Imolecule::Imolecule(const Imolecule& rhs_molecule)		:
	name(rhs_molecule.name)								,
	num_of_atoms(rhs_molecule.num_of_atoms)				,
	num_of_ao(rhs_molecule.num_of_ao)					,
	f_chg(rhs_molecule.f_chg)							,
	mol_charge(rhs_molecule.mol_charge)					,
	num_of_electrons(rhs_molecule.num_of_electrons)		,
	molar_mass(rhs_molecule.molar_mass)					,
	energy_tot(rhs_molecule.energy_tot)					,
	homo_energy(rhs_molecule.homo_energy)				,
	lumo_energy(rhs_molecule.lumo_energy)				, 
	heat_of_formation(rhs_molecule.heat_of_formation)	,
	MOnmb(rhs_molecule.MOnmb)							,
	MOnmb_beta(rhs_molecule.MOnmb_beta)					,
	homoN(rhs_molecule.homoN)							,
	lumoN(rhs_molecule.lumoN)							,
	normalized(rhs_molecule.normalized)					,
	bohr(rhs_molecule.bohr)								,
	betad(rhs_molecule.betad)							,
	orb_energies(rhs_molecule.orb_energies)				,
	orb_energies_beta(rhs_molecule.orb_energies_beta)	,
	coeff_MO(rhs_molecule.coeff_MO)						,
	coeff_MO_beta(rhs_molecule.coeff_MO_beta)			,
	m_dens(rhs_molecule.m_dens)							,
	beta_dens(rhs_molecule.beta_dens)					,
	m_overlap(rhs_molecule.m_overlap)					,
	occupied(rhs_molecule.occupied)						,
	occupied_beta(rhs_molecule.occupied_beta)			,
	atoms(rhs_molecule.atoms)							{
	
	for(int i = 0 ; i < 3; i++) {
		dipole_moment[i]	= rhs_molecule.dipole_moment[i];
		ver_inf[i]			= rhs_molecule.ver_inf[i];
		ver_sup[i]			= rhs_molecule.ver_sup[i];
	}
}
/***************************************************************************************/
Imolecule& Imolecule::operator=(const Imolecule& rhs_molecule){
	if(this != &rhs_molecule) {	
		name				= rhs_molecule.name;
		num_of_atoms		= rhs_molecule.num_of_atoms;
		num_of_electrons	= rhs_molecule.num_of_electrons;
		num_of_ao			= rhs_molecule.num_of_ao;
		f_chg				= rhs_molecule.f_chg;
		molar_mass			= rhs_molecule.molar_mass;
		energy_tot			= rhs_molecule.energy_tot;
		homo_energy			= rhs_molecule.homo_energy;
		lumo_energy			= rhs_molecule.lumo_energy;
		heat_of_formation	= rhs_molecule.heat_of_formation;
		MOnmb				= rhs_molecule.MOnmb;
		MOnmb_beta			= rhs_molecule.MOnmb_beta;
		homoN				= rhs_molecule.homoN;
		lumoN				= rhs_molecule.lumoN;
		normalized			= rhs_molecule.normalized;
		bohr				= rhs_molecule.bohr;
		betad				= rhs_molecule.betad;
		orb_energies		= rhs_molecule.orb_energies;
		orb_energies_beta	= rhs_molecule.orb_energies_beta;
		mol_charge			= rhs_molecule.mol_charge;
		coeff_MO			= rhs_molecule.coeff_MO;
		coeff_MO_beta		= rhs_molecule.coeff_MO_beta;
		m_dens				= rhs_molecule.m_dens;
		beta_dens			= rhs_molecule.beta_dens;
		m_overlap			= rhs_molecule.m_overlap;
		occupied			= rhs_molecule.occupied;
		occupied_beta		= rhs_molecule.occupied_beta;
		atoms				= rhs_molecule.atoms;
		
		for(int i = 0 ; i < 3; i++) {
			dipole_moment[i] 	= rhs_molecule.dipole_moment[i];
			ver_inf[i]			= rhs_molecule.ver_inf[i];
			ver_sup[i]			= rhs_molecule.ver_sup[i];
		}
	}
	return *this;
}
/***************************************************************************************/
Imolecule::Imolecule(Imolecule&& rhs_molecule) noexcept			:
	name(move (rhs_molecule.name) )								,
	num_of_atoms(rhs_molecule.num_of_atoms)						,
	num_of_electrons(rhs_molecule.num_of_electrons)				,
	num_of_ao(rhs_molecule.num_of_ao)							,
	f_chg(rhs_molecule.f_chg)									,
	molar_mass(rhs_molecule.molar_mass)							,
	mol_charge(rhs_molecule.mol_charge)							,
	energy_tot(rhs_molecule.energy_tot)							,
	homo_energy(rhs_molecule.homo_energy)						,
	lumo_energy(rhs_molecule.lumo_energy)						, 
	heat_of_formation(rhs_molecule.heat_of_formation)			,
	MOnmb(rhs_molecule.MOnmb)									,
	MOnmb_beta(rhs_molecule.MOnmb_beta)							,
	homoN(rhs_molecule.homoN)									, 
	lumoN(rhs_molecule.lumoN)									,
	normalized(rhs_molecule.normalized)							,
	bohr(rhs_molecule.bohr)										,
	betad(rhs_molecule.betad)									,
	orb_energies( move(rhs_molecule.orb_energies) )				,
	orb_energies_beta( move(rhs_molecule.orb_energies_beta) )	,
	coeff_MO( move(rhs_molecule.coeff_MO) )						,
	coeff_MO_beta( move(rhs_molecule.coeff_MO_beta) )			,
	m_dens( move(rhs_molecule.m_dens) )							,
	beta_dens( move(rhs_molecule.beta_dens) )					,
	m_overlap( move(rhs_molecule.m_overlap) )					,
	occupied( move(rhs_molecule.occupied))						,
	occupied_beta( move(rhs_molecule.occupied_beta) )			,
	atoms( move(rhs_molecule.atoms) )							{
	
	for(int i = 0 ; i < 3; i++) {
		dipole_moment[i]= rhs_molecule.dipole_moment[i];
		ver_inf[i]		= rhs_molecule.ver_inf[i];
		ver_sup[i]		= rhs_molecule.ver_sup[i];
	}
}
/***************************************************************************************/
Imolecule& Imolecule::operator=(Imolecule&& rhs_molecule) noexcept {
	if(this != &rhs_molecule) {	
		name				= move(rhs_molecule.name);
		num_of_atoms		= rhs_molecule.num_of_atoms;
		num_of_electrons	= rhs_molecule.num_of_electrons;
		num_of_ao			= rhs_molecule.num_of_ao;
		f_chg				= rhs_molecule.f_chg; 
		mol_charge			= rhs_molecule.mol_charge;
		molar_mass			= rhs_molecule.molar_mass;
		energy_tot			= rhs_molecule.energy_tot;
		homo_energy			= rhs_molecule.homo_energy;
		lumo_energy			= rhs_molecule.lumo_energy;
		heat_of_formation	= rhs_molecule.heat_of_formation;
		MOnmb				= rhs_molecule.MOnmb;
		MOnmb_beta			= rhs_molecule.MOnmb_beta;
		homoN				= rhs_molecule.homoN;
		lumoN				= rhs_molecule.lumoN;
		normalized			= rhs_molecule.normalized;
		bohr				= rhs_molecule.bohr;
		betad				= rhs_molecule.betad;
		orb_energies		= move(rhs_molecule.orb_energies);
		orb_energies_beta	= move(rhs_molecule.orb_energies_beta);
		coeff_MO			= move(rhs_molecule.coeff_MO);
		coeff_MO_beta		= move(rhs_molecule.coeff_MO_beta);
		m_dens				= move(rhs_molecule.m_dens);
		beta_dens			= move(rhs_molecule.beta_dens);
		m_overlap			= move(rhs_molecule.m_overlap);
		occupied			= move(rhs_molecule.occupied);
		occupied_beta		= move(rhs_molecule.occupied_beta);
		atoms				= move(rhs_molecule.atoms);
		
		for(int i = 0 ; i < 3; i++) {
			dipole_moment[i]= rhs_molecule.dipole_moment[i];
			ver_inf[i]		= rhs_molecule.ver_inf[i];
			ver_sup[i]		= rhs_molecule.ver_sup[i];
		}
	}
	return *this;
}
/***************************************************************************************/
void Imolecule::add_atom(double x	,
						 double y	, 
						 double z	, 
						 string typ){
							 
	atoms.emplace_back( x,y,z,move(typ) );
	num_of_atoms = atoms.size();
}
/***************************************************************************************/
void Imolecule::add_atom(Iatom atom) {
	atoms.emplace_back( move(atom) ) ; 
	num_of_atoms = atoms.size();
}
/***************************************************************************************/
void Imolecule::print_coordinates() {
	for(int i=0;i<num_of_atoms;i++){
			 cout 	<< "x # " << i << ": " << atoms[i].xcoord << " "
					<< "y # " << i << ": " << atoms[i].ycoord << " "
					<< "z # " << i << ": " << atoms[i].zcoord << endl;
	}
}
/***************************************************************************************/
void Imolecule::write_xyz() {
	string xyz_name;
	xyz_name = name + ".xyz";
	std::ofstream xyz_file(xyz_name.c_str());
	xyz_file << num_of_atoms << "\n";
	xyz_file << "xyz file generated by Imolecule class by Barden/2020" << "\n";
	m_log->input_message("Writing xyz!");
	for(int i=0;i<num_of_atoms;i++){ 
		xyz_file.precision(7);
		xyz_file << std::fixed;
		xyz_file << std::setw(3) << std::left << atoms[i].element
		<< "  " 
	    << std::setw(11) << std::right << atoms[i].xcoord 
	    << "  "
	    << std::setw(11) << std::right << atoms[i].ycoord
	    << "  "
	    << std::setw(11) << std::right << atoms[i].zcoord
	    << "\n";
	}
}
/***************************************************************************************/
void Imolecule::mol_vert_up(){
	ver_sup[0] = atoms[0].xcoord;
	ver_sup[1] = atoms[0].ycoord;
	ver_sup[2] = atoms[0].zcoord;
	
	ver_inf[0] = atoms[0].xcoord;
	ver_inf[1] = atoms[0].ycoord;
	ver_inf[2] = atoms[0].zcoord;
	
	m_log->input_message("Definind the grid vertices.");
	
	for (unsigned int i=0;i<atoms.size();i++){
		if ( atoms[i].xcoord > ver_sup[0] ) ver_sup[0] = atoms[i].xcoord;
		if ( atoms[i].ycoord > ver_sup[1] ) ver_sup[1] = atoms[i].ycoord;
		if ( atoms[i].zcoord > ver_sup[2] ) ver_sup[2] = atoms[i].zcoord;
		if ( atoms[i].xcoord < ver_inf[0] ) ver_inf[0] = atoms[i].xcoord;
		if ( atoms[i].ycoord < ver_inf[1] )	ver_inf[1] = atoms[i].ycoord;
		if ( atoms[i].zcoord < ver_inf[2] )	ver_inf[2] = atoms[i].zcoord;
	}
}
/***************************************************************************************/
vector<double> Imolecule::extract_MO(int MO,bool beta){
	int nMO 			= 0;
	if ( !beta ) nMO 	= MOnmb; 
	else         nMO 	= MOnmb_beta;
	vector<double> res_mo(nMO);

	for (int i=0;i<nMO;i++){
		if ( !beta ) res_mo[i] = coeff_MO[MO*MOnmb+i];
		else         res_mo[i] = coeff_MO_beta[MO*MOnmb_beta+i];
	}
	return res_mo;
}
/***************************************************************************************/
void Imolecule::print(){
	std::cout << "Molecule's name: " << name << std::endl;
	this->print_coordinates();
	std::cout << "The coeff are normalized :"	<< normalized	<< std::endl;
	std::cout << "LUMO energy  : "				<< lumo_energy 	<< std::endl;
	std::cout << "HOMO energy  : "				<< homo_energy	<< std::endl;
	std::cout << "Total energy : "				<< energy_tot	<< std::endl;
	std::cout << "Printing the first ten values of the vector in this object" << std::endl;
	
	for (unsigned int o=0;o<num_of_atoms;o++) atoms[o].print();
	for (unsigned int i=0;i<orb_energies.size();i++) {
		std::cout 	<< orb_energies[i]	<< " " 
						<< coeff_MO[i]	<< " " 
						<< m_overlap[i] << " "
						<< occupied[i] 	<< std::endl;
	}
}
/***************************************************************************************/
void Imolecule::ang_to_bohr(){
	const double angtobohr = 1.0/0.52917726;
	for(unsigned int i=0;i<num_of_atoms;i++){
		atoms[i].xcoord = atoms[i].xcoord*angtobohr;
		atoms[i].ycoord = atoms[i].ycoord*angtobohr;
		atoms[i].zcoord = atoms[i].zcoord*angtobohr;
	}
	bohr = true;
	m_log->input_message("Cartesian molecular coordinates converted to bohr.");
}
/***************************************************************************************/
void Imolecule::bohr_to_ang(){
	const double bohrtoang = 0.52917726;
	for(unsigned int i=0;i<num_of_atoms;i++){
		atoms[i].xcoord = atoms[i].xcoord*bohrtoang;
		atoms[i].ycoord = atoms[i].ycoord*bohrtoang;
		atoms[i].zcoord = atoms[i].zcoord*bohrtoang;
	}
	bohr = false;
	m_log->input_message("Cartesian molecular coordinates converted to angstron.");
}
/***************************************************************************************/
void Imolecule::norm_orbs(){
	for (unsigned int i=0;i<atoms.size();i++){
		for (unsigned int j=0;j<atoms[i].orbitals.size();j++){
			atoms[i].orbitals[j].normalize();
		}
	}
}
/***************************************************************************************/
int Imolecule::get_ao_number(){
	num_of_ao = 0;
	for(unsigned int i=0;i<atoms.size();i++){
		for (unsigned int j=0;j<atoms[i].orbitals.size();j++) num_of_ao++;
	}
	return num_of_ao;
}
/***************************************************************************************/
double Imolecule::get_homo(){
	int homo_nalfa = homoN;
	int homo_nbeta = homoN;
	for (unsigned int i=0;i<occupied.size();i++){ 
		if ( occupied[i]>=1 )
			homo_nalfa = i;
	}
	for (unsigned int j=0;j<occupied_beta.size();j++) { 
		if ( occupied_beta[j]>=1 )
			homo_nbeta = j;
	}
	
	double homo_alfa  = orb_energies[homo_nalfa];
	double homo_beta  = -1000.0;
	
	if ( orb_energies_beta.size() > 0 ){
		if ( homo_nbeta != 0 ) homo_beta  = orb_energies_beta[homo_nbeta]; 
	}
	
	if ( homo_alfa >= homo_beta ) {
		homoN 			= homo_nalfa;
		homo_energy = orb_energies[homo_nalfa];
	}
	else if ( homo_beta >  homo_alfa ){
		homoN 			= homo_nbeta;
		homo_energy = orb_energies_beta[homo_nbeta];
	}
	return homo_energy;
}
/***************************************************************************************/
double Imolecule::get_lumo(){
	int lumo_nalfa		= lumoN;
	int lumo_nbeta 		= lumoN;
	double lumo_energy_b= 1000;
	unsigned int i,j 	= 0;
	vector<double> energy_unn;
	vector<double> energy_unn_beta;

	for (i=0;i<MOnmb;i++){ 
		if ( occupied[i] < 1) {
			if ( lumo_nalfa == 0 ) { lumo_nalfa = i; }
			energy_unn.push_back(orb_energies[i]);
		}
	}
	for ( j=0;j<MOnmb_beta;j++) {
		if ( occupied_beta[j] < 1 ) {
			if ( lumo_nbeta == 0 ) { lumo_nbeta = j; }
			energy_unn_beta.push_back(orb_energies_beta[j]);
		}
	}
	
	lumo_energy = min_dvec(energy_unn);
	
	if ( betad ) {lumo_energy_b	= min_dvec(energy_unn_beta); }
	
	if ( lumo_energy_b < lumo_energy ){
		lumo_energy = lumo_energy_b;
		lumoN 			= lumo_nbeta;
	}else{
		lumoN = lumo_nalfa;
	}
	return lumo_energy;
}
/***************************************************************************************/
void Imolecule::print_basis(){
	for(unsigned int i=0;i<atoms.size();i++){
		for(unsigned int j=0;j<atoms[i].orbitals.size();j++)	atoms[i].orbitals[j].print();
	}
}
/***************************************************************************************/
void Imolecule::update(){
	for(unsigned int i=0;i<atoms.size();i++) {
		mol_charge	+= atoms[i].charge;
		molar_mass	+= atoms[i].atomic_mass;
	}
	
	//--------------------------------------------
	if ( num_of_electrons == 0 ){
		m_log->input_message("None electrons found in for the molecule!");
		m_log->input_message("Starting to count based on atomic numbers!!");
		for(unsigned int i=0;i<atoms.size();i++){
			num_of_electrons += atoms[i].atomicN;
		}
	}
	//------------------------------------------
	f_chg			= std::round(mol_charge);
	num_of_atoms	= atoms.size();
	this->get_ao_number();
	unsigned int noe= num_of_electrons;
	int i;
	if ( !betad ){
		noe					+=-1*f_chg;
		occupied.resize(num_of_ao);
		for(i=0;i<noe/2;i++){ occupied[i] = 2; }
	}else{
		if ( noe % 2 == 0 ) {
			noe +=-1*f_chg;
		}
		occupied.resize(num_of_ao);
		occupied_beta.resize(num_of_ao);
		for(i=0;i<(noe+1)/2;i++) occupied[i] = 1;
		for(i=0;i<(noe-1)/2;i++) occupied_beta[i] = 1;
	}
	this->get_homo();
	this->get_lumo();
	this->norm_orbs();
	m_log->input_message("Formal Charge: ");
	m_log->input_message(int(f_chg));
}
/***************************************************************************************/
void Imolecule::clear(){
	vector<double>().swap(orb_energies);      
	vector<double>().swap(orb_energies_beta);
	vector<double>().swap(coeff_MO);
	vector<double>().swap(coeff_MO_beta);
	vector<double>().swap(m_dens);
	vector<double>().swap(beta_dens);
	vector<double>().swap(m_overlap);
	vector<int>().swap(occupied);
	vector<int>().swap(occupied_beta);
}
/***************************************************************************************/
bool Imolecule::check(){
	
	bool all_ok = true;
	
	m_log->input_message("Checking molecular information!");
	if ( num_of_atoms == 0 || atoms.size() == 0 ){
		m_log->input_message("Information of atoms not stored"); 
		all_ok = false;
	}
	if ( num_of_electrons == 0 ){
		m_log->input_message("Zero electrons counted for the molecule!"); 
		all_ok = false;
	}
	if ( num_of_ao == 0 ){
		m_log->input_message("Zero atomic orbitals counted for the molecule!"); 
		all_ok = false;
	}
	if ( coeff_MO.size() == 0 ){
		m_log->input_message("Molecular orbitals coeffcients not stored"); 
		all_ok = false;
	}
	if ( orb_energies.size() == 0 || MOnmb == 0) {
		m_log->input_message("Molecular orbitals energies not stored");  
		all_ok = false;
	}
	if ( m_overlap.size() == 0 || m_overlap[0] == 0 ){
		m_log->input_message("Overlap matrix coeffcients not stored.");  
		m_log->input_message("If the QM output is from MOPAC, orca or gamess your files may contain problems. ");  
	}
	if ( occupied.size() == 0 ){
		m_log->input_message("Molecular Orbital occupation information was not stored.");
		all_ok = false;
	}else if ( occupied[0] == 0){
		m_log->input_message("Molecular Orbital occupation information was not stored.");
		all_ok = false;
	}
	
	if ( homo_energy == 0.0000 || homoN == 0 ){
		m_log->input_message("HOMO energy may bee wrong!");
		all_ok = false;
	}
	if ( lumo_energy == 0.0000 || lumoN == 0 ){
		m_log->input_message("LUMO energy may bee wrong!");
		all_ok = false;
	}
	double mol_charge_test = 0;
	for( int i=0;i<atoms.size()/2;i++){
		mol_charge_test+=atoms[i].charge;
	}
	if ( mol_charge_test == 0.000 ){
		m_log->input_message("Charges stored in atoms objects are not consistent!");
		all_ok = false;
	}
	
	if ( betad ){
		if ( coeff_MO_beta.size() == 0 || MOnmb_beta == 0 ){
			m_log->input_message("Your QM output is said to contain beta set molecular orbitals, but they were not store.");
			all_ok = false;
		}
		if ( orb_energies_beta.size() == 0 ){
			m_log->input_message("Your QM output is said to contain beta set molecular orbitals, but their energies  were not store.");
			all_ok = false;
		}
	}	
	if ( !all_ok) {
		cout << "May have errors in the parsing processes of QM output files! " << endl;
		cout << "Please, run the program with -log/-verbose flag to search the source of possible errors. " << endl;
	}
	return all_ok;
}
/***************************************************************************************/
double Imolecule::check_ed(){
	double ed = 0.0;
	for (int i=0;i<MOnmb;i++){
		if ( occupied[i] > 0 )
			for ( int mu=0;mu<num_of_ao;mu++)
				for ( int nu=0;nu<num_of_ao;nu++)
					ed += coeff_MO[i*num_of_ao+mu]*
							coeff_MO[i*num_of_ao+nu];
	}
	//cout << "total electron density: " << ed << endl;
	ed = 0;
	for ( int nu=0;nu<num_of_ao;nu++){
		ed += coeff_MO[num_of_ao*1+nu]*
				coeff_MO[num_of_ao*1+nu];
	}
	return ed;
}
/***************************************************************************************/
void Imolecule::center_coord(){
	this->mol_vert_up();
	double m_x = ver_sup[0] - ver_inf[0];
	double m_y = ver_sup[1] - ver_inf[1];
	double m_z = ver_sup[2] - ver_inf[2];
	for(unsigned int i=0;i<num_of_atoms;i++){
		atoms[i].xcoord = atoms[i].xcoord - m_x;
		atoms[i].ycoord = atoms[i].ycoord - m_y;
		atoms[i].zcoord = atoms[i].zcoord - m_z;
	}
}
/***************************************************************************************/
Imolecule::~Imolecule(){
	//m_log->input_message("molecule going out of scope, name: ");
	//m_log->input_message(name);
}
/***************************************************************************************/
//================================================================================
//END OF FILE
//================================================================================