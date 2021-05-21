// source file for the Iprotein and Iresidue class 
// Iprotein.cpp

// include statements from c++ library
#include <iostream>
#include <string> 
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <map>
// include statements from PRIMORDiA-libs
#include "../include/common.h"
#include "../include/Iprotein.h"
#include "../include/Ibuffer.h"
#include "../include/Iline.h"
#include "../include/Imolecule.h"
#include "../include/Iatom.h"
//---------------------------------------
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::move;
using std::stod;

std::map<string, double> res_vol = {
	{"ALA", 88.6}, {"ARG",173.4},  {"ASN",114.1 }, {"ASP", 111.1}, 
	{"CYS",108.5 }, {"GLN",143.8}, {"GLU",138.4 }, {"GLY",60.1 }, 
	{"HIS",153.2 }, {"ILE",166.7}, {"LEU",166.7},  {"LYS",168.6 },
	{"MET",162.9 }, {"PHE",189.9}, {"PRO",112.7 }, {"SER", 89.0}, 
	{"THR",116.1 }, {"TRP",227.8}, {"TYR",193.6 }, {"VAL", 140.0},
	{"HID",153.2 }, {"HIE",153.2}, {"HIP",153.2 }, {"UKW",0.001}
};

/*****************************************************/
Iresidue::Iresidue()	:
	name("none")		,
	type("UKW")			,
	ligand(false)		,
	molar_vol(0.001)	,
	atom_s(0)			{
}
/*****************************************************/
Iresidue::Iresidue(const Iresidue& res)		:
	name(res.name)							,
	type(res.type)							,
	ligand(false)							,
	atom_type(res.atom_type)				,
	molar_vol(res.molar_vol)				,
	atom_s(res.atom_s)						,
	atoms_index(res.atoms_index)			,
	side_chain_atoms(res.side_chain_atoms)	,
	back_bone_atoms(res.back_bone_atoms)	{
}
/*****************************************************/
Iresidue::Iresidue(Iresidue&& res) noexcept		:
	name( move(res.name) )						,
	type( move(res.type) )						,
	ligand(false)								,
	atom_type( move(res.atom_type) )			,
	molar_vol(res.molar_vol)					,
	atom_s(res.atom_s)							,
	atoms_index( move(res.atoms_index) )		,
	side_chain_atoms( move(res.side_chain_atoms) ),
	back_bone_atoms( move(res.back_bone_atoms) ){
	
}
/*****************************************************/
Iresidue& Iresidue::operator=(const Iresidue& res){
	if ( this!=&res){
		name			= res.name;
		type			= res.type;
		ligand			= res.ligand;
		atom_type		= res.atom_type;
		molar_vol		= res.molar_vol;
		atom_s			= res.atom_s;
		atoms_index		= res.atoms_index;
		side_chain_atoms= res.side_chain_atoms;
		back_bone_atoms	= res.back_bone_atoms;
	}
	return *this;
}
/*****************************************************/
Iresidue& Iresidue::operator=(Iresidue&& res) noexcept{
	if ( this!=&res){
		name			= move(res.name);
		type			= move(res.type);
		ligand			= res.ligand;
		atom_type		= move(res.atom_type);   
		molar_vol		= res.molar_vol;
		atom_s			= res.atom_s;                   
		atoms_index		= move(res.atoms_index);        
		side_chain_atoms= move(res.side_chain_atoms);
		back_bone_atoms = move(res.back_bone_atoms);
	}
	return *this;
}
/*****************************************************/
Iresidue::~Iresidue(){}
/*****************************************************/
Iprotein::Iprotein()	:
	name("none")		,
	ligand(false)		,
	remark("none")		,
	title("none")		,
	num_of_res(0)		{
}
/*****************************************************/
Iprotein::Iprotein (const char* pdb_name):
	num_of_res(0)									,
	remark("none")									,
	title("none")										,
	ligand(false)										{
	
	if ( !check_file_ext(".pdb",pdb_name) )	{
		cout << "Warning! The file has wrong extension name!" << endl;
		m_log->write_warning("Warning! The file has wrong etension name!");
	}
	name = remove_extension(pdb_name);
	
	vector<string> resn;
	vector<string> atomtype;
	vector<string> element;
	vector<string> resname;
	
	string old_resi 	= "0";
	int cnt 			= 0;
	string temp 		= "";
	
	char pdb_line[600];
	
	if ( IF_file( pdb_name ) ){
		std::ifstream buf (pdb_name) ;
		while( !buf.eof() ){
			buf.getline(pdb_line,600);
			string word(pdb_line,0,6);
			if ( word == "ATOM  " || word == "HETATM") {
				atomtype.emplace_back( pdb_line,12,4  );
				resname.emplace_back( pdb_line,17,3 );
				resn.emplace_back( pdb_line,23,3  );
				string temp2d_a(pdb_line,31,6);
				xcoord.push_back( stod(temp2d_a) );
				string temp2d_b(pdb_line,39,6);
				ycoord.push_back( stod(temp2d_b) ) ;
				string temp2d_c(pdb_line,47,6);
				zcoord.push_back( stod(temp2d_c) ) ;
				if (  (unsigned)strlen(pdb_line) >85 ) {
					string word2(pdb_line,85,2);
					if ( word2 == "P1" ) 
						p1.push_back(cnt);
					else if ( word2 == "P2" ){
						p2.push_back(cnt);
					}
				}
				
				if ( resn[cnt] != old_resi ){
					Iresidue res;
					res.type = resname[cnt];
					if ( res.type == "LIG" || res.type == "MOL" ) { 
						ligand = true;
						res.ligand = true;
						lig_num = num_of_res;
						res.molar_vol = res_vol[res.type];
					}
					residues.emplace_back(res);
					num_of_res++;
					old_resi = resn[cnt];
				}
				residues[num_of_res-1].atom_type.push_back(atomtype[cnt]);
				residues[num_of_res-1].atoms_index.push_back(cnt);
				residues[num_of_res-1].atom_s++;
				cnt++;
			}
		}
	}else{
		m_log->write_error("Not possible to open pdb! Verify if the file is on the path inficated in the input file!");
		cout << "Can't open pdb!" << endl;
	}
}
/*****************************************************/
Iprotein::Iprotein(const Iprotein& prot_rhs):
	name(prot_rhs.name)						,
	num_of_res(prot_rhs.num_of_res)			,
	remark(prot_rhs.remark)					,
	title(prot_rhs.title)					,
	ligand(prot_rhs.ligand)					,
	residues(prot_rhs.residues)				,
	xcoord(prot_rhs.xcoord)					,
	ycoord(prot_rhs.ycoord)					,
	zcoord(prot_rhs.zcoord)					,
	b_factor(prot_rhs.b_factor)				,
	p1(prot_rhs.p1)							,
	p2(prot_rhs.p2)							{
}
/*****************************************************/
Iprotein::Iprotein(Iprotein&& prot_rhs) noexcept:
	name( move(prot_rhs.name) )					,
	num_of_res( move(prot_rhs.num_of_res) )		,
	ligand(prot_rhs.ligand)						,
	residues( move(prot_rhs.residues)  )		,
	remark( move(prot_rhs.remark) )				,
	title(	move(prot_rhs.title) )				,
	xcoord( move(prot_rhs.xcoord) )				,
	ycoord( move(prot_rhs.ycoord) )				,
	zcoord( move(prot_rhs.zcoord) )				,
	b_factor( move(prot_rhs.b_factor) )			,
	p1( move(prot_rhs.p1) )						,
	p2( move(prot_rhs.p2) )						{
}
/*****************************************************/
Iprotein& Iprotein::operator=(const Iprotein& prot_rhs){
	if ( this!=&prot_rhs){
		name		= prot_rhs.name;                    
		num_of_res 	= prot_rhs.num_of_res;  
		ligand		= prot_rhs.ligand;
		residues	= prot_rhs.residues;
		remark		= prot_rhs.remark;
		title		= prot_rhs.title;
		xcoord		= prot_rhs.xcoord;
		ycoord		= prot_rhs.ycoord;
		zcoord		= prot_rhs.zcoord;
		b_factor	= prot_rhs.b_factor;
		p1			= prot_rhs.p1;
		p2			= prot_rhs.p2;
	}
	return *this;
}
/*****************************************************/
Iprotein& Iprotein::operator=(Iprotein&& prot_rhs) noexcept{
	if ( this!=&prot_rhs){
		name		= move(prot_rhs.name);                    
		num_of_res	= move(prot_rhs.num_of_res);  
		residues	= move(prot_rhs.residues);
		remark		= move(prot_rhs.remark);
		title		= move(prot_rhs.title);
		xcoord		= move(prot_rhs.xcoord);
		ycoord		= move(prot_rhs.ycoord);
		zcoord		= move(prot_rhs.zcoord);
		b_factor	= move(prot_rhs.b_factor);
		p1 			= move(prot_rhs.p1);
		p2 			= move(prot_rhs.p2);
	}
	return *this;	
}
/*****************************************************/
void Iprotein::load_b_column( vector<double>& b_fact ){
	if (b_factor.size() > 0 ) { b_factor.clear(); }
	copy( b_fact.begin(),b_fact.end(),back_inserter(b_factor) );
}
/*****************************************************/
void Iprotein::print(){
	for( unsigned int i=0;i<residues.size();i++ ){
		cout << residues[i].type << endl;
		for ( unsigned int j=0;j<residues[i].atom_s;j++){
			cout << residues[i].atom_type[j] << endl;
		}
	}
}
/*****************************************************/
Iprotein::~Iprotein(){}
/*****************************************************/
//====================================
pdb::pdb()			:
	nModels(0)		,
	name("noname")	{
}
/*****************************************************/
pdb::pdb(const pdb& rhs_pdb)	:
	nModels(rhs_pdb.nModels)	,
	name(rhs_pdb.name)			,
	models(rhs_pdb.models)		{
}
/*****************************************************/
pdb::pdb(pdb&& rhs_pdb) noexcept	:
	nModels( move(rhs_pdb.nModels) ),
	name( move(rhs_pdb.name) )		,
	models( move(rhs_pdb.models) )	{
}
/*****************************************************/
pdb& pdb::operator=(pdb& rhs_pdb){
	if ( this!=&rhs_pdb ){
		nModels	= rhs_pdb.nModels;
		name	= rhs_pdb.name;
		models	= rhs_pdb.models;
	}
	return *this;
}
/*****************************************************/
pdb& pdb::operator=(pdb&& rhs_pdb) noexcept{
	if ( this!=&rhs_pdb ){
		nModels	= rhs_pdb.nModels;
		name	= move(rhs_pdb.name);
		models	= move(rhs_pdb.models);
	}
	return *this;
}
/*****************************************************/
pdb::pdb(const char* file_name){
	//para ler multipdbs
}
/*****************************************************/
pdb::pdb(Iprotein& prot){
	nModels++;
	models.emplace_back(prot);
	name = prot.name +".pdb";
}
/*****************************************************/
void pdb::add_protein(Iprotein& prot){
	nModels++;
	models.emplace_back(prot);
}
/*****************************************************/
Iprotein& pdb::get_model(unsigned int i){
	return models[i];
}
/*****************************************************/
void pdb::write_models(string path){
	unsigned int i,j,k;
	for ( k=0; k<models.size(); k++){
		std::ofstream pdb_file;
		string fname = path + "/"+ models[k].title + ".pdb"; 
		pdb_file.open( fname.c_str() );
		pdb_file << std::fixed;
		pdb_file.precision(3);
		
		int cont = 0;
		for( i=0; i<models[k].residues.size(); i++ ){
			for ( j=0; j<models[k].residues[i].atom_s; j++ ){
				pdb_file	<< std::setw(6) << std::left  << "ATOM" 
							<< " "
							<< std::setw(4) << std::right  << (cont+1) 
							<< " "
							<< std::setw(4) << models[k].residues[i].atom_type[j] 
							<< " "
							<< std::left << std::setw(4) << models[k].residues[i].type 
							<< " "
							<< std::right << std::setw(4) << (i+1)
							<< std::setw(5) << " "
							<< std::setw(7) << models[k].xcoord[cont] 
							<< " "
							<< std::setw(7) << models[k].ycoord[cont] 
							<< " "
							<< std::setw(7) << models[k].zcoord[cont] 
							<< " "
							<< std::setw(5)  << "1.00"
							<< " "
							<< std::setw(5) << models[k].b_factor[cont]
							<< "\n";
							cont++;
			}
		}
		pdb_file.close();
	}
}
/*****************************************************/
void pdb::write_pdb(string fname){
	
	std::ofstream pdb_file;
	pdb_file.open( fname.c_str() );
	pdb_file << std::fixed;
	pdb_file.precision(3);
	
	for (unsigned int k=0;k<models.size();k++){
		pdb_file << "REMARK  " << models[k].remark << " \n"
				 << "TITLE " << models[k].title << " \n"
				 << "MODEL " << k << endl;
		int cont = 0;
		for( unsigned int i=0;i<models[k].residues.size();i++ ){
			for ( unsigned int j=0;j<models[k].residues[i].atom_s;j++){
				pdb_file	<< std::setw(6) << std::left  << "ATOM" 
							<< " "
							<< std::setw(4) << std::right  << (cont+1) 
							<< " "
							<< std::setw(4) << models[k].residues[i].atom_type[j] 
							<< " "
							<< std::left << std::setw(4) << models[k].residues[i].type 
							<< " "
							<< std::right << std::setw(4) << (i+1)
							<< std::setw(5) << " "
							<< std::setw(7) << models[k].xcoord[cont] 
							<< " "
							<< std::setw(7) << models[k].ycoord[cont] 
							<< " "
							<< std::setw(7) << models[k].zcoord[cont] 
							<< " "
							<< std::setw(5)  << "1.00"
							<< " "
							<< std::setw(5) << models[k].b_factor[cont]
							<< "\n";
							cont++;
			}
		}
		pdb_file << "ENDMDL" << endl;
	}
	pdb_file.close();
}
/*****************************************************/