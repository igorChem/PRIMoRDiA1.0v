//Icube.cpp

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

//Including C++ headers
#include <iostream> 
#include <cstring>
#include <string> 
#include <sstream>
#include <vector> 
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <algorithm>

//Including PRIMoRDiA headers
#include "../include/common.h"
#include "../include/Ibuffer.h"
#include "../include/Iline.h"
#include "../include/log_class.h"
#include "../include/Imolecule.h"
#include "../include/Iatom.h"
#include "../include/Icube.h"

using std::vector;
using std::string;
using std::move;
using std::unique_ptr;
using std::cout;
using std::endl;

/***************************************************************************/
Icube::Icube()			:
	name("dummyname")	,
	elec_dens(true)		,
	MOn(0)				,
	voxelN(1)			,
	header("")			{
	
	for( int i=0; i<3; i++ ){
		origin[i]	= 0.0;
		gridsides[i]= 0.0;
		grid[i]		= 0.0;
	}
}
/***************************************************************************/
Icube::Icube(const char* file_nam)	:
	name(file_nam)					,
	elec_dens(true)					,
	MOn(0)							,
	voxelN(1)						,
	header("")						{
		
	m_log->input_message("Trying to open cube file.\n");
	
	if (!check_file_ext(".cube",file_nam)) { 
		cout << "Warning! The extension is not .cube \n";
		m_log->write_warning("The extension is not  .cube \n");
	}
	
	int nof = 0;
	
	name = name.substr(0,name.size()-4);
	
	if ( IF_file(file_nam) ){
		std::ifstream cube_file(file_nam);
		if( cube_file.is_open() ) {
			int line_len = 0;
			while( !cube_file.eof() ){
				string line;
				getline(cube_file,line);
				++line_len;
			}
			cube_file.clear();
			cube_file.seekg(0,std::ios::beg);
			
			m_log->input_message( std::to_string(line_len) );
			
			for(int i=0;i<line_len;i++){
				string line;
				double dummy;
				double temp;
				getline(cube_file,line);
				std::stringstream stream(line);
			
				if( i==0||i==1 ){
					header += line;
					header += "\n";
				}else if( i==2 ){
					stream  >> nof 	>> origin[0]  >> origin[1] >> origin[2];
				}else if( i==3 ){
					stream  >> grid[0]	>> gridsides[0];
				}else if( i==4 ){
					stream >> grid[1] >> dummy	>> gridsides[1]; 
				}else if( i==5 ){
					stream  >> grid[2] >> dummy  >> dummy >> gridsides[2];
				}else if( i>5 && i<(6+nof) ){
					Iatom  atom; 
					stream >> atom.atomicN   >> dummy   >> atom.xcoord  >> atom.ycoord  >> atom.zcoord;
					molecule.add_atom( atom );
				}else if( i>(6+nof) ){
					break;
				}
			}
			voxelN = grid[0]*grid[1]*grid[2];
			scalar.resize(voxelN);
			cube_file.clear();
			cube_file.seekg(0,std::ios::beg);
			unsigned int count = 0;
			for(unsigned int j=0;j<line_len;j++){
				std::string line;
				getline(cube_file,line);
				if ( j >= 6+molecule.num_of_atoms){
					double temp;
					std::stringstream stream(line);
					while( stream >> temp ){
						scalar[count++] = temp;
					}
				}
			}
		}else{
			cout << "File is not open, Icube instance not initialized" << endl;
		}
	}else{
		cout << "Error in openning the cube file: " << file_nam << endl;
	}
}
/***************************************************************************/
Icube::Icube(const Icube& rhs_cube)	:
		name(rhs_cube.name)				,
		elec_dens(rhs_cube.elec_dens)	,
		MOn(rhs_cube.MOn)				,
		voxelN(rhs_cube.voxelN)			,
		header(rhs_cube.header)			,
		molecule(rhs_cube.molecule)		,
		scalar(rhs_cube.scalar)			{
			
		for(int i=0;i<3;i++){
			origin[i]		= rhs_cube.origin[i];
			gridsides[i]	= rhs_cube.gridsides[i];
			grid[i]			= rhs_cube.grid[i];
		}
}
/***************************************************************************/
Icube& Icube::operator=(const Icube& rhs_cube){	
	if ( this != &rhs_cube ){
		name		= rhs_cube.name;
		elec_dens 	= rhs_cube.elec_dens;
		MOn			= rhs_cube.MOn;
		voxelN		= rhs_cube.voxelN;
		header		= rhs_cube.header;
		scalar		= rhs_cube.scalar;
		molecule	= rhs_cube.molecule;
		
		for(int i=0;i<3;i++){
			origin[i]    = rhs_cube.origin[i];
			gridsides[i] = rhs_cube.gridsides[i];
			grid[i]      = rhs_cube.grid[i];
		}
	}
	return *this;
}
/***************************************************************************/
Icube::Icube(Icube&& rhs_cube) noexcept	:
		name( move(rhs_cube.name) )		,
		elec_dens(rhs_cube.elec_dens)	,
		MOn(rhs_cube.MOn)				, 
		voxelN(rhs_cube.voxelN)			,
		header(rhs_cube.header)			,
		molecule( move(rhs_cube.molecule) ),
		scalar( move(rhs_cube.scalar) )	{

		for(int i=0;i<3;i++){
			origin[i]	= rhs_cube.origin[i];
			gridsides[i]= rhs_cube.gridsides[i];
			grid[i]		= rhs_cube.grid[i];
		}
}
/***************************************************************************/
Icube& Icube::operator=(Icube&& rhs_cube) noexcept {
	if ( this != &rhs_cube ){
		name		= move ( rhs_cube.name );
		elec_dens	= rhs_cube.elec_dens;
		MOn			= rhs_cube.MOn;
		voxelN		= rhs_cube.voxelN;
		header		= rhs_cube.header;
		scalar		= move(rhs_cube.scalar);
		molecule	= move(rhs_cube.molecule);
		
		for(int i=0;i<3;i++){
			origin[i]	= rhs_cube.origin[i];
			gridsides[i]= rhs_cube.gridsides[i];
			grid[i]		= rhs_cube.grid[i];
		}
	}
	return *this;
}
/***************************************************************************/
bool operator==(const Icube& lhs_cube,const Icube& rhs_cube){
	if ( &lhs_cube == &rhs_cube ){
		return true;
	}
	bool result = true;
	for(int i=0;i<3;i++){
		if ( lhs_cube.origin[i]		!= rhs_cube.origin[i]		) result =  false; 
		if ( lhs_cube.grid [i]		!= rhs_cube.grid[i]			) result =  false; 
		if ( lhs_cube.gridsides[i]	!= rhs_cube.gridsides[i]	) result =  false; 
	}
	if ( result == false ) {
		cout << lhs_cube.origin[0] << endl;
		cout << rhs_cube.origin[0] << endl;
	}
	return result;
}  
/***************************************************************************/
Icube operator-(const Icube& lhs_cube, const Icube& rhs_cube){
	Icube Result(lhs_cube);
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] -= rhs_cube.scalar[x]; }
	return Result;
} 
/***************************************************************************/
Icube operator+(const Icube& lhs_cube, const Icube& rhs_cube){	
	Icube Result(lhs_cube);
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] += rhs_cube.scalar[x]; }
	return Result;
}
/***************************************************************************/
Icube operator*(const Icube& lhs_cube, const  Icube& rhs_cube){
	Icube Result(lhs_cube);
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] *= rhs_cube.scalar[x] ;}
	return Result;
}
/***************************************************************************/
Icube operator/(const Icube& lhs_cube, const Icube& rhs_cube){
	Icube Result(lhs_cube);
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) {
		if ( rhs_cube.scalar[x] == 0.00 ) Result.scalar[x] = 1000.0;
		else Result.scalar[x] /= rhs_cube.scalar[x]; 
	}
	return Result;
}
/***************************************************************************/
Icube operator-(const Icube& lhs_cube, double value){
	Icube Result(lhs_cube);
	
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] -=  value; } 
	return Result;	
}
/***************************************************************************/
Icube operator+(const Icube& lhs_cube, double value){
	Icube Result(lhs_cube);
	
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] +=  value; }
	return Result;	
}
/***************************************************************************/
Icube operator*(const Icube& lhs_cube, double value){
	Icube Result(lhs_cube);
	for(unsigned int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] *=  value; }
	return Result;	
}
/***************************************************************************/
Icube operator/(const Icube& lhs_cube, double value){
	Icube Result(lhs_cube);
	for(int x=0;x<lhs_cube.voxelN;x++) { Result.scalar[x] /= value; }
	return Result;
}
/***************************************************************************/
Icube Icube::scale_cube(double val){
	Icube Result(*this);	
	for(int x=0;x<voxelN;x++) { Result.scalar[x] = pow(scalar[x],val); }
	return Result;
}
/***************************************************************************/
Icube Icube::SQ(){ 
	Icube Result(*this);
	for(int x=0;x<voxelN;x++) { Result.scalar[x] = this->scalar[x]*this->scalar[x]; }
	return Result;
}
/***************************************************************************/
double Icube::calc_cube_integral(){
	double integral = 0;
	for (int i=0;i<voxelN;i++) { integral += scalar[i]; }
	integral *= std::abs(gridsides[0]*gridsides[1]*gridsides[2]);
	return integral;
}
/***************************************************************************/
void Icube::normalize( double norm ){
	double inte = this->calc_cube_integral();
	for( int x=0; x<voxelN; x++ ) { 
		scalar[x] /= inte;
	}
	if ( norm > 0 ){
		for( int x=0; x<voxelN; x++ ) { 
			scalar[x] *= norm;
		}
	}
}
/***************************************************************************/
double Icube::diff_integral(const Icube& cube){
	Icube res (*this - cube ) ;
	res.write_cube(name +"_diff_" +cube.name);
	return res.calc_cube_integral();
}
/***************************************************************************/
double Icube::similarity_index(Icube& cube,string type){
	double result = 0.0;
	if ( type == "default" ){
		Icube res = *this;
		res = res*cube;	
		double res_int = res.calc_cube_integral();
		Icube cub1  = this->SQ();
		Icube cub2  = cube.SQ();
		double int1 = cub1.calc_cube_integral();
		double int2 = cub2.calc_cube_integral();
		result = res_int/sqrt(int1*int2);		
	}
	return result;
}
/***************************************************************************/
Icube Icube::log_cube(){
	Icube result(*this);
	omp_set_num_threads(NP);
	#pragma omp parallel for
	for(unsigned int i=0;i<voxelN;i++) { result.scalar[i] = log(this->scalar[i]);}
	return result;
}
/***************************************************************************/
Icube Icube::calculate_complement(bool clog ){
	Icube complement(*this);
	complement.name   = name + "complement_log";
	complement.header =" Cube with the scalar field for the electron density complement of the " 
					   + name +" electron density \n generated by PRIMoRDiA \n";
	for(unsigned int i=0;i<voxelN;i++) { complement.scalar[i] = 1.0000000000; }
	complement = complement - *this;
	if ( clog ) { complement = complement.log_cube(); }
	return complement;
}
/******************************************************************************************/
double Icube::get_scalar(int x, int y, int z){ return scalar[x*x*grid[0] + y*grid[1] + z]; }
/*******************************************************************************************/
double Icube::get_scalar(double x, double y, double z){	
	int xx = (x - origin[0])/gridsides[0];
	int yy = (y - origin[1])/gridsides[1];
	int zz = (z - origin[2])/gridsides[2];
	return get_scalar(xx,yy,zz);
}
/***************************************************************************/
void Icube::add_data(vector<vector<vector<double> > >& data){
	int count = 0;
	for (unsigned int i=0;i<grid[0];i++){ 
		for (unsigned int j=0;j<grid[1];j++){
			for (unsigned int k=0;k<grid[2];k++) scalar[count++] = data[i][j][k];
		}
	}
}
/***************************************************************************/
void Icube::write_cube(string cubeName){
	std::ofstream cube_file;
	cube_file.open(cubeName.c_str());
	cube_file.precision(6);
	cube_file << std::fixed;
	
	string cube_type = " ";
	if ( elec_dens ) cube_type = "Total Electronic Density ";
	else cube_type = "Molecular Orbital ";
	
	if ( header == "" ){
		header =  "Cube file written by Icube class created by Igor Barden Grillo igorChem on github\n";
		header +=  cube_type + "for "  + name + "\n";
	}
		
	cube_file << header
			  << std::setw(5)  << std::right << molecule.num_of_atoms
			  << "  " 
			  << std::setw(12) << std::right << origin[0] 
			  << "  " 
			  << std::setw(12) << std::right << origin[1] 
			  << "  " 
			  << std::setw(12) << std::right << origin[2] 
			  << "\n"
			  << std::setw(5)  << std::right << grid[0]
			  << "  "  
			  << std::setw(12) << std::right << gridsides[0] 
			  << "  "  
			  << std::setw(12) << std::right << " 0.000000" 
			  << "  "
			  << std::setw(12) << std::right << " 0.000000" 
			  << "\n"
			  << std::setw(5)  << std::right << grid[1]
			  << "  "  
			  << std::setw(12) << std::right << "0.000000" 
			  << "  " 
			  << std::setw(12) << std::right << gridsides[1] 
			  << "  "
			  << std::setw(12) << std::right << "0.000000"
			  << "\n"
			  << std::setw(5)  << std::right << grid[2]
			  << "  "  
			  << std::setw(12) << std::right << "0.000000" 
			  << "  "
			  << std::setw(12) << std::right << "0.000000" 
			  << "  "
			  << std::setw(12) << std::right << gridsides[2] 
			  << "\n";
	
	if ( elec_dens == false ) cube_file << "    1    " << MOn << "\n"; 

	for (int i=0; i<molecule.num_of_atoms;i++){
		cube_file << std::setw(5)  << std::right << molecule.atoms[i].atomicN 
				  << "   " 
				  << std::setw(4)  << std::right << molecule.atoms[i].atomicN
				  << "   " 
				  << std::setw(12) << std::right << molecule.atoms[i].xcoord
				  << "   " 
				  << std::setw(12) << std::right << molecule.atoms[i].ycoord
				  << "   " 
				  << std::setw(12) << std::right << molecule.atoms[i].zcoord
				  << "\n";
	}
	
	cube_file.precision(6);
	cube_file << std::scientific;
	
	for(unsigned int x=0;x<voxelN;x++){
		cube_file << std::setw(16) << scalar[x] << " " ;
		if ( x%6==5 ) cube_file << "\n";
	}	
	cube_file.close();
}
/***************************************************************************/
void Icube::get_cube_stats(double& mean, double& sum, double& min, double& max){
	mean	= mean_dvec(scalar);
	sum		= sum_dvec(scalar);
	min		= min_dvec(scalar);
	max		= max_dvec(scalar);
}
/***************************************************************************/
void Icube::print(){
	cout << name << endl;
	cout << molecule.num_of_atoms << endl;
	for(int i=0;i<3;i++){
		cout << "Origin for "			<< i << " dimension: " << origin[i] << endl; 
		cout << "Number of points for "	<< i << " dimension: " << grid[i]   << endl; 
		cout << "Gridsides for "		<< i << " dimension: " << endl; 
	}
	std::cout << "Printing the first fifty voxel values" << endl;
	for(int j=0;j<50;j++) cout << scalar[j] << endl;
}
//**************************************************/
Icube::~Icube(){

}

//**************************************************/
cube_diffs::cube_diffs(const char* name){	
	unique_ptr<Ibuffer> file ( new Ibuffer(name,true) );
	unsigned int size = stoi( file->lines[0].words[0] );
	for ( unsigned int i = 1; i<size+1; i++ ) {
		const char* cub1name = file->lines[i].words[0].c_str();
		const char* cub2name = file->lines[i].words[1].c_str();
		Icube cub1(cub1name);
		Icube cub2(cub2name);
		double val = cub1.similarity_index(cub2,"default");
		diffs.push_back(val);
		labels.push_back(cub1.name + "__" + cub2.name);
	}
}
//**************************************************/
void cube_diffs::write(){
	std::ofstream report("reportdiffs");
	for(unsigned int i = 0; i<diffs.size();i++){
		report << labels[i] << " " << diffs[i] << " \n";
	}
	report.close();
}
/****************************************************/
cube_diffs::~cube_diffs(){}
//================================================================================
//END OF FILE
//================================================================================
