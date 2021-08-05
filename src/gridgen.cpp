//gridgen.cpp
// source file for the gridgen class

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
//---------------------------------------------
//including c++ headers
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
//---------------------------------------------
//including PRIMoRDiA headers 
#include "../include/log_class.h"
#include "../include/common.h"
#include "../include/Iaorbital.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/Icube.h"
#include "../include/gridgen.h" 
//-----------------------------------------------
using std::unique_ptr;
using std::move;
using std::string;
using std::cout;
using std::endl;
/***********************************************************************/
gridgen::gridgen()	:
	name("nonamed")	,
	orbital(false)	,
	Norb(0)			{ 
		
	for( int i=0; i<3; i++ ){
		origin[i]		= 0.0;
		top_corner[i]	= 0.0;
		grid_sides[i]	= 0.0;
		grid_len[i]		= 0;
	}
}
/***********************************************************************/
gridgen::gridgen(int grd					, 
				Imolecule&&	 mol)	noexcept:
	name(mol.name)							,
	orbital(false)							,
	Norb(0)									,
	molecule( move (mol) )					{
	
	unsigned int i,j,k;
		
	if ( !molecule.bohr ) molecule.ang_to_bohr();
	molecule.mol_vert_up();
	double ver_inf	= molecule.ver_inf[0];
	double ver_sup	= molecule.ver_sup[0];

	for ( i=0; i<3; i++){
		//cout << molecule.ver_inf[i] << " " << molecule.ver_sup[i] << endl;
		if ( ver_inf > molecule.ver_inf[i] )	ver_inf	= molecule.ver_inf[i];
		if ( ver_sup < molecule.ver_sup[i] )	ver_sup = molecule.ver_sup[i];
	}

	for( i=0;i<3;i++){
		origin[i]		= ver_inf - 5.0;
		top_corner[i]	= ver_sup + 5.0;
		grid_len[i]		= grd;
		grid_sides[i]	= ( top_corner[i] - origin[i] ) /(grid_len[i]-1);
	}
	
	for ( j=0; j<molecule.num_of_atoms; j++ ){
		for ( k=0; k<molecule.atoms[j].norb; k++ ){
 			orbs.emplace_back( molecule.atoms[j].orbitals[k] );
			AOxcoords.push_back( molecule.atoms[j].xcoord ); 
			AOycoords.push_back( molecule.atoms[j].ycoord );
			AOzcoords.push_back( molecule.atoms[j].zcoord );
		}
	}
	if ( orbs.size() <= molecule.num_of_atoms ){
		m_log->input_message("Your atomic basis may not be loaded properly.\nVerify your QM data and/or use -log option to run PRIMoRDiA\n.");
	}
	psi.resize(grid_len[0]);
	for ( j=0; j<grid_len[0]; j++){
		psi[j].resize(grid_len[1]);
		for( k =0; k<grid_len[1]; k++) psi[j][k].resize(grid_len[2]);
	}
	
	density.voxelN = grd*grd*grd;
	density.scalar.resize(density.voxelN);
	for ( i=0; i<3; i++){
		density.origin[i]	= origin[i];
		density.gridsides[i]= grid_sides[i];
		density.grid[i]		= grid_len[i];
	}
	density.molecule.atoms			= molecule.atoms;
	density.molecule.num_of_atoms	= molecule.num_of_atoms;
	density.molecule,name			= molecule.name;
}
/***********************************************************************/
void gridgen::set_lim(double* Min, double* gridSides, int *gridSize){
	for (int i=0;i<3;i++){
		origin[i]		= Min[i];
		grid_sides[i]	= gridSides[i];
		grid_len[i]		= gridSize[i];
		top_corner[i]	= grid_sides[i]*grid_len[i] + origin[i];
	}
	psi.resize(grid_len[0]);
	for (unsigned int j=0;j<grid_len[0];j++){
		psi[j].resize(grid_len[1]);
		for(unsigned int k=0;k<grid_len[1];k++) psi[j][k].resize(grid_len[2]);
	}
	density.voxelN = grid_len[0]*grid_len[1]*grid_len[2];
	density.scalar.resize(density.voxelN);
	for (int i=0;i<3;i++){
		density.origin[i]	= origin[i];
		density.gridsides[i]= grid_sides[i];
		density.grid[i]		= grid_len[i];
	}
}
/***********************************************************************/
void gridgen::redefine_lim(int atom,double size){
	double xc = molecule.atoms[atom-1].xcoord;
	double yc = molecule.atoms[atom-1].ycoord;
	double zc = molecule.atoms[atom-1].zcoord;
	
	origin[0] = xc - size/2;
	origin[1] = yc - size/2;
	origin[2] = zc - size/2;
	
	top_corner[0] = xc + size/2;
	top_corner[1] = yc + size/2;
	top_corner[2] = zc + size/2;
	
	for(int i=0;i<3;i++){
		grid_sides[i] = ( top_corner[i] - origin[i] ) /(grid_len[i]-1);
	}
	
	density.voxelN = grid_len[0]*grid_len[1]*grid_len[2];
	density.scalar.resize(density.voxelN);
	density.name = molecule.name;
	for (int i=0;i<3;i++){
		density.origin[i]	= origin[i];
		density.gridsides[i]= grid_sides[i];
		density.grid[i]		= grid_len[i];
	}
}
/***********************************************************************/
void gridgen::redefine_lim(double xc, double yc, double zc,int size){
	
	origin[0] = xc*1.889726 - size/2;
	origin[1] = yc*1.889726 - size/2;
	origin[2] = zc*1.889726 - size/2;
	
	top_corner[0] = xc*1.889726 + size/2;
	top_corner[1] = yc*1.889726 + size/2;
	top_corner[2] = zc*1.889726 + size/2;
	
	for(int i=0;i<3;i++){
		grid_sides[i] = ( top_corner[i] - origin[i] ) /(grid_len[i]-1);
	}
	
	density.voxelN = grid_len[0]*grid_len[1]*grid_len[2];
	density.scalar.resize(density.voxelN);	
	for (int i=0;i<3;i++){
		density.origin[i]    = origin[i];
		density.gridsides[i] = grid_sides[i];
		density.grid[i]      = grid_len[i];
	}
}
/***********************************************************************/
double gridgen::calc_slater_orb(int i, int x, int y, int z){
	double value	= 1e-14;
	double xx		= x*grid_sides[0] + origin[0] - AOxcoords[i];
	double yy		= y*grid_sides[1] + origin[1] - AOycoords[i];
	double zz		= z*grid_sides[2] + origin[2] - AOzcoords[i];
	double r		= sqrt( xx*xx + yy*yy + zz*zz ); 
	double prr		= orbs[i].alpha*r;
	if ( prr > 40 )	return value;
	else{
		prr			= exp(-prr);
		double xxx	= pow(xx,orbs[i].powx);
		double yyy	= pow(yy,orbs[i].powy);
		double zzz	= pow(zz,orbs[i].powz);
		int    ll	= orbs[i].powx + orbs[i].powy + orbs[i].powz;
		double rr	= pow(r,orbs[i].shell-1-ll);
		value		= orbs[i].n_factor*xxx*yyy*zzz*rr*prr;
		return value;
	}
}
/***********************************************************************/
double gridgen::calc_gauss_orb(int i, int x, int y, int z){
	double value= 0.0;
	double xi	= x*grid_sides[0] + origin[0] - AOxcoords[i];
	double yi	= y*grid_sides[1] + origin[1] - AOycoords[i];
	double zi	= z*grid_sides[2] + origin[2] - AOzcoords[i];
	double r	= xi*xi + yi*yi + zi*zi; 
	double xx	= pow(xi,orbs[i].powx);
	double yy	= pow(yi,orbs[i].powy);
	double zz	= pow(zi,orbs[i].powz);
	for(unsigned int k=0;k<orbs[i].gtos.size();k++){
		double dr	= r*orbs[i].gtos[k].exponent;
		dr			= exp(-dr);
		value		+= orbs[i].gtos[k].n_fact*xx*yy*zz*dr;
	}
	return value;
}
/***********************************************************************/
double gridgen::calc_orca_sphe(int i, int x, int y, int z){
	double value = 0.0;
	double xi	= x*grid_sides[0] + origin[0] - AOxcoords[i];
	double yi	= y*grid_sides[1] + origin[1] - AOycoords[i];
	double zi	= z*grid_sides[2] + origin[2] - AOzcoords[i];
	double r	= xi*xi + yi*yi + zi*zi; 
	double dg	= 0;
	double xx	= pow(xi,orbs[i].powx);
	double yy	= pow(yi,orbs[i].powy);
	double zz	= pow(zi,orbs[i].powz);
	dg =  xx*yy*zz;
	if		( orbs[i].symmetry == "D0"	)	dg = 3*zi*zi - r;
	else if	( orbs[i].symmetry == "D1p"	)	dg = xi*zi;
	else if	( orbs[i].symmetry == "D1n"	)	dg = yi*zi;
	else if	( orbs[i].symmetry == "D2p"	)	dg = xi*xi - yi*yi;
	else if	( orbs[i].symmetry == "D2n"	)	dg = xi*yi;
	else if	( orbs[i].symmetry == "f0"	)	dg = -3*xi*xi*zi - 3*yi*yi*zi + 2*zi*zi*zi;
	else if	( orbs[i].symmetry == "f1p"	)	dg = -xi*xi*xi - yi*yi*xi  + 4*zi*zi*xi;
	else if	( orbs[i].symmetry == "f1n"	)	dg = -xi*xi*yi - yi*yi*yi + 4*zi*zi*yi;
	else if	( orbs[i].symmetry == "f2p"	)	dg = -yi*yi*zi + xi*xi*zi;
	else if	( orbs[i].symmetry == "f2n"	)	dg = xi*yi*zi;
	else if	( orbs[i].symmetry == "f3p" )	dg = -xi*xi*x + 3*yi*yi*xi;
	else if	( orbs[i].symmetry == "f3n"	)	{ dg = -3*xi*xi*yi + yi*yi*yi; }
	for(unsigned int k=0;k<orbs[i].gtos.size();k++){
		double dr	= r*orbs[i].gtos[k].exponent;
		dr			= exp(-dr);
		value		+= orbs[i].gtos[k].n_fact*dg*dr;
	}
	return value;
}
/***********************************************************************/
double gridgen::calc_aorb(int i, int x, int y, int z){
	if ( orbs[i].gto ) return this->calc_gauss_orb(i,x,y,z);
	else  return this->calc_slater_orb(i,x,y,z);
}
/***********************************************************************/
double gridgen::calc_orb_voxel(int nmo,int x,int y,int z,bool beta){
	double orb_value= 0.0;
	unsigned int aos= orbs.size();
	if ( !beta ){
		for(unsigned int i=0;i<aos;i++){
			if ( molecule.coeff_MO[aos*nmo + i] > 1e-8 ) {
				orb_value += calc_aorb(i,x,y,z)*molecule.coeff_MO[aos*nmo + i];
			}
		}
	}else if ( beta ){
		for(unsigned int i=0;i<orbs.size();i++){
			if ( molecule.coeff_MO_beta[molecule.MOnmb_beta*nmo + i] > 1e-8 ) {
				orb_value += calc_aorb(i,x,y,z)*molecule.coeff_MO_beta[aos*nmo + i];
			}
		}
	}
	return orb_value;
}
/***********************************************************************/
double gridgen::calc_orb_voxel_orca(int nm,int x,int y,int z,bool beta){
	double orb_value = 0.0;
	unsigned int aos = orbs.size();
	if ( !beta ){
		for(unsigned int i=0;i<aos;i++){
			if ( molecule.coeff_MO[aos*nm + i] > 1e-8 ) {
				orb_value +=calc_orca_sphe (i,x,y,z)*molecule.coeff_MO[aos*nm + i];
			}
		}
	}else if ( beta ){
		for(unsigned int i=0;i<orbs.size();i++){
			if ( molecule.coeff_MO_beta[molecule.MOnmb_beta*nm + i] > 1e-8 ) {
				orb_value += calc_orca_sphe(i,x,y,z)*molecule.coeff_MO_beta[aos*nm + i];
			}
		}
	}
	return orb_value;
}
/***********************************************************************/
void gridgen::calculate_orb(int Nmo, bool beta){
	orbital = true;
	Norb    = Nmo;
	unsigned int x,y,z;
	omp_set_num_threads(NP);
	#pragma omp parallel for collapse(3) default(shared) private(x,y,z) 
	for (x=0;x<grid_len[0];x++){
		for (y=0;y<grid_len[1];y++){
			for (z=0;z<grid_len[2];z++) {
				psi[x][y][z] = calc_orb_voxel(Nmo,x,y,z,beta);
			}
		}
	}
	density.add_data(psi);
	density.name = name;
}
/***********************************************************************/
void gridgen::calculate_orb_orca(int Nmo, bool beta){
	orbital = true;
	Norb    = Nmo;
	unsigned int x,y,z;
	omp_set_num_threads(NP);
	#pragma omp parallel for collapse(3) default(shared) private(x,y,z) 
	for (x=0;x<grid_len[0];x++){
		for (y=0;y<grid_len[1];y++){
			for (z=0;z<grid_len[2];z++) {
				psi[x][y][z] = calc_orb_voxel_orca(Nmo,x,y,z,beta);
			}
		}
	}
	density.add_data(psi);
	density.name = name;
}

/***********************************************************************/
double gridgen::electron_density_mo(int x, int y, int z){
	double valuealfa	= 0.0;
	double valuebeta	= 0.0;
	double phiK			= 0.0;
	std::vector<double> phialpha;
	std::vector<double> phibeta;
	unsigned int i,j,k,l,m,MOn,MOnb,AOn;
	MOn		= molecule.MOnmb;
	MOnb	= molecule.MOnmb_beta;
	phialpha.resize(MOn);
	phibeta.resize(MOnb);
	AOn		= orbs.size();

	for (i=0;i<AOn;i++){
		phiK = calc_aorb(i,x,y,z);
		for (j=0;j<MOn;j++){
			if( molecule.occupied[j] > 1e-08 ) {
				phialpha[j] += molecule.coeff_MO[AOn*j+i]*phiK;
			}
		}
		for (k=0;k<MOnb;k++){
			if( molecule.occupied[k] > 1e-08 ) {
				phibeta[k] += molecule.coeff_MO_beta[AOn*k+i]*phiK;
			}
		}
	}
	for (l=0;l<MOn;l++){
		valuealfa += molecule.occupied[l]*phialpha[l]*phialpha[l];
	}
	for (m=0;m<MOnb;m++){
		valuebeta += molecule.occupied_beta[m]*phibeta[m]*phibeta[m];
	}
	return valuealfa + valuebeta;
}
/***********************************************************************/
double gridgen::electron_density_mo_orca(int x, int y, int z){
	double valuealfa	= 0.0;
	double valuebeta	= 0.0;
	double phiK			= 0.0;
	std::vector<double> phialpha;
	std::vector<double> phibeta;
	unsigned int i,j,k,l,m,MOn,MOnb,AOn;
	MOn		= molecule.MOnmb;
	MOnb	= molecule.MOnmb_beta;
	phialpha.resize(MOn);
	phibeta.resize(molecule.MOnmb_beta);

	AOn		= orbs.size();

	for (i=0;i<AOn;i++){
		phiK = calc_orca_sphe(i,x,y,z);
		for (j=0;j<MOn;j++){
			if( molecule.occupied[j] > 1e-08 ) {
				phialpha[j] += molecule.coeff_MO[AOn*j+i]*phiK;
			}
		}
		for (k=0;k<MOnb;k++){
			if( molecule.occupied[k] > 1e-08 ) {
				phibeta[k] += molecule.coeff_MO_beta[AOn*k+i]*phiK;
			}
		}
	}
	for (l=0;l<MOn;l++){
		valuealfa += molecule.occupied[l]*phialpha[l]*phialpha[l];
	}
	for (m=0;m<MOnb;m++){
		valuebeta += molecule.occupied_beta[m]*phibeta[m]*phibeta[m];
	}
	return valuealfa + valuebeta;
}
/***********************************************************************/
double gridgen::electron_density(int x,int y, int z){
	int k 			= 0;
	double value = 0.0;
	unsigned int mu,nu;
	double occ	= 2.0;
	double occ2  = 1.0;
	
	if (molecule.betad) occ		= 1.0;
	if (molecule.betad) occ2	= 0.5;
	
	for(mu = 0;mu<molecule.MOnmb;mu++){
		for(nu = 0;nu<=mu;nu++){
			if ( mu == nu) value += occ2*calc_aorb(mu,x,y,z)*calc_aorb(nu,x,y,z)*molecule.m_dens[k++];
			else           value += occ*calc_aorb(mu,x,y,z)*calc_aorb(nu,x,y,z)*molecule.m_dens[k]*molecule.m_overlap[k++];
		}
	}
	
	for(mu = 0;mu<molecule.MOnmb_beta;mu++){
		for(nu =0;nu<=mu;nu++){
			if ( mu == nu) value += occ2*calc_aorb(mu,x,y,z)*calc_aorb(nu,x,y,z)*molecule.beta_dens[k++];
			else           value += occ*calc_aorb(mu,x,y,z)*calc_aorb(nu,x,y,z)*molecule.beta_dens[k]*molecule.m_overlap[k++];
		}
	}
	
	return value;
}
/***********************************************************************/
void gridgen::calculate_density(){
	unsigned int x,y,z;
	omp_set_num_threads(NP);
	//chronometer.reset();
	#pragma omp parallel for collapse (3) default(shared) private(x,y,z) 
	for (unsigned int x=0;x<grid_len[0];x++){
		for (unsigned int y=0;y<grid_len[1];y++){
			for (unsigned int z=0;z<grid_len[2];z++){
				psi[x][y][z] = electron_density_mo(x,y,z);
			}
		}
	} 
	//cout << "Execution time: " << chronometer.return_wall_time() << " seconds" << endl;
	density.add_data(psi);
	density.name = name;
}
/***********************************************************************/
void gridgen::calculate_density_orca(){
	unsigned int x,y,z;
	omp_set_num_threads(NP);
	//chronometer.reset();
	#pragma omp parallel for collapse (3) default(shared) private(x,y,z) 
	for (unsigned int x=0;x<grid_len[0];x++){
		for (unsigned int y=0;y<grid_len[1];y++){
			for (unsigned int z=0;z<grid_len[2];z++){
				psi[x][y][z] = electron_density_mo_orca(x,y,z);
			}
		}
	} 
	//cout << "Execution time: " << chronometer.return_wall_time() << " seconds" << endl;
	density.add_data(psi);
	density.name = name;
}
/***********************************************************************/
void gridgen::calculate_mep_from_charges(){
	unsigned int x,y,z,i;
	double xi,yi,zi,xj,yj,zj,r,invR,v = 0.0;
	double precision 				= 1e-13;
	
	omp_set_num_threads(NP);
	#pragma omp parallel for collapse(3) shared(precision) private(xi,yi,zi,xj,yj,zj,r,invR,x,y,z,i) reduction(+:v)
	for ( x=0;x<grid_len[0];x++ ){
		for ( y=0;y<grid_len[1];y++ ){
			for ( z=0;z<grid_len[2];z++ ){
				xi	=	x*grid_sides[0] + origin[0];
				yi	=	y*grid_sides[1] + origin[1];
				zi	=	z*grid_sides[2] + origin[2] ;
				for ( i=0;i<molecule.atoms.size();i++ ){
					xj	=	xi -	molecule.atoms[i].xcoord;
					yj	=	yi -	molecule.atoms[i].ycoord;
					zj	=	zi -	molecule.atoms[i].zcoord;
					r	=	sqrt(xj*xj +  yj*yj + zj*zj);
					invR	= 1/(r+precision);
					v	+= invR*molecule.atoms[i].charge;
				}
				psi[x][y][z] = v;
				v= 0;
			}
		}
	} 
	
	density.add_data(psi);
	density.name = name;
}
/***********************************************************************/
Icube& gridgen::get_cube(){ return density; }
/***********************************************************************/
Icube& gridgen::calc_HOMO(){
	if ( orbs[0].spherical ) this->calculate_orb_orca(molecule.homoN,false);
	else this->calculate_orb(molecule.homoN,false);
	return density;
}
/***********************************************************************/
Icube& gridgen::calc_LUMO(){
	if ( orbs[0].spherical ) this->calculate_orb_orca(molecule.lumoN,false);
	else this->calculate_orb(molecule.lumoN,false);
	return density;
}
/***********************************************************************/
Icube& gridgen::calc_HOMO_density(){
	if ( orbs[0].spherical ) this->calculate_orb_orca(molecule.homoN,false);
	else { this->calculate_orb(molecule.homoN,false); }
	density = move( density.SQ() );
	return density;
}
/***********************************************************************/
Icube& gridgen::calc_LUMO_density(){
	if ( orbs[0].spherical ) this->calculate_orb_orca(molecule.lumoN,false);
	else { this->calculate_orb(molecule.lumoN,false); }
	density = move( density.SQ() );
	return density;
}
/***********************************************************************/
void gridgen::write_grid(){
	string cb_name;
	string cb_name2;
	if ( orbital ) {
		std::string nmo = std::to_string(Norb);
		cb_name  = molecule.name + "_ph1_" + nmo + ".cube";
		cb_name2 = molecule.name + "_ph2_" + nmo + ".cube";
		density.write_cube(cb_name);
		density.write_cube(cb_name2);
	}else{
		cb_name =  molecule.name + ".cube"; 
		density.write_cube(cb_name);
	}
}
/***********************************************************************/
Icube gridgen::calc_HOMO_band(int bandn){
	Icube temp(density);
	for ( int i=(molecule.homoN-bandn);i<=molecule.homoN;i++){
		if ( molecule.orb_energies[i] >= (molecule.homo_energy-energy_crit) ){
			this->calculate_orb(i,false);
			temp = temp + density;
		}
	}
	return temp;
}
/***********************************************************************/
Icube gridgen::calc_LUMO_band(int bandn){
	Icube temp(density);
	for ( int i=molecule.lumoN;i<=molecule.lumoN+bandn;i++){
		if ( molecule.orb_energies[i] <= (molecule.lumo_energy+energy_crit) ){
			this->calculate_orb(i,false);
			temp = temp + density;
		}
	}
	return temp;
}
/***********************************************************************/
Icube gridgen::calc_band_EAS(int bandn){
	Icube temp = density;
	temp = temp*0.0;
	int cnt = 0;
	for ( int i=(molecule.homoN-bandn+1);i<=molecule.homoN;i++){
		if ( molecule.orb_energies[i] >= (molecule.homo_energy-energy_crit) ){
			this->calculate_orb(i,false);
			temp = temp + density.SQ();
			cnt++;
		}
	}
	temp = temp/cnt;
	//cout << "used " << cnt << " MO" << endl;
	return temp;
}
/***********************************************************************/
Icube gridgen::calc_band_NAS(int bandn){
	Icube temp(density);
	temp = temp*0.0;
	int cnt = 0;
	for ( int i=molecule.lumoN;i<molecule.lumoN+bandn;i++){
		if ( molecule.orb_energies[i] <= (molecule.lumo_energy+energy_crit) ){
			this->calculate_orb(i,false);
			temp = temp + density.SQ();
			cnt++;
		}
	}
	temp = temp/cnt;
	return temp;
}
/***********************************************************************/
Icube gridgen::calc_EBLC_EAS(){
	Icube temp = density;
	double coefficient = 0.0;
	temp 	= temp*0.0;
	for ( int i=0;i<=molecule.homoN;i++){
		coefficient = exp(-abs(molecule.orb_energies[i]-molecule.homo_energy ) );
		if ( coefficient > 0.36 ){
			this->calculate_orb(i,false);
			temp = temp + density.SQ()*coefficient;
		}
	}
	return temp;
} 
/***********************************************************************/
Icube gridgen::calc_EBLC_NAS( ){
	Icube temp = density;
	double coefficient = 0.0;
	temp = temp*0.0;
	for ( int i=molecule.lumoN;i<=molecule.orb_energies.size();i++){
		coefficient = exp(-abs(molecule.orb_energies[i]-molecule.lumo_energy ) );
		if ( coefficient > 0.36 ){
			this->calculate_orb(i,false);
			temp = temp + density.SQ()*coefficient;
			//cnt++;
		}
	}
	return temp;
} 
/***********************************************************************/
Icube& gridgen::grid_from_atoms(std::vector<double> values){
	unsigned int x,y,z,i;
	double xi,yi,zi,xj,yj,zj,r,invR,v = 0.0;
	double precision 				= 1e-13;
	
	if ( values.size() == molecule.num_of_atoms ){
		omp_set_num_threads(NP);
		#pragma omp parallel for collapse(3) shared(precision) private(xi,yi,zi,xj,yj,zj,r,invR,x,y,z,i) reduction(+:v)
		for ( x=0;x<grid_len[0];x++ ){
			for ( y=0;y<grid_len[1];y++ ){
				for ( z=0;z<grid_len[2];z++ ){
					xi	=	x*grid_sides[0] + origin[0];
					yi	=	y*grid_sides[1] + origin[1];
					zi	=	z*grid_sides[2] + origin[2];
					for ( i=0;i<molecule.atoms.size();i++ ){
						xj	=	xi -molecule.atoms[i].xcoord;
						yj	=	yi -molecule.atoms[i].ycoord;
						zj	=	zi -molecule.atoms[i].zcoord;
						r	=	sqrt(xj*xj +  yj*yj + zj*zj);
						invR	= 1/(r+precision);
						v	+= invR*values[i];
					}
					psi[x][y][z] = v;
					v= 0;
				}
			}
		} 
		density.add_data(psi);
		density.name = name;
	}
	return density;
}
/***********************************************************************/
gridgen::~gridgen(){}

//================================================================================
//END OF FILE
//================================================================================