//Iaorbital.cpp

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

//inclding c++ headers
#include <iostream> 
#include <string>
#include <vector>
#include <cstring>
#define _USE_MATH_DEFINES
#include <cmath>
//including PRIMoRDiA headers
#include "../include/common.h"
#include "../include/Iaorbital.h"

const double m_pi = 3.14159265358979323846;

using std::move;
using std::string;
using std::cout;
using std::endl;

/*****************************************************************************/
Iprimitive::Iprimitive():
	n_fact(0.0)			,
	exponent(1.0)		,
	c_coef(0.0)			{
}
/*****************************************************************************/
Iprimitive::Iprimitive(double exp			,
					   double contrac_coeff):
	n_fact(0.0)								,
	exponent(exp)							,
	c_coef(contrac_coeff)					{
}
/*****************************************************************************/
Iprimitive::Iprimitive(const Iprimitive& prim_rhs):
	n_fact(prim_rhs.n_fact)							,
	exponent(prim_rhs.exponent)						,
	c_coef(prim_rhs.c_coef)							{
}
/*****************************************************************************/
Iprimitive::Iprimitive(Iprimitive&& prim_rhs) noexcept	:
	n_fact( move(prim_rhs.n_fact) )						,
	exponent( move(prim_rhs.exponent) )					,
	c_coef( move(prim_rhs.c_coef) )						{
}
/*****************************************************************************/
Iprimitive& Iprimitive::operator=(const Iprimitive& prim_rhs){
	if ( this == &prim_rhs ){
		n_fact		= prim_rhs.n_fact;
		exponent	= prim_rhs.exponent;
		c_coef		= prim_rhs.c_coef;
	}
	return *this;
}
/*****************************************************************************/
Iprimitive& Iprimitive::operator=(Iprimitive&& prim_rhs) noexcept{
	if ( this == &prim_rhs ){
		n_fact		= move(prim_rhs.n_fact);
		exponent	= move(prim_rhs.exponent);
		c_coef		= move(prim_rhs.c_coef);
	}
	return *this;
}
/*****************************************************************************/
void Iprimitive::print(){
	cout << "	exponent: "			<< exponent << "\n"
		 << "	c. coefficient: "	<< c_coef   << "\n"
		 << "	Normalization f.: " << n_fact   << "\n\n"; 
}
/*****************************************************************************/
Iprimitive::~Iprimitive(){}
/*****************************************************************************/
//END OF Iprimitive Class
///////////////////////////////////////////////////////////////////////////////
/*****************************************************************************/
Iaorbital::Iaorbital()	:
	shell(1)			,
	gto(false)			,
	spherical(false)	,
	symmetry("S")		,
	n_factor(0.0)		,
	alpha(0.0)			,
	powx(0)				,
	powy(0)				,
	powz(0)				{
}
/*****************************************************************************/
Iaorbital::Iaorbital(unsigned int level	,
							string sym	,
							double coef):
	shell(level)						,
	gto(false)							,
	symmetry( move(sym) )				,
	spherical(false)					,
	n_factor(0.0)						,
	alpha(coef)							,
	powx(0)								,
	powy(0)								,
	powz(0)								{
	
	if (symmetry == "S"){
		powx = 0; 
		powy = 0;
		powz = 0;
	}else if ( symmetry == "PX"){
		powx = 1; 
		powy = 0;
		powz = 0;
	}else if ( symmetry == "PY"){
		powx = 0; 
		powy = 1;
		powz = 0;
	}else if ( symmetry == "PZ"){
		powx = 0; 
		powy = 0;
		powz = 1;
	}else if (symmetry == "XX"){
		powx = 2;
		powy = 0;
		powz = 0;
	}else if (symmetry == "YX"){
		powx = 0;
		powy = 2;
		powz = 0;
	}else if (symmetry == "ZZ"){
		powx = 0;
		powy = 0;
		powz = 2;
	}else if (symmetry == "XZ"){
		powx = 1;
		powy = 0;
		powz = 1;
	}else if (symmetry == "YZ"){
		powx = 0;
		powy = 1;
		powz = 1;
	}else if (symmetry == "XY"){
		powx = 1;
		powy = 1;
		powz = 0;
	}
}
/*****************************************************************************/
Iaorbital::Iaorbital(const Iaorbital& rhs_orb)	:
		shell(rhs_orb.shell)					,
		gto(rhs_orb.gto)						,
		symmetry(rhs_orb.symmetry)				,
		spherical(rhs_orb.spherical)			,
		n_factor(rhs_orb.n_factor)				,
		alpha(rhs_orb.alpha)					, 
		powx(rhs_orb.powx)						,
		powy(rhs_orb.powy)						, 
		powz(rhs_orb.powz)						,
		gtos(rhs_orb.gtos)						{
}
/*****************************************************************************/
Iaorbital& Iaorbital::operator=(const Iaorbital& rhs_orb){
	if( this!=&rhs_orb ){
		shell		= rhs_orb.shell;
		gto			= rhs_orb.gto;
		symmetry	= rhs_orb.symmetry;
		spherical	= rhs_orb.spherical;
		n_factor	= rhs_orb.n_factor;
		alpha		= rhs_orb.alpha;
		powx		= rhs_orb.powx;
		powy		= rhs_orb.powy;
		powz		= rhs_orb.powz;
		gtos		= rhs_orb.gtos;
	}
	return *this;
}
/****************************************************************************************/
Iaorbital::Iaorbital(Iaorbital&& rhs_orb) noexcept :
	shell(rhs_orb.shell)							,
	gto(rhs_orb.gto)								,
	symmetry( move(rhs_orb.symmetry) )				,
	spherical( move(rhs_orb.spherical) )			,
	n_factor(rhs_orb.n_factor)						,
	alpha(rhs_orb.alpha)							,
	powx(rhs_orb.powx)								,
	powy(rhs_orb.powy)								,
	powz(rhs_orb.powz)								,
	gtos( move(rhs_orb.gtos) )						{
}
/****************************************************************************************/
Iaorbital& Iaorbital::operator=(Iaorbital&& rhs_orb) noexcept {
	if( this!=&rhs_orb ){
		shell		= rhs_orb.shell;
		gto			= rhs_orb.gto;
		symmetry	= move(rhs_orb.symmetry);
		spherical	= rhs_orb.spherical;
		n_factor	= rhs_orb.n_factor;
		alpha		= rhs_orb.alpha;
		powx		= rhs_orb.powx;
		powy		= rhs_orb.powy;
		powz		= rhs_orb.powz;
		gtos		= move(rhs_orb.gtos);
	}
	return *this;
}
/****************************************************************************************/
void Iaorbital::add_primitive(double expo, double contrac){
	gto = true;
	gtos.emplace_back(expo,contrac);
}
/****************************************************************************************/
void Iaorbital::add_primitive(Iprimitive gorb){
	gto = true;
	gtos.emplace_back( move(gorb) );
}
/****************************************************************************************/
bool Iaorbital::normalize(){
	//m_log->input_message("Normalizing atomic orbitals!");
	if (!gto){
		double ff = 1.000000;
		n_factor  = 0;
		if 		( symmetry == "PX" || symmetry == "PY" || symmetry == "PZ" ) { ff = 3.000000; }
		else if	( symmetry == "XX" || symmetry == "YY" || symmetry == "ZZ" ) { ff = 5.0000000/4.00000; } 
		else if	( symmetry == "XY" || symmetry == "YZ" || symmetry == "XZ" ) { ff = 15.000000; } 
		double part1 = pow(2.0*alpha,(shell + 0.50));
		double part2 = sqrt((ff/ (4.0*m_pi)));
		double part3 = sqrt( factorial(2.0*shell) );
		n_factor		= part1*part2/part3;
	}else{
		double tmp		= (8./(m_pi*m_pi*m_pi));
		double tmp1		= (2./(m_pi*m_pi*m_pi));
		double tmp_f0	= 4./(pow(15.,.5));
		double tmp_f1	= 4./(pow(5.,.5));
		double tmp_f2p	= 4.;
		double tmp_f2n	= 8.;
		double tmp_f3	= 4./(pow(3.,.5));
		
		for(unsigned int i=0;i<gtos.size();i++){
			if 		( symmetry == "S" ) { 
				gtos[i].n_fact = gtos[i].c_coef * pow(gtos[i].exponent,0.75) * 0.71270547;
			}else if ( symmetry == "PX" || symmetry == "PY" || symmetry == "PZ" ){ 
				gtos[i].n_fact =  gtos[i].c_coef * pow( gtos[i].exponent,1.25 ) * 1.425410941;
			}else if ( symmetry == "XX" || symmetry == "YY" || symmetry == "ZZ"){
				gtos[i].n_fact = gtos[i].c_coef * pow( gtos[i].exponent,1.75 ) * 1.645922781;
			}else if ( symmetry == "XY" || symmetry == "YZ" || symmetry == "XZ"){
				gtos[i].n_fact = gtos[i].c_coef * pow( gtos[i].exponent,1.75 ) * 2.850821881;
			}else if ( symmetry == "XXX" || symmetry == "YYY" || symmetry == "ZZZ" ){
				gtos[i].n_fact = gtos[i].c_coef * pow( gtos[i].exponent,2.25 ) * 1.47215808929;
			}else if ( symmetry=="XXY" || symmetry=="YYX" || symmetry=="YYZ" || symmetry=="XXZ" || symmetry=="XYY" || symmetry=="XZZ"){
				gtos[i].n_fact = gtos[i].c_coef * pow( gtos[i].exponent,2.25 ) * 3.2918455612;
			}else if ( symmetry == "XYZ"){
				gtos[i].n_fact = gtos[i].c_coef * pow( gtos[i].exponent,2.25 ) * 5.701643762;
			}else if ( symmetry == "D0") {
				gtos[i].n_fact = gtos[i].c_coef * pow(2048 * pow(gtos[i].exponent, 7.0)/(9.0 *m_pi*m_pi*m_pi),0.25);
			}else if ( symmetry == "D1p" || symmetry == "D1n" || symmetry == "D2n" ){
				gtos[i].n_fact = gtos[i].c_coef * pow(2048 * pow(gtos[i].exponent, 7.0)/(m_pi*m_pi*m_pi),0.25);
			}else if ( symmetry == "D2p" ){
				gtos[i].n_fact = gtos[i].c_coef * pow(128 * pow(gtos[i].exponent, 7.0)/(m_pi*m_pi*m_pi),0.25);
			}else if ( symmetry == "F7"){
				gtos[i].n_fact = gtos[i].c_coef * pow(gtos[i].exponent, 2.5) * 1.4721580892990935;
			}else if ( symmetry == "F0" ){
				gtos[i].n_fact = gtos[i].c_coef *  tmp_f0 * pow(gtos[i].exponent,2) * pow( (tmp * gtos[i].exponent),0.25 ); 
			}else if ( symmetry == "F1p" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f1 * pow(gtos[i].exponent,2)  * pow( (tmp1 * gtos[i].exponent ),0.25);
			}else if ( symmetry == "F1n" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f1 * pow(gtos[i].exponent,2)  * pow( (tmp1 * gtos[i].exponent ),0.25);
			}else if ( symmetry == "F2p" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f2p * pow(gtos[i].exponent,2)  * pow( (tmp * gtos[i].exponent),0.25);
			}else if ( symmetry == "F2n" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f2n * pow(gtos[i].exponent,2)  * pow( (tmp * gtos[i].exponent),0.25);
			}else if ( symmetry == "F3p" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f3 * pow(gtos[i].exponent,2)  * pow( (tmp1 * gtos[i].exponent ),0.25);
			}else if ( symmetry == "F3n" ){
				gtos[i].n_fact = gtos[i].c_coef * tmp_f3 * pow(gtos[i].exponent,2)  * pow( (tmp1 * gtos[i].exponent),0.25);
			}
		}
	}
	return true;
}
/****************************************************************************************/ 
void Iaorbital::print(){
	 std::cout 	<< "Shell level :"			<< shell	<< "  \n" 
				<< "orbital symmetry :"		<< symmetry	<< "  \n"
				<< "Normalization factor :"	<< n_factor	<< "  \n"
				<< "Zeta coefficient :"		<< alpha	<< "  \n"
				<< "x sym :"				<< powx		<< "  " 
				<< "y sym :"				<< powy		<< "  "
				<< "z sym :"				<< powz		<< "  \n"
				<< "gaussian  primitives: "	<< gtos.size() 	<< "\n"
				<< std::endl;
	for (int i=0;i<gtos.size();i++) gtos[i].print();
}
/*******************************************************************************************/
Iaorbital::~Iaorbital(){} 
//================================================================================
//END OF FILE
//================================================================================
