//pair.cpp

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

#include <string> 
#include <vector>

#include "../include/pair_rd.h"
#include "../include/local_rd.h"
#include "../include/primordia.h"
#include "../include/Iatom.h"
#include "../include/Imolecule.h"
#include "../include/local_rd_cnd.h"


using std::move;
using std::vector;
using std::to_string;

/***********************************************************************/
pair_rd::pair_rd()	:
	label("H-H")	{	
}
/***********************************************************************/
pair_rd::pair_rd(vector<primordia>& lrd, int ax, int ay){
	
	double tmp_dist		= 0.000;
	double ct_tmp1		= 0.0;
	double ct_tmp2		= 0.0;
	unsigned int a1, a2 = 0;
	a1 = ax - 1;
	a2 = ay - 1;	
	
	label = lrd[0].mol_info.atoms[a1].element;
	label +=to_string(ax);
	label +="_";
	label +=lrd[0].mol_info.atoms[a2].element;
	label +=to_string(ay);
	
	rd_labels.push_back(label);
	rd_labels[0] += "_NCTP";
	rd_labels.push_back(label);
	rd_labels[1] += "_HPI_Vee";
	rd_labels.push_back(label);
	rd_labels[2] += "_HPI_LCP";
	rd_labels.push_back(label);
	rd_labels[3] += "_HPI_FP";
	rd_labels.push_back(label);
	rd_labels[4] += "_SPI";
	rd_labels.push_back(label);
	rd_labels[5] += "_EEP";
	
	CTP.resize( lrd.size() );
	HPI_A.resize( lrd.size() );
	HPI_B.resize( lrd.size() );
	HPI_C.resize( lrd.size() );
	SPI.resize( lrd.size() );
	EEP.resize( lrd.size() );

	for( int k=0; k<lrd.size(); k++ ){
		tmp_dist= lrd[k].mol_info.calc_dist( a1 , a2 );
		ct_tmp1 = lrd[k].lrdCnd.lrds[0][ a1 ] - lrd[k].lrdCnd.lrds[1][ a2 ];
		ct_tmp2	= lrd[k].lrdCnd.lrds[0][ a2 ] - lrd[k].lrdCnd.lrds[1][ a1 ];
		CTP[k]	= 2*( (ct_tmp1 - ct_tmp2)/tmp_dist );
		HPI_A[k]= lrd[k].lrdCnd.lrds[4][a1]*lrd[k].lrdCnd.lrds[4][a2];
		HPI_B[k]= lrd[k].lrdCnd.lrds[5][a1]*lrd[k].lrdCnd.lrds[5][a2];
		HPI_C[k]= lrd[k].lrdCnd.lrds[6][a1]*lrd[k].lrdCnd.lrds[6][a2];
		SPI[k]	= lrd[k].lrdCnd.lrds[17][a1]*lrd[k].lrdCnd.lrds[17][a2];
		EEP[k]	= lrd[k].lrdCnd.lrds[11][a1]*ct_tmp2 + lrd[k].lrdCnd.lrds[11][a2]*ct_tmp1;
	}
}
/***********************************************************************/
pair_rd::~pair_rd(){}
/***********************************************************************/
pair_rd::pair_rd(const pair_rd& rhs):
	a_indices(rhs.a_indices)		,
	label(rhs.label)				,
	rd_labels(rhs.rd_labels)		,
	CTP(rhs.CTP)					,
	HPI_A(rhs.HPI_A)				,
	HPI_B(rhs.HPI_B)				,
	HPI_C(rhs.HPI_C)				,
	SPI(rhs.SPI)					,
	EEP(rhs.EEP)					{
}
/***********************************************************************/
pair_rd& pair_rd::operator=(const pair_rd& rhs){
	if ( this != &rhs ){
		a_indices	= rhs.a_indices;
		label		= rhs.label;
		rd_labels	= rhs.rd_labels;
		CTP			= rhs.CTP;
		HPI_A		= rhs.HPI_A;
		HPI_B		= rhs.HPI_B;
		HPI_C		= rhs.HPI_C;
		SPI			= rhs.SPI;
		EEP			= rhs.EEP;
	}
	return *this;
}
/***********************************************************************/
pair_rd::pair_rd(pair_rd&& rhs) noexcept:
	a_indices( move(rhs.a_indices) )	,
	label( move(rhs.label) )			,
	rd_labels( move(rhs.rd_labels) )	,
	CTP( move(rhs.CTP) )				,
	HPI_A( move(rhs.HPI_A) )			,
	HPI_B( move(rhs.HPI_B) )			,
	HPI_C( move(rhs.HPI_C) )			,
	SPI( move(rhs.SPI) )				,
	EEP( move(rhs.EEP) )				{
	
}
/***********************************************************************/
pair_rd& pair_rd::operator=(pair_rd&& rhs) noexcept{
	if ( this != &rhs ){
		a_indices	= move(rhs.a_indices);
		label		= move(rhs.label);
		rd_labels	= move(rhs.rd_labels);
		CTP			= move(rhs.CTP);
		HPI_A		= move(rhs.HPI_A);
		HPI_B		= move(rhs.HPI_B);
		HPI_C		= move(rhs.HPI_C);
		SPI			= move(rhs.SPI);
		EEP			= move(rhs.EEP);
	}
	return *this;
}
//========================================================================