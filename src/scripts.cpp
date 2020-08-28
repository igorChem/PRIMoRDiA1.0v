//scripts.cpp

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
#include <cstring>
#include <string>
#include <fstream>

#include "../include/common.h"
#include "../include/log_class.h"
#include "../include/scripts.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;

/**************************************************/
scripts::scripts()		:
	file_name("noname")	,
	s_type("notype")	{

}
/**************************************************/
scripts::scripts(const char* Nm):
	file_name(Nm)				,
	s_type("notype")			{
	
}
/**************************************************/
scripts::~scripts(){
	
}
/**************************************************/
void scripts::write_r_dos(vector<double>& energies){
	string Name =  file_name;
	Name 		+= ".DOS";
	std::ofstream dos_file( Name.c_str() );
	std::ofstream dos_file_gnu( (Name +".R").c_str() );
	dos_file << "Energy\n";
	for(int i=0;i<energies.size();i++) {dos_file << energies[i] << endl; }
	
	m_log->input_message("Outputing Density of States information to file.");
	
	dos_file_gnu << "require(ggplot2) \n"
				 << "dos = read.table( '" << Name << "',header=T)\n"
				 << "attach(dos) \n"
				 << "p <-ggplot(dos, aes( x=Energy) )+\n"
				 << "geom_density(fill='blue',bw=1) \n"
				 << "png('"<< Name << ".png',width = 5, height = 3.5, units ='"  << "in', res = 400)\n"
				 << "p\n dev.off()";
	dos_file.close();
	dos_file_gnu.close();
}
/**************************************************/
void scripts::write_pymol(){
	
}
//================================================================================
//END OF FILE
//================================================================================
	/*
	if ( M_R ){
		r_script_lrd.open((name+".R").c_str() );
		r_script_lrd << "#!/usr/bin/env Rscript\n"
					 << "library(pheatmap)\n"
					 << "data1 <-read.table('" << names << "',header=T)\n"
					 << "attach(data1)\n"
					 << "ras_s <-scale(RAS,scale=T,center=F)\n"
					 << "hard_s <-scale(Hardness,scale=T,center=F)\n"
					 << "d = ras_s + hard_s\n"
					 << "d_vec <-c(d)\n"
					 << "data2 <-data.frame(data1,d_vec)\n"
					 << "detach(data1)\n"
					 << "attach(data2)\n"
					 << "data3 <-subset(data2,d_vec>1.5*mean(d_vec))\n"
					 << "detach(data2)\n"
					 << "attach(data3)\n"
					 << "pro1 <-data.matrix(data3[2:7])\n"
					 << "pro1 <-scale(pro1,scale=T,center=F)\n"
					 << "rownames(pro1) <-res\n"
					 << "detach(data3)\n"
					 << "pheatmap(pro1,color = colorRampPalette(c('navy','white','red'))(100),border_color=NA,cluster_cols=F,fontsize=8,main='\n',filenam='"
					 << names << "heatmap.png'"
					 << ",width=4.3,height=5.3)\n"	
					 << "data1 <-read.table('" << names1 << "',header=T)\n"
					 << "attach(data1)\n"
					 << "ras_s <-scale(RAS,scale=T,center=F)\n"
					 << "hard_s <-scale(Hardness,scale=T,center=F)\n"
					 << "d = ras_s + hard_s\n"
					 << "d_vec <-c(d)\n"
					 << "data2 <-data.frame(data1,d_vec)\n"
					 << "detach(data1)\n"
					 << "attach(data2)\n"
					 << "data3 <-subset(data2,d_vec>1.5*mean(d_vec))\n"
					 << "detach(data2)\n"
					 << "attach(data3)\n"
					 << "pro1 <-data.matrix(data3[2:7])\n"
					 << "pro1 <-scale(pro1,scale=T,center=F)\n"
					 << "rownames(pro1) <-res\n"
					 << "detach(data3)\n"
					 << "pheatmap(pro1,color = colorRampPalette(c('navy','white','red'))(100),border_color=NA,cluster_cols=F,fontsize=8,main='\n',filenam='"
					 << names1 << "heatmap.png'"
					 << ",width=4.3,height=5.3)\n"
					 << "data1 <-read.table('" << names2 << "',header=T)\n"
					 << "attach(data1)\n"
					 << "ras_s <-scale(RAS,scale=T,center=F)\n"
					 << "hard_s <-scale(Hardness,scale=T,center=F)\n"
					 << "d = ras_s + hard_s\n"
					 << "d_vec <-c(d)\n"
					 << "data2 <-data.frame(data1,d_vec)\n"
					 << "detach(data1)\n"
					 << "attach(data2)\n"
					 << "data3 <-subset(data2,d_vec>0.8*mean(d_vec))\n"
					 << "detach(data2)\n"
					 << "attach(data3)\n"
					 << "pro1 <-data.matrix(data3[2:7])\n"
					 << "pro1 <-scale(pro1,scale=T,center=F)\n"
					 << "rownames(pro1) <-res\n"
					 << "pheatmap(pro1,color = colorRampPalette(c('navy','white','red'))(100),border_color=NA,cluster_cols=F,fontsize=8,main='\n',filenam='"
					 << names2 << "heatmap.png'"
					 << ",width=4.3,height=5.3)\n";
		r_script_lrd.close();	
	}
	*/
