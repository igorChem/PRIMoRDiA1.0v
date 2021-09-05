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
#include "../include/local_rd.h"
#include "../include/Icube.h"
#include "../include/ReactionAnalysis.h"
#include "../include/reaction_coord.h"
#include "../include/pos_traj.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::to_string;

/************************************************************************/
scripts::scripts()		:
	file_name("noname")	,
	s_type("notype")	{
		
	script_file.open("default.R");

}
/************************************************************************/
scripts::scripts(	string Nm	,
					string _type 	):
	file_name(Nm)					,
	s_type(_type)					{
	
	string fname = file_name;
	fname += "_";
	fname += s_type;
	if ( s_type == "pymols" || s_type == "pymols_pdb" ){
		fname += ".pym";
	}else{
		fname += ".R";
	}
	script_file.open( fname.c_str() );
	
	if( s_type == "pymols" || s_type == "pymols_pdb" ){
		script_file << "preset.publication(selection='all')\n"
					<< "set sphere_scale, 0.2\n"
					<< "set bg_rgb, white \n"
					<< "set stick_radius, 0.18\n";
	}else{
		script_file << "#!/usr/bin/env Rscript\n";
	}	
}
/**********************************************************/
scripts::scripts(const scripts& rhs):
	file_name(rhs.file_name)     	,
	s_type(rhs.s_type)				{
	
}
/***********************************************************/
scripts& scripts::operator=(const scripts& rhs){
	if ( this != &rhs){
		file_name	= rhs.file_name;
		s_type		= rhs.s_type;
	}
	return *this;
}
/**************************************************************************/
scripts::~scripts(){
	script_file.close();
}
/*************************************************************************/
void scripts::write_r_dos(vector<double>& energies){
	
	string Name = file_name;
	Name 		+= ".DOS";
	std::ofstream dos_file( Name.c_str() );
	dos_file << "Energy\n";
	for( unsigned i=0; i<energies.size(); i++ ) { dos_file << energies[i] << endl; }
	
	m_log->input_message("Outputing Density of States information to file.");
	
	script_file	 << "require(ggplot2) \n"
				 << "dos = read.table( '" << Name << "',header=T)\n"
				 << "attach(dos) \n"
				 << "p <-ggplot(dos, aes( x=Energy) )+\n"
				 << "\tgeom_density(fill='blue',bw=1) \n"
				 << "png('"<< Name << ".png',width = 5, height = 3.5, units ='"  << "in', res = 400)\n"
				 << "p\ndev.off()";
}
/****************************************************************************/
void scripts::write_pymol_cube(local_rd& lrdVol){
	
	if ( lrdVol.lrds[5].voxelN > 0 ){
		string pdb_name = file_name + ".pdb";
		std::string typestr;
		std::string typestr2;
		if ( lrdVol.FD ){ 
			typestr = lrdVol.name+"FD_";
			typestr2 = get_file_name( typestr.c_str() );
		}
		else{ 
			typestr = lrdVol.name+"FOA_";
			typestr2 = get_file_name( typestr.c_str() );
		}
		
		double mean,max,sum,min = 0.0;
		
		//isos used in the fukui functions and MOs
		double iso1 = 0.0005;
		double iso2 = 0.005;
		double iso3 = -0.0005;
		double iso4 = -0.005;
		
		//iso for LH LCP
		lrdVol.lrds[9].get_cube_stats(mean,min,max);
		double iso5 = mean*2;
		double iso6 = iso5*30;
		
		//iso for LH Vee
		lrdVol.lrds[10].get_cube_stats(mean,min,max);
		double iso7 = mean;
		double iso8 = iso7*10;
		
		//iso for potential Fukui
		lrdVol.lrds[11].get_cube_stats(mean,min,max);
		double iso9 = mean;
		double iso10 = iso9*10;
		
		
		//iso for LH int 
		lrdVol.lrds[14].get_cube_stats(mean,min,max);
		double iso11 = mean;
		double iso12 = iso11*10;
		
		
		lrdVol.lrds[15].get_cube_stats(mean,min,max);
		double iso13 = min*0.1;
		double iso14 = min*0.3;
		double iso15 = max*0.1;
		double iso16 = max*0.3;
		
		//iso for right potential Fukui
		lrdVol.lrds[12].get_cube_stats(mean,min,max);
		double iso17 = mean*4;
		double iso18 = iso17*10;
		
		//iso for zero potential Fukui
		lrdVol.lrds[13].get_cube_stats(mean,min,max);
		double iso19 = mean*2;
		double iso20 = iso19*15;
		
		//iso for 
		lrdVol.lrds[16].get_cube_stats(mean,min,max);
		double iso21 = mean;
		double iso22 = iso21*10;
		
		//iso for 
		lrdVol.lrds[17].get_cube_stats(mean,min,max);
		double iso23 = 0.0002;
		double iso24 = iso23*5;
		
		//iso for multiphilicity
		lrdVol.lrds[18].get_cube_stats(mean,min,max);
		double iso25 = -max*0.3;
		double iso26 = -max*0.6;
		double iso27 = max*0.3;
		double iso28 = max*0.6;
		
		//iso for MEP
		lrdVol.lrds[19].get_cube_stats(mean,min,max);
		double iso29 = 0.05;
		double iso30 = 0.5;
		double iso31 = -0.5;
		double iso32 = -0.05;
		
		//iso for 
		lrdVol.lrds[20].get_cube_stats(mean,min,max);
		double iso33 = max*0.01;
		double iso34 = max*0.2;
		
		
		script_file	<< "load "	<< pdb_name									<< " \n"
					<< "load "	<< typestr << lrdVol.lrds[0].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[0].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[1].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[1].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[5].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[6].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[7].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[8].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[8].name			<< "_ph2.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[9].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[10].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[11].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[14].name			<< ".cube\n"
					<< "load "	<< typestr << lrdVol.lrds[15].name			<< "_ph1.cube\n"
					<< "load "	<< typestr << lrdVol.lrds[15].name			<< "_ph2.cube\n"
					<< "volume " << typestr2 << lrdVol.lrds[0].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[0].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[0].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[0].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[1].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[1].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[1].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[1].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[5].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[5].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[6].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[6].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[7].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[7].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[8].name			<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[8].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[8].name			<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[8].name	<<	"_ph2"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[9].name			<< "_volume, "		<< typestr2 << lrdVol.lrds[9].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[10].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[10].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[11].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[11].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[14].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[14].name				<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[15].name		<< "_ph1_volume, "	<< typestr2 << lrdVol.lrds[15].name	<<	"_ph1"	<< " \n" 
					<< "volume " << typestr2 << lrdVol.lrds[15].name		<< "_ph2_volume, "	<< typestr2 << lrdVol.lrds[15].name	<<	"_ph2"	<< " \n" 
					<< "volume_color " << typestr2 << lrdVol.lrds[0].name	<< "_ph1_volume, "	<< iso1		<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[0].name	<< "_ph2_volume, "	<< iso4		<< " red 0.04 "		<< iso3		<< " orange 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[1].name	<< "_ph1_volume, "	<< iso1 	<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[1].name	<< "_ph2_volume, "	<< iso4		<< " red 0.04 "		<< iso3		<< " orange 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[5].name	<< "_volume, "		<< iso1		<< " cyan 0.02 "	<< iso2		<< " blue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[6].name	<< "_volume, "		<< iso1		<< " pink 0.02 "	<< iso2		<< " red  0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[7].name	<< "_volume, "		<< iso1		<< " yellow 0.02 "	<< iso2		<< " green 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[8].name	<< "_ph1_volume, "	<< iso1		<< " pink 0.03 "	<< iso2		<< " red 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[8].name	<< "_ph2_volume, "	<< iso4		<< " blue 0.05 "	<< iso3		<< " cyan 0.02 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[9].name	<< "_volume, "		<< iso5		<< " limon 0.02 "	<< iso6		<< " purpleblue 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[10].name	<< "_volume, "		<< iso7		<< " brightorange 0.001 "	<< iso8	<< " purple 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[11].name	<< "_volume, "		<< iso9		<< " aquamarine 0.002 "		<< iso10 << " purple 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[14].name	<< "_volume, "		<< iso11	<< " salmon 0.02 "			<< iso12	<< " slate 0.05 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[15].name	<< "_ph1_volume, "	<< iso14	<< " skyblue  0.04 "		<< iso13	<<" greencyan 0.01 \n"
					<< "volume_color " << typestr2 << lrdVol.lrds[15].name	<< "_ph2_volume, "	<< iso15	<< " yelloworange 0.02 "	<< iso16	<< " warmpink 0.05 \n";
		
		if ( extra_RD ){
			script_file	<< "load "	<< typestr << lrdVol.lrds[12].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[13].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[16].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[17].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[18].name			<< "_ph1.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[18].name			<< "_ph2.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[19].name			<< "_ph1.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[19].name			<< "_ph2.cube\n"
						<< "load "	<< typestr << lrdVol.lrds[20].name			<< ".cube\n"
						<< "load "	<< typestr << lrdVol.lrds[21].name			<< ".cube\n"
						<< "volume " << typestr2 << lrdVol.lrds[12].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[12].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[13].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[13].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[16].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[16].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[17].name		<< "_volume, "		<< typestr2 << lrdVol.lrds[17].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[18].name		<< "_ph1_volume, "	<< typestr2	<< lrdVol.lrds[18].name	<< "_ph1"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[18].name		<< "_ph2_volume, "	<< typestr2	<< lrdVol.lrds[18].name	<< "_ph2"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[19].name		<< "_ph1_volume, "	<< typestr2	<< lrdVol.lrds[19].name	<< "_ph1"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[19].name		<< "_ph2_volume, "	<< typestr2	<< lrdVol.lrds[19].name	<< "_ph2"	<<	" \n" 
						<< "volume " << typestr2 << lrdVol.lrds[20].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[20].name	<< " \n" 
						<< "volume " << typestr2 << lrdVol.lrds[21].name		<< "_volume, "		<< typestr2	<< lrdVol.lrds[21].name	<< " \n" 
						<< "volume_color " << typestr2 << lrdVol.lrds[12].name	<< "_volume, "		<< iso17	<< " cyan 0.0 "				<< iso18	<< " blue 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[13].name	<< "_volume, "		<< iso19	<< " yellow 0.0 "			<< iso20	<< " orange 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[16].name	<< "_volume, "		<< iso21 	<< " cyan 0.0 "				<< iso22	<< " blue 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[17].name	<< "_volume, "		<< iso23	<< " yellow 0.02 "			<< iso24	<< " orange 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[18].name	<< "_ph1_volume, "	<< iso26	<< " blue 0.04 "			<< iso25	<< " cyan 0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[18].name	<< "_ph2_volume, "	<< iso27	<< " red 0.05 "				<< iso28	<< " pink  0.02 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[19].name	<< "_ph1_volume, "	<< iso29	<< " pink 0.0 "				<< iso30	<< " red  0.04 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[19].name	<< "_ph2_volume, "	<< iso31	<< " blue 0.04 "			<< iso32 << " aquamarine  0.0 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[20].name	<< "_volume, "		<< iso33	<< " greencyan 0.02 "		<< iso34	<< " green  0.05 \n"
						<< "volume_color " << typestr2 << lrdVol.lrds[21].name	<< "_volume, "		<< iso7		<< " brightorange 0.001 "	<< iso8 	<< " purple 0.05 \n";
		}
	}
	
	
}
/**************************************************************************/
void scripts::write_pymol_pdb(){
	
	
	//std::ofstream script_f;
	string fname = get_file_name( file_name.c_str() );
	fname += ".pym";
	//script_f.open( fname.c_str() );

	string fname2 = fname.substr( 0,fname.size()-4 );

	string load_pdb_basic = "load "+ fname2 + "_PDB_RD/" + fname2;
	
	script_file	<< load_pdb_basic  << "_nucleophilicity.pdb\n"	//0
				<< load_pdb_basic  << "_electrophilicity.pdb\n"	//1
				<< load_pdb_basic  << "_radicality.pdb\n"		//2
				<< load_pdb_basic  << "_netphilicity.pdb\n"		//3
				<< load_pdb_basic  << "_hardness_Vee.pdb\n"		//4
				<< load_pdb_basic  << "_hardness_lcp.pdb\n"		//5
				<< load_pdb_basic  << "_fukui_pot_left.pdb\n"	//6
				<< load_pdb_basic  << "_fukui_pot_right.pdb\n"	//7
				<< load_pdb_basic  << "_fukui_pot_zero.pdb\n"	//8
				<< load_pdb_basic  << "_softness_dual.pdb\n"	//9
				<< load_pdb_basic  << "_hyper_softness.pdb\n"	//10
				<< load_pdb_basic  << "_fukushima.pdb\n"		//11
				<< load_pdb_basic  << "_mep.pdb\n"				//12
				<< load_pdb_basic  << "_hardness_TFD.pdb\n"		//13
				<< load_pdb_basic  << "_softness_avg.pdb\n"		//14
				<< load_pdb_basic  << "_hardness_int.pdb\n"		//15
				<< load_pdb_basic  << "_multiphilicity.pdb\n"	//16
				<< load_pdb_basic  << "_charge.pdb\n"			//17
				<< load_pdb_basic  << "_electron_density.pdb\n"	//18
				<< "spectrum b, blue_white_red, minimum=-0.3, maximum=0.3\n"
				<< "spectrum b, white_yellow_orange_red_black, minimum=0.1, maximum=0.5\n"
				<< "spectrum b, white_cyan_blue, minimum=0, maximum=0.1\n"
				<< "spectrum b, white_pink_red, minimum=0, maximum=0.1\n";
}
/***************************************************************************/
void scripts::write_r_heatmap(vector< vector<double> > rd_numerical	, 
							vector<string> rds						,
							vector<string> residues)				{
								
	string sname = file_name;
	sname+="residuesRD_4_R";
	std::ofstream data_txt(sname.c_str() );
	
	data_txt << "residue ";
	for(unsigned i=0;i<rds.size();i++){
		data_txt << rds[i] << " ";
	}
	data_txt << "\n";
	
	for(unsigned i=0;i<residues.size();i++){
		data_txt << residues[i] << " ";
		for(unsigned j=0;j<rds.size();j++){
			data_txt << rd_numerical[i][j] << " "; 
		}
		data_txt << "\n";
	}
	
	script_file << "library(pheatmap)\n"
				<< "data1 <-read.table('" << sname << "',header=T)\n"
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
				<< sname << "heatmap.png'";
}
/*********************************************************************/
void scripts::write_r_residuos_barplot(){
	script_file << "require(ggpubr)\n"
				<< "data_stat <-read.table('residues_data_stat',header=T) \n"
				<< "avg <-subset(data_stat,type=='AVG')\n"
				<< "sd <-subset(data_stat,type=='SD')\n"
				<< "#-------------------------------------\n" 
				<< "p1 <-ggplot(avg, aes(x=res, y=EAS)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='blue', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Nucleophilicity')+\n"
				<< "xlab('Residues')+\n"				
				<< "geom_errorbar(aes(ymin=NAS-sd$EAS,ymax=EAS+sd$EAS), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p2 <-ggplot(avg, aes(x=res, y=NAS)) +\n" 
				<< "geom_bar(stat='identity', color='black',fill='red', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Electrophilicity')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=NAS-sd$NAS,ymax=NAS+sd$NAS), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p3 <-ggplot(avg, aes(x=res, y=Netphilicity)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='lightgreen', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Netphilicity')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=Netphilicity-sd$Netphilicity,ymax=Netphilicity+sd$Netphilicity), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p4 <-ggplot(avg, aes(x=res, y=Hardness_A)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='orange', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Hardness (LCP)')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=Hardness_A-sd$Hardness_A,ymax=Hardness_A+sd$Hardness_A), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p5 <-ggplot(avg, aes(x=res, y=Hardness_B)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='orange', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Hardness (MEP EE)')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=Hardness_B-sd$Hardness_B,ymax=Hardness_B+sd$Hardness_B), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p6 <-ggplot(avg, aes(x=res, y=Hardness_C)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='orange', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Hardness (Fukui Potential)')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=Hardness_C-sd$Hardness_C,ymax=Hardness_C+sd$Hardness_C), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "p7 <-ggplot(avg, aes(x=res, y=Electron_Density)) +\n"
				<< "geom_bar(stat='identity', color='black',fill='orange', position=position_dodge() ) +\n"
				<< "theme_minimal()+\n"
				<< "ylab('Electron Density')+\n"
				<< "xlab('Residues')+\n"
				<< "geom_errorbar(aes(ymin=Electron_Density-sd$Electron_Density,ymax=Electron_Density+sd$Electron_Density), width=.2, position=position_dodge(.9))\n"
				<< "#-------------------------------------\n" 
				<< "png('EAS_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p1\n"
				<< "dev.off()\n"				
				<< "png('NAS_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p2\n"
				<< "dev.off()\n"				
				<< "png('NET_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p3\n"
				<< "dev.off()\n"
				<< "png('Hardness_A_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p4\n"
				<< "dev.off()\n"
				<< "png('Hardness_B_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p5\n"
				<< "dev.off()\n"
				<< "png('Hardness_C_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p6\n"
				<< "dev.off()\n"
				<< "png('Electron_dens_res.png',units='in',res=1000,width=4.5,height=4)\n"
				<< "p7\n"
				<< "dev.off()\n"
				<< "#=====================================================\n" 				
				<< "dat <-read.table('residues_data_frames',header=T)\n"
				<< "#-------------------------------------\n"	
				<< "pm1<-ggplot(dat, aes(x = frame, y =NAS ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Electrophiliicty')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm2<-ggplot(dat, aes(x = frame, y =EAS ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Nucleophilicity')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm3<-ggplot(dat, aes(x = frame, y =Netphilicity ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"				
				<< "ylab('Netphilicity')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm4<-ggplot(dat, aes(x = frame, y =Hardness_A ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Hardness_A')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm5<-ggplot(dat, aes(x = frame, y =Hardness_B ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Hardness_B')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm6<-ggplot(dat, aes(x = frame, y =Hardness_C ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Hardness_C')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm7<-ggplot(dat, aes(x = frame, y =Softness ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Softness')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "pm8<-ggplot(dat, aes(x = frame, y =Electron_Density ))+\n"
				<< "theme_minimal()+\n"
				<< "geom_point(aes(color=res))+\n"
				<< "geom_smooth(aes(color=res),method='loess')+\n"
				<< "ylab('Electron Density')+\n"
				<< "xlab('Frame')\n"
				<< "#-------------------------------------\n"
				<< "wdht = 5\n"
				<< "png('eas_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm1\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('nas_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm2\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"				
				<< "png('net_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm3\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('hard_A_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm4\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('hard_B_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm5\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('hard_C_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm6\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"				
				<< "png('Softness_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm7\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n"
				<< "png('ED_mov_avg.png',units='in',res=1000,width=wdht,height=4)\n"
				<< "pm8\n"
				<< "dev.off()\n"
				<< "#-------------------------------------\n";			
}
/*********************************************************************/
void scripts::write_r_reaction_analysis(traj_rd& path_rd			,
										vector<string>& pair_labels	,
										ReactionAnalysis& r_info	,
										string& nameb				){
											
	
	vector<string> gl_names		= {"HOF","ECP","Hardness","Softness","Electrophilicity","Energy"};									
	vector<string> gl_legends	= {"HOF \\n(kCal/mol)",
									"ECP (eV)",
									"Hardness \\n(eV)",
									"Softness \\n(eV)",
									"Electrophilicity \\n(eV)",
									"Energy \\n(eV)"};									
	
	string delim = "#====================================================\n";
	string delim2 = "#---------------------------------------------------\n";
	
	script_file		<< delim
					<< "require(ggpubr)\n"
					<< "jet.colors <-colorRampPalette(c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow','#FF7F00','red','#7F0000'))\n\n"
					<< delim
					<< "atom_lrd='" << nameb <<  "'\n"
					<< "rc1_name='" << r_info.RCs[0].rc_label << "'\n"
					<< "df1 <-read.table(atom_lrd,header=T)\n";
					if ( r_info.nrcs == 2 ){
						script_file << "rc2_name='" << r_info.RCs[1].rc_label << "'\n";
					}
	script_file << delim;

	string pr_obj;
	
	unsigned int a = 0;
	unsigned int c = 0;	
	
	if ( r_info.ndim == 1 ){
		for( unsigned i=0; i<gl_names.size(); i++){
			script_file << "grc1_" << to_string(i) << " <-ggplot(df1,aes(x=RC1,y="
						<< gl_names[i] << ") )+\n "
						<< "\tgeom_point()+ \n"
						<< "\tgeom_line()+ \n"
						<< "\ttheme_minimal()+ \n"
						<< "\tylab('" << gl_legends[i] << "') + \n"
						<< "\txlab(rc1_name)\n\n"
						<< delim2;
		}
	
		script_file << "png('global_rc1.png',width=6.5,height=5,units='in',res=1000)\n";
		script_file << "ggarrange(grc1_0,grc1_1,grc1_2,grc1_3,grc1_4,grc1_5,ncol=3,nrow=2)\n";
		script_file << "dev.off()\n";
		script_file << delim;
		
		for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
			pr_obj = "la_";
			script_file<< pr_obj
						<< to_string(i) << " <-ggplot(df1,aes(x=RC1,y="
						<< path_rd.rds_labels[i] << ") )+\n "
						<< "\tgeom_point()+ \n"
						<< "\tgeom_line()+ \n"
						<< "\ttheme_minimal()+ \n"
						<< "\tylab('" << path_rd.atoms_labels[i] << "') + \n"
						<< "\txlab(rc1_name)  \n"
						<< delim2;
			c++;
			if( c % 9 == 0 ){
				script_file << "png('" <<pr_obj << to_string(a++) << ".png',width=8,height=7,units='in',res=1000)\n";
				script_file << "ggarrange(" 
							<< pr_obj << to_string(c-9) << ","
							<< pr_obj << to_string(c-8) << ","
							<< pr_obj << to_string(c-7) << ","
							<< pr_obj << to_string(c-6) << ","
							<< pr_obj << to_string(c-5) << ","
							<< pr_obj << to_string(c-4) << ","
							<< pr_obj << to_string(c-3) << ","
							<< pr_obj << to_string(c-2) << ","
							<< pr_obj << to_string(c-1) << ","
							<< "ncol=3,nrow=3)\n";
				script_file << "dev.off()\n";
				script_file << delim;
			}		
		}			
	
		c = 0;
		for( unsigned i=0; i<pair_labels.size(); i++){
			pr_obj = "pr_";
			pr_obj += to_string(c);
			pr_obj += "_c1_";
			int j = i+1;
			script_file<< pr_obj
						<< to_string(j) << "<-ggplot(df1,aes(x=RC1,y="
						<< pair_labels[i] << ") )+\n "
						<< "\tgeom_point()+ \n"
						<< "\tgeom_line()+ \n"
						<< "\ttheme_minimal()+ \n"
						<< "\tylab('" << pair_labels[i] << "') + \n"
						<< "\txlab(rc1_name)\n"
						<< delim2;
			if ( j % 6 == 0 ) { 
				c++;		
				script_file << "png('" << pr_obj << to_string(j) << ".png',width=7,height=6,units='in',res=1000)\n";
				script_file << "ggarrange(" << pr_obj << to_string(j-5) << ","
							<< pr_obj << to_string(j-4) << ","
							<< pr_obj << to_string(j-3) << ","
							<< pr_obj << to_string(j-2) << ","
							<< pr_obj << to_string(j-1) << ","
							<< pr_obj << to_string(j)   << ","
							<< "ncol=3,nrow=2)\n";					 
				script_file << "dev.off()\n";		
				script_file << delim;
			}
		}

		if ( r_info.nrcs == 2 ){
			c = 0;
			a = 0;
			for( unsigned i=0;i<gl_names.size(); i++){
				script_file << " grc2_" << to_string(i) << " <-ggplot(df1,aes(x=RC2,y="
							<< gl_names[i] << ") )+\n "
							<< "\tgeom_point()+ \n"
							<< "\tgeom_line()+ \n"
							<< "\ttheme_minimal()+ \n"
							<< "\tylab('" << gl_legends[i] << "') + \n"
							<< "\txlab(rc2_name) \n"
							<< delim2;
			}
			script_file << "png('global_rc2.png',width=6.5,height=5,units='in',res=1000)\n"
						<< "ggarrange(grc2_0,grc2_1,grc2_2,grc2_3,grc2_4,grc2_5,ncol=3,nrow=2)\n"
						<< "dev.off()\n"
						<< delim;
		
			for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
				pr_obj = "la_";
				script_file << pr_obj
							<< to_string(i) << " <-ggplot(df1,aes(x=RC2,y="
							<< path_rd.rds_labels[i] << ") )+\n "
							<< "\tgeom_point()+ \n"
							<< "\tgeom_line()+ \n"
							<< "\ttheme_minimal()+ \n"
							<< "\tylab('" << path_rd.atoms_labels[i] << "') + \n"
							<< "\txlab(rc2_name)  \n"
							<< delim2;
				c++;
				if( c % 9 == 0 ){
					script_file	<< "png('" <<pr_obj << "rc2_" << to_string(a++) <<".png',width=8,height=7,units='in',res=1000)\n";
					script_file	<< "ggarrange(" 
								<< pr_obj << to_string(c-9) << ","
								<< pr_obj << to_string(c-8) << ","
								<< pr_obj << to_string(c-7) << ","
								<< pr_obj << to_string(c-6) << ","
								<< pr_obj << to_string(c-5) << ","
								<< pr_obj << to_string(c-4) << ","
								<< pr_obj << to_string(c-3) << ","
								<< pr_obj << to_string(c-2) << ","
								<< pr_obj << to_string(c-1) << ","					 
								<< "ncol=3,nrow=3)\n";					 
					script_file<< "dev.off()\n";		
					script_file<< delim;
				}		
			}	
			c = 0;
			for( unsigned i=0; i<pair_labels.size(); i++){
				pr_obj = "pr_";
				pr_obj += to_string(c);
				pr_obj += "_c2_";
				int j = i+1;
				script_file<< pr_obj
							<< to_string(j) << "<-ggplot(df1,aes(x=RC2,y="
							<< pair_labels[i] << ") )+\n "
							<< "\tgeom_point()+ \n"
							<< "\tgeom_line()+ \n"
							<< "\ttheme_minimal()+ \n"
							<< "\tylab('" << pair_labels[i] << "') + \n"
							<< "\txlab(rc2_name)\n"
							<< delim2;
				if ( j % 6 == 0 ) { 
					c++;		
					script_file<< "png('" << pr_obj << "rc2_" << to_string(j) << ".png',width=7,height=6,units='in',res=1000)\n";
					script_file<< "ggarrange(" << pr_obj << to_string(j-5) << ","
								<< pr_obj << to_string(j-4) << ","
								<< pr_obj << to_string(j-3) << ","
								<< pr_obj << to_string(j-2) << ","
								<< pr_obj << to_string(j-1) << ","
								<< pr_obj << to_string(j)
								<< ",ncol=3,nrow=2)\n";
					script_file << "dev.off()\n";
					script_file << delim;
				}
			}		
		}	
	}else if ( r_info.ndim == 2 ){
		for( unsigned i=0; i<gl_names.size(); i++){
			script_file << "grc1_" << to_string(i) << " <-ggplot(df1,aes(x=rc1,y=rc2,z="
						<< gl_names[i] << ") )+\n "
						<< "\tstat_contour(geom='polygon', aes(fill = ..level..))+ \n" 
						<< "\tgeom_tile(aes(fill =" << gl_names[i] << " ))+ \n"
						<< "\tstat_contour(bins = 6,color='black')+ \n"
						<< "\tscale_fill_gradientn(colours=jet.colors(70))+ \n"
						<< "\tylab(rc1_name) + \n"
						<< "\txlab(rc2_name) + \n"
						<< "\tguides(fill = guide_colorbar('" << gl_legends[i] <<"'))\n"
						<< delim2;
		}
	
		script_file << "png('global_2d.png',width=8,height=8,units='in',res=1000)\n";
		script_file << "ggarrange(grc1_0,grc1_1,grc1_2,grc1_3,grc1_4,grc1_5,ncol=2,nrow=3)\n";
		script_file << "dev.off()\n";
		script_file << delim;
		
		for ( unsigned i=0; i<path_rd.rds_labels.size(); i++){
				pr_obj = "la_";
				script_file << pr_obj
							<< to_string(i) << " <-ggplot(df1,aes(x=rc1,y=rc2,z="
							<< path_rd.rds_labels[i] << ") )+\n "
							<< "\tstat_contour(geom='polygon', aes(fill = ..level..))+ \n" 
							<< "\tstat_contour(bins = 6,color='black')+ \n"
							<< "\tscale_fill_gradientn(colours=jet.colors(7))+\n"
							<< "\tgeom_tile(aes(fill =" << path_rd.rds_labels[i] << "))+\n"
							<< "\tylab(rc2_name) + \n"
							<< "\txlab(rc1_name)  +\n"
							<< "\tguides(fill = guide_colorbar(title ='" << path_rd.atoms_labels[i] << "'))\n"
							<< delim2;
				c++;
			if( c % 9 == 0 ){
				script_file	<< "png('" <<pr_obj << "rc2_" << to_string(a++) <<".png',width=11,height=7,units='in',res=1000)\n";
				script_file	<< "ggarrange(" 
							<< pr_obj << to_string(c-9) << ","
							<< pr_obj << to_string(c-8) << ","
							<< pr_obj << to_string(c-7) << ","
							<< pr_obj << to_string(c-6) << ","
							<< pr_obj << to_string(c-5) << ","
							<< pr_obj << to_string(c-4) << ","
							<< pr_obj << to_string(c-3) << ","
							<< pr_obj << to_string(c-2) << ","
							<< pr_obj << to_string(c-1) << ","
							<< "ncol=3,nrow=3)\n";
					script_file<< "dev.off()\n";
					script_file<< delim;
			}		
		}
		for( unsigned i=0; i<pair_labels.size(); i++){
				pr_obj = "pr_";
				pr_obj += to_string(c);
				pr_obj += "_c2_";
				int j = i+1;
				script_file<< pr_obj
							<< to_string(j) << "<-ggplot(df1,aes(x=rc1,y=rc2,z="
							<< pair_labels[i] << ") )+\n "
							<< "\tstat_contour(geom='polygon', aes(fill = ..level..))+ \n" 
							<< "\tstat_contour(bins = 6,color='black')+ \n"
							<< "\tscale_fill_gradientn(colours=jet.colors(7))+\n"
							<< "\tgeom_tile(aes(fill =" << pair_labels[i] << "))+\n"
							<< "\tylab(rc1_name) + \n"
							<< "\txlab(rc2_name) + \n"
							<< "\tguides(fill = guide_colorbar(title ='" << path_rd.atoms_labels[i] << "'))\n"
							<< delim2;
				if ( j % 6 == 0 ) { 
					c++;		
					script_file<< "png('" << pr_obj << "rc2_" << to_string(j) << ".png',width=11,height=7,units='in',res=1000)\n";
					script_file<< "ggarrange(" << pr_obj << to_string(j-5) << ","
								<< pr_obj << to_string(j-4) << ","
								<< pr_obj << to_string(j-3) << ","
								<< pr_obj << to_string(j-2) << ","
								<< pr_obj << to_string(j-1) << ","
								<< pr_obj << to_string(j)
								<< ",ncol=3,nrow=2)\n";
					script_file << "dev.off()\n";
					script_file << delim;
			}
		}
	}
}

//================================================================================
//END OF FILE
//================================================================================

