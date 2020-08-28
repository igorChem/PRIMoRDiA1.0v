//include c++ headers 
#include <iostream>
#include <string>
#include <fstream>
#include <memory>
//include primordia headers
#include "../include/common.h"
#include "../include/Imolecule.h"
#include "../include/Icube.h"
#include "../include/global_rd.h"
#include "../include/local_rd.h"
#include "../include/local_rd_cnd.h"
#include "../include/gridgen.h"
#include "../include/QMparser.h"
#include "../include/primordia.h"
#include "../include/Iprotein.h"
#include "../include/scripts.h"
#include "../include/residue_lrd.h"

using std::unique_ptr;
using std::string;
using std::cout;
using std::endl;
using std::move;
/*************************************************************************************/
primordia::primordia()			:
	name("nonamed")				,
	grd( new global_rd() )		,
	lrdVol( new local_rd() )	,
	lrdCnd( new local_rd_cnd() ){
}
/*************************************************************************************/
primordia::primordia(const primordia& pr_rhs)	:
	name(pr_rhs.name)							,
	grd( new global_rd (*pr_rhs.grd) )			,
	lrdVol( new local_rd(*pr_rhs.lrdVol) )		,
	lrdCnd( new local_rd_cnd(*pr_rhs.lrdCnd) ){
}
/***************************************************************************************/
primordia& primordia::operator=(const primordia& pr_rhs){
	if ( this != &pr_rhs ) {
		name	= pr_rhs.name;
		*grd	= *pr_rhs.grd;
		*lrdVol	= *pr_rhs.lrdVol;
		*lrdCnd	= *pr_rhs.lrdCnd;
	}
	return *this;
}
/***************************************************************************************/
primordia::primordia(primordia&& pr_rhs) noexcept:
	name( move(pr_rhs.name) )					,
	grd( move(pr_rhs.grd) )						,
	lrdVol( move(pr_rhs.lrdVol) )				,
	lrdCnd( move(pr_rhs.lrdCnd) )				{
}
/***************************************************************************************/
primordia& primordia::operator=(primordia&& pr_rhs) noexcept{
	if ( this != &pr_rhs ) {
		name	= move(pr_rhs.name);
		grd		= move(pr_rhs.grd);
		lrdVol	= move(pr_rhs.lrdVol);
		lrdCnd	= move(pr_rhs.lrdCnd);
	}
	return *this;
}
/***************************************************************************************/
primordia operator-(const primordia& pr_lhs, const primordia& pr_rhs){
	primordia Result(pr_lhs);
	*Result.grd		= *pr_lhs.grd		- *pr_rhs.grd;
	*Result.lrdVol	= *pr_lhs.lrdVol	- *pr_rhs.lrdVol;
	*Result.lrdCnd	= *pr_lhs.lrdCnd	- *pr_rhs.lrdCnd;
	return Result;
}
/*************************************************************************************/
void primordia::init_FOA(const char* file_neutro,
							   int grdN			,
							   string loc_hard	,
							   bool mep			, 
							   string Program)	{

	name = remove_extension(file_neutro);
	//-------------------------------------------------------------------
	// log messages
	m_log->input_message("molecule name "+name);
	m_log->input_message("Parameters for calculating the Reactivity descriptors: ");
	m_log->input_message("Approximation: Frozen Orbital");
	m_log->input_message("Grid size: "+int_to_string(grdN) );
	m_log->input_message("Program QM output: "+Program);
	m_log->input_message("local hardness method: "+loc_hard);
	//-------------------------------------------------------------------
	
	unique_ptr<QMparser> qmfile( new QMparser(file_neutro,Program) ); 
	Imolecule molecule( qmfile->get_molecule() ); 
	qmfile.reset(nullptr);
	if ( molecule.name == "empty"){
		m_log->input_message("Skipping reactivity descriptors calculations for:");
		m_log->input_message(name);
	}else{
		if ( dos ){
			scripts dos(molecule.name.c_str());
			dos.write_r_dos(molecule.orb_energies);
		}
	
		grd.reset( new global_rd(molecule) );
		grd->calculate_rd();
		grd->write_rd();
		lrdCnd.reset( new local_rd_cnd( move(molecule) ) );
		lrdCnd->calculate_Fukui();
		lrdCnd->calculate_RD(*grd);
		lrdCnd->calculate_Hardness(*grd);
		lrdCnd->write_LRD();

		if ( grdN  > 0 ){
			unique_ptr<gridgen> grid1( new gridgen(grdN,move(lrdCnd->molecule) ) );
			Icube homo_cub	= grid1->calc_HOMO_density();
			Icube lumo_cub	= grid1->calc_LUMO_density();
			if ( loc_hard == "mepEE" || loc_hard == "LCP" ) {
				if ( Program == "orca" ) grid1->calculate_density_orca();
				else { grid1->calculate_density(); }
				Icube dens =  grid1->density;
				lrdVol.reset( new local_rd( dens,homo_cub,lumo_cub ) );
			}else { lrdVol.reset( new local_rd( homo_cub,lumo_cub )  ); }
		
			lrdVol->calculate_fukui();
			lrdVol->calculate_RD(*grd);
			if ( loc_hard !="not" ){ lrdVol->calculate_hardness(*grd,loc_hard); }
			lrdVol->write_LRD();
			if ( mep ){
				grid1->calculate_mep_from_charges();
				grid1->density.write_cube(grid1->name+"_mep.cube");
				grid1->density.write_cube(grid1->name+"_mep2.cube");
				m_log->input_message("MEP grid calculations required ");
			}
		}
		if ( pymol_script ) { 
			this->write_pymol("none",true); 
		}
	}
	
}
/*************************************************************************************/
void primordia::init_FD(const char* file_neutro			,
							  const char* file_cation	,
							  const char* file_anion	, 
							  const int grdN			, 
							  int charge				,
							  bool mep					,
							  string loc_hard			, 
							  string Program)			{

	name	= remove_extension(file_neutro);
	
	//-------------------------------------------------------------------
	// log messages
	m_log->input_message("molecule name "+name);
	m_log->input_message("Parameters for calculating the Reactivity descriptors:");
	m_log->input_message("Approximation: Finite Differences");
	m_log->input_message("Grid size : "+int_to_string(grdN) );
	m_log->input_message("Program QM output: "+Program);
	m_log->input_message("local hardness method "+loc_hard);
	//----------------------------------------------------------------------
	
	unique_ptr<QMparser> qmfile1 ( new QMparser(file_neutro,Program) );
	Imolecule molecule_a ( qmfile1->get_molecule() ) ;
	qmfile1.reset(nullptr);
	unique_ptr<QMparser> qmfile2 ( new QMparser(file_cation,Program) );
	Imolecule molecule_b (  qmfile2->get_molecule() ) ;
	qmfile2.reset(nullptr);
	unique_ptr<QMparser> qmfile3 ( new QMparser(file_anion,Program) );
	Imolecule molecule_c ( qmfile3->get_molecule() ) ;
	qmfile3.reset(nullptr);
	if ( molecule_a.name == "empty"){
		m_log->input_message("Skipping reactivity descriptors calculations for:");
		m_log->input_message(file_neutro);
	}else if (molecule_b.name == "empty") {
		m_log->input_message("Skipping reactivity descriptors calculations for:");
		m_log->input_message(file_cation);
	}else if (molecule_c.name =="empty"){
		m_log->input_message("Skipping reactivity descriptors calculations for:");
		m_log->input_message(file_anion);
	}else{
		grd.reset ( new global_rd( molecule_a,molecule_b,molecule_c) );
		grd->calculate_rd();
		grd->write_rd();
	
		lrdCnd.reset( new local_rd_cnd(molecule_a, molecule_b, molecule_c) );
	
		lrdCnd->calculate_RD(*grd);
		lrdCnd->calculate_Hardness(*grd);
		lrdCnd->write_LRD();
	
		if ( grdN > 0 ){
			unique_ptr<gridgen> grid1 ( new gridgen( grdN,move(molecule_a) ) );
			if ( Program == "orca" ) grid1->calculate_density_orca();
			else {	grid1->calculate_density();	}
	 
			unique_ptr<gridgen> grid2 ( new gridgen( grdN,move(molecule_b) ) );
			if ( Program == "orca" ) grid2->calculate_density_orca();
			else { grid2->calculate_density();	}
	
			unique_ptr<gridgen> grid3 ( new gridgen( grdN, move(molecule_c) ) );
			if ( Program == "orca" ) grid3->calculate_density_orca();
			else { grid3->calculate_density();	}
		
			lrdVol.reset ( new local_rd(grid1->density,grid2->density,grid3->density,charge) );	
			lrdVol->calculate_fukui();
			lrdVol->calculate_RD(*grd);
			lrdVol->calculate_hardness(*grd,loc_hard);
			lrdVol->write_LRD();
			if ( mep ){
				grid1->calculate_mep_from_charges();
				grid1->density.write_cube(name+"_mep.cube");
				grid1->density.write_cube(name+"_mep2.cube");
				m_log->input_message("MEP grid calculations required");
			}
		}
		if ( pymol_script ) { this->write_pymol("none", true); }
	}
}
/*************************************************************************************/
void primordia::init_protein_RD(const char* file_name	,
								string locHardness		,
								int gridN				,
								int bandgap				,
								double* ref_atom		,
								int size				,
								const char* pdb			,
								bool mep				,
								string bt				,
								string Program			){
	band		= bandgap;
	int band2	= bandgap;
	name = remove_extension(file_name);
	
	//-------------------------------------------------------------------
	// log messages
	m_log->input_message("molecule name "+name);
	m_log->input_message("Parameters for calculating the band Reactivity descriptors:");
	m_log->input_message("Approximation: Frozen Orbital");
	m_log->input_message("Grid size: "+int_to_string(gridN) );
	m_log->input_message("Program QM output: "+Program);
	m_log->input_message("local hardness method "+locHardness);
	m_log->input_message("Band Reactivity Descriptors method "+bt);
	m_log->input_message("Band  molecular orbitals max size "+bandgap);
	m_log->input_message("Energy criteria "+double_to_string(energy_crit)+" (eV)" );
	//---------------------------------------------------------------------
	
	unique_ptr<QMparser> fileQM ( new QMparser(file_name,Program) );
	Imolecule molecule ( fileQM->get_molecule()  );
	fileQM.reset(nullptr);
	if ( molecule.name == "empty"){
		m_log->input_message("Skipping reactivity descriptors calculations for:");
		m_log->input_message(file_name);
	}else{
		name = remove_extension(file_name);
		Iprotein pdbfile(pdb);
		if ( dos ){
			scripts dos(molecule.name.c_str());
			dos.write_r_dos(molecule.orb_energies);
		}
	
		grd.reset( new global_rd( molecule) );
		lrdCnd.reset( new local_rd_cnd( move(molecule) )  );
		lrdCnd->mep		= mep;
		grd->calculate_rd();
	
		if ( band > 0 ) {
			if 	( bt == "EW" )	lrdCnd->calculate_Fukui_EW(band);
			else{ 				lrdCnd->calculate_Fukui_band(band); }
			lrdCnd->calculate_RD(*grd);
			lrdCnd->calculate_Hardness(*grd);
			lrdCnd->rd_protein(pdb);
			lrdCnd->write_rd_protein_pdb(pdb);
			lrdCnd->write_LRD();
		}else{
			lrdCnd->calculate_Fukui();
			lrdCnd->calculate_RD(*grd);
			lrdCnd->calculate_Hardness(*grd);
			lrdCnd->rd_protein(pdb);
			lrdCnd->write_rd_protein_pdb(pdb);
			lrdCnd->write_LRD();
		}
		if ( gridN > 0 ){
			unique_ptr<gridgen> grid ( new gridgen(gridN,move(lrdCnd->molecule) ) );
			if ( size > 0 ) { 
				//cout << ref_atom[0] << " " << ref_atom[1] << " " << ref_atom[2] << endl;
				grid->redefine_lim(ref_atom[0],ref_atom[1],ref_atom[2],size);
			}
			Icube EAS;
			Icube NAS;
			if ( band  > 0 ) {
				if ( bt == "BD" ){
					EAS  = grid->calc_band_EAS(bandgap);
					NAS  = grid->calc_band_NAS(bandgap);
				}else if ( bt == "EW" ){
					EAS = grid->calc_EBLC_EAS();
					NAS = grid->calc_EBLC_NAS();
				}
				//EAS.print();
				if ( locHardness == "LCP" || locHardness == "mepEE" ){
					grid->calculate_density();
					lrdVol.reset( new local_rd(grid->density,EAS,NAS) );
				}else{ 
					lrdVol.reset( new local_rd(EAS,NAS) );
				}
				lrdVol->name = name;
				lrdVol->calculate_fukui();
				lrdVol->calculate_RD(*grd);
				lrdVol->calculate_hardness(*grd,locHardness);
				lrdVol->write_LRD();
			}else{
				EAS = grid->calc_HOMO_density() ;
				NAS = grid->calc_LUMO_density() ;
				if ( locHardness == "LCP" || locHardness == "mepEE" ){ 
					grid->calculate_density();
					lrdVol.reset( new local_rd(grid->density,EAS,NAS) );
				}else{ 
					lrdVol.reset( new local_rd(EAS,NAS) );
				}
				lrdVol->name = name;
				lrdVol->calculate_fukui();
				lrdVol->calculate_RD(*grd);
				lrdVol->calculate_hardness(*grd,locHardness);
				lrdVol->write_LRD();
			}
			if ( mep ){
				grid->calculate_mep_from_charges();
				grid->density.write_cube(name+"_mep.cube");
				grid->density.write_cube(name+"_mep2.cube");
			}
		}else{	lrdCnd->molecule.clear();	}
		if ( pymol_script ) { this->write_pymol(pdbfile.name, false); }
	}
}
/************************************************************************************/
void primordia::write_pymol(const std::string& pdb, bool fixed){
		
	if ( lrdVol->EAS.voxelN > 0 ){
		std::string typestr;
		if ( lrdVol->finite_diff )	{ typestr = "_FD";}
		else{ typestr = "_FOA";	}
		string sname = name + typestr +".pymol";
		std::ofstream pymols(sname.c_str() );
		
		double mean,max,sd,min = 0.0;
		double iso1 = 0.0005;
		double iso2 = 0.005;
		double iso3 = 0.0005;
		double iso4 = 0.005;
		if  ( !fixed ) {	
			lrdVol->EAS.get_cube_stats(mean,max,sd,min);
			iso1 = mean;
			iso2 = iso1*20;
			lrdVol->NAS.get_cube_stats(mean,max,sd,min);
			iso3 = mean;
			iso4 = iso3*20;
		}
		lrdVol->Hardness.get_cube_stats(mean,max,sd,min);
		double iso5 = mean;
		double iso6 = iso5*20;

		pymols << "load "		<< pdb						<<  " \n"
				<< "bg_color white\n" 
				<< "load "   << lrdVol->EAS.name				<<  ".cube\n"
				<< "load "   << lrdVol->NAS.name				<<  ".cube\n"
				<< "load "   << lrdVol->RAS.name				<<  ".cube\n"
				<< "load "   << lrdVol->Dual.name				<<  ".cube\n"
				<< "load "   << lrdVol->Hardness.name			<<  ".cube\n"
				<< "bg_color white"								<< "\n"
				<< "volume " << lrdVol->EAS.name				<< "_vol, "	<< lrdVol->EAS.name    << " \n" 
				<< "volume " << lrdVol->NAS.name				<< "_vol, " << lrdVol->NAS.name     << " \n" 
				<< "volume " << lrdVol->RAS.name				<< "_vol, " << lrdVol->RAS.name << " \n" 
				<< "volume " << lrdVol->Dual.name				<< "_vol, " << lrdVol->Dual.name     << " \n" 
				<< "volume " << lrdVol->Dual.name				<< "2_vol, "<< lrdVol->Dual.name     << " \n" 
				<< "volume " << lrdVol->Hardness.name			<< "_vol, " << lrdVol->Hardness.name  << " \n" 
				<< "volume_color " <<	lrdVol->EAS.name		<< "_vol, " << iso1 << " cyan 0.02 "   << iso2 << " blue 0.07 \n"
				<< "volume_color " <<	lrdVol->NAS.name		<< "_vol, " << iso3 << " pink 0.02 "   << iso4 << " red  0.07 \n"
				<< "volume_color " <<	lrdVol->RAS.name		<< "_vol, " << iso1 << " yellow 0.02 " << iso2 << " green 0.07 \n"
				<< "volume_color " <<	lrdVol->Dual.name		<< "_vol, " << iso1 << " pink 0.03 "   << iso2 << " red 0.07 \n"
				<< "volume_color " <<	lrdVol->Dual.name		<< "2_vol, "<< (-iso4) << " blue 0.07 "    << (-iso3) << " cyan 0.02 \n"
				<< "volume_color " <<	lrdVol->Hardness.name	<< "_vol, " << iso5 << " orange 0.003 "<< iso6 << " purple 0.07 \n";
		pymols.close();	
	}
	
	if ( pdb != "none" ){
		string sname = name + "_EAS.pymol";
		std::ofstream py_eas(sname.c_str());
		sname = name  +"_NAS.pymol";
		std::ofstream py_nas(sname.c_str());
		sname = name  +"_RAS.pymol";
		std::ofstream py_ras(sname.c_str());
		sname = name  + "_Dual.pymol";
		std::ofstream py_dual(sname.c_str());
		sname = name  + "_Hard.pymol";
		std::ofstream py_hard(sname.c_str());
	
		double min = min_dvec(lrdCnd->EAS);
		double mean = mean_dvec(lrdCnd->EAS);
		double max = 0.0;
		py_eas  << "load "   << name    << "_EAS.pdb\n"
				<< "bg_color white\n"
				<< "spectrum b, white_blue, minimum= " << min << "," << "maximum=" << mean <<  "\n"
				<< "as cartoon\n";
		py_eas.close();
	
		min = min_dvec(lrdCnd->NAS);
		mean = mean_dvec(lrdCnd->NAS);
		py_nas  << "load "   << name    << "_NAS.pdb\n"
				<< "bg_color white\n"
				<< "spectrum b, white_red, minimum=" << min << "," << "maximum=" <<mean <<  "\n"
				<< "as cartoon\n";
		py_nas.close();
	
		min = min_dvec(lrdCnd->RAS);
		mean = mean_dvec(lrdCnd->RAS);
		py_ras  	<< "load "   << name    << "_RAS.pdb\n"
						<< "bg_color white\n"
						<< "spectrum b, white_green, minimum="<< min << "," << "maximum=" <<mean <<  "\n"
						<< "as cartoon\n";
		py_ras.close();
	
		min = min_dvec(lrdCnd->dual);
		max = max_dvec(lrdCnd->dual);
		min /=2;
		max /=2;
		py_dual 	<< "load "   << name    << "_dual.pdb\n"
						<< "bg_color white\n"
						<< "spectrum b, blue_white_red, minimum="<< min << "," << "maximum=" <<max <<  "\n"
						<< "as cartoon\n";
		py_dual.close();
	
		min = min_dvec(lrdCnd->hardness_A);
		max = max_dvec(lrdCnd->hardness_A);
		py_hard  << "load "   << name    << "_hardness.pdb\n"
						<< "bg_color white\n"
						<< "spectrum b, white_yellow_orange_red, minimum="<< min << "," << "maximum=" <<max <<  "\n"
						<< "as cartoon\n";
		py_hard.close();
	
	}
	
}
/***************************************************************************************/
primordia::~primordia(){};
//================================================================================
//END OF FILE
//================================================================================