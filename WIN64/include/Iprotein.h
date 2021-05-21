#ifndef PROTEIN
#define PROTEIN 
//------------------------
#include <vector>
#include <string>
//------------------------
#include "../include/common.h"
//------------------------
class Imolecule;
//=====================================================
/**
 * @class Iresidue
 * @author igor
 * @date 09/03/20
 * @file Iprotein.h
 * @brief 
 */
class Iresidue{
	public:
		std::string name;
		std::string type;
		std::vector<std::string> atom_type;
		unsigned int atom_s;
		unsigned int num;
		bool ligand;
		double molar_vol;
		std::vector<int> atoms_index;
		std::vector<int> side_chain_atoms;
		std::vector<int> back_bone_atoms;
		//--------------------------------------------------------------------
		Iresidue();
		Iresidue(const Iresidue& res);
		Iresidue& operator=(const Iresidue& res);
		Iresidue(Iresidue&& res) noexcept;
		Iresidue& operator=(Iresidue&& res) noexcept;
		~Iresidue();
		
};
//===================================================
/**
 * @class Iprotein
 * @author igor
 * @date 09/03/20
 * @file Iprotein.h
 * @brief 
 */
class Iprotein{
	public:
		std::string name;
		std::string remark;
		std::string title;
		unsigned int num_of_res;
		unsigned int lig_num;
		bool ligand;
		std::vector<int> p1;
		std::vector<int> p2;
		std::vector<Iresidue> residues;
		std::vector<double> xcoord;
		std::vector<double> ycoord;
		std::vector<double> zcoord;
		std::vector<double> b_factor;
		//---------------------------------------------------------
		Iprotein();
		Iprotein (const char* pdb_name);
		Iprotein(const Iprotein& prot_rhs);
		Iprotein(Iprotein&& prot_rhs) noexcept;
		Iprotein& operator=(const Iprotein& prot_rhs);
		Iprotein& operator=(Iprotein&& prot_rhs) noexcept;
		//---------------------------------------------------------
		void print();
		void load_b_column(std::vector<double>& b_fact);
		~Iprotein();
	
};
//====================================================
/**
 * @class pdb
 * @author igor
 * @date 09/03/20
 * @file Iprotein.h
 * @brief 
 */
class pdb{
	public:
		unsigned int nModels;
		std::string name;
		std::vector<Iprotein> models;
		//---------------------------------------------------------
		pdb();
		pdb(const pdb& rhs_pdb);
		pdb(pdb&& rhs_pdb) noexcept;
		pdb& operator=(pdb& rhs_pdb);
		pdb& operator=(pdb&& rhs_pdb) noexcept;
		pdb(const char* file_name); // to read multi pdbs
		pdb(Iprotein& prot);
		~pdb(){};	
		//---------------------------------------------------------
		void add_protein(Iprotein& prot);
		Iprotein& get_model(unsigned int i);
		void write_models(std::string path);
		void write_pdb(std::string fname);
};


#endif 