#include <iostream>
#include <iomanip>
#include <armadillo>
#include <map>
#include <string>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <float.h>
#include <omp.h>

#include "galaxy.h"
#include "neighbours.h"
#include "util.h"
#include "definitions.h"

//git

using namespace std;

galaxy::galaxy(string pfile)
{
	iname = pfile;
	Read_params(pfile);
	

	MAX_SFR_INDEX=10; // length of SFH storage
	RSFH = fmat(prm.at("Rings_i").i,MAX_SFR_INDEX,fill::zeros);
	
	MAX_OTFL_INDEX=10; // length of SFH storage
	ROTFL_GAS = fmat(prm.at("Rings_i").i,MAX_OTFL_INDEX,fill::zeros);
	ROTFL_METALS = fmat(prm.at("Rings_i").i,MAX_OTFL_INDEX,fill::zeros);

	GLOB_OUTFLOW_MGAS=0;
	GLOB_OUTFLOW_METALS=0;
	
	
	//==================================================================
	// GLX RADIAL PHOT
	//------------------------------------------------------------------
	//prm["PEG_RINGS_i"].i=1; // <- all galaxy is one ring
	//GlxRadPhot = fmat(prm.at("PEG_RINGS_i").i, gp.size(),fill::zeros);
	//==================================================================
	if(prm["OUTPUT_TYPES_n"].dict.at("0d"))
	{
		GLX["SFR"]=0;
		GLX["ACCRETION"]=0;
		GLX["SPONTANEOUS_EVENTS"]=0;
		GLX["TRIGGERED_EVENTS"]=0;
	}
	// mappings of GlxRings columns
	gr["RadPc"]=0;  //radius in kpc
	gr["AziCor"]=1; //azimuthal coordinate of first cell in Rad units
	gr["RotRad"]=2; //rotation speed in rad per timestep
	gr["Mgas"]=3;   //gas content in the cell
	gr["Zgas"]=4;   //gas metallicity in the cell
	gr["CellW"]=5;  // CellWidth
	gr["TMgas"]=6;	// Total Mgas
	gr["AccEf"]=7;	// normalized accretion efficiency
	gr["SF_event"]=8;
	gr["ACTIVE"]=9;
	gr["SPONT_event"]=10;
	gr["SFR"]=11;
	gr["EqGrad"]=12;
	gr["O_mgas"]=13;
	gr["O_metals"]=14;

	//mappings for GlxCells columns
	gc["Mgas"]=0;
	gc["Zgas"]=1;
	gc["RefTime"]=2;
	gc["Diff_gas_rez"]=3;	//diffusion reservour for gas
	gc["Diff_zgas_rez"]=4;	//diffusion reservour for zgas
	gc["SPONT_BUFFER"]=5;	//to speed up spontaneous probability evaluation
	gc["REF_TIME_BUFFER"]=6;//to prevent multiple reevaluation of refractory time
	gc["SFR_TRESH"]=7;//avreaged over neighbours last delta-SF  time
	gc["SFE"]=8;
	gc["last_Mstr"]=9;
	gc["Metals"]=10;
	gc["TVEL"]=11;
	gc["rho"]=12;
	gc["O_mgas"]=13;
	gc["DiffN"]=14;
	gc["TrigTime"]=15;
	gc["OflN0"]=16;
	

	// initialise random number generator
	gsl_type =gsl_rng_mt19937;
	random_number_gsl = gsl_rng_alloc (gsl_type);
	gsl_rng_set(random_number_gsl, prm.at("RND_SEED_i").i);

	NCells=Col<unsigned short>(prm.at("Rings_i").i);
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		NCells(i)=i*6;
	}
	NCells(0)=1;
	//==================================================================
	// GLX CELLS
	//------------------------------------------------------------------
	GlxCells = field<mat>(prm.at("Rings_i").i);
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxCells(i)=mat(NCells(i),gc.size(),fill::zeros);
		GlxCells(i).col(gc.at("RefTime")).fill(1000);
		GlxCells(i).col(gc.at("TrigTime")).fill(1000);
	}
	//==================================================================
	// GLX PHOT
	//------------------------------------------------------------------
	GlxPhot = field<fmat>(prm.at("Rings_i").i);
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxPhot(i)=fmat(NCells(i), gp.size(),fill::zeros);
	}
	//==================================================================
	// GLX STARS
	//------------------------------------------------------------------
	GlxStars = field<fcube>(prm.at("Rings_i").i);
	GlxStarsBuff=field<fcube>(prm.at("Rings_i").i);
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxStars(i)=fcube(prm.at("SSP_matrix_age_i").i,prm.at("SSP_matrix_z_i").i, NCells(i), fill::zeros);
		GlxStarsBuff(i)=fcube(prm.at("SSP_matrix_age_i").i,prm.at("SSP_matrix_z_i").i, NCells(i), fill::zeros);
	}
	StarsFlowSpeed=fmat(prm.at("SSP_matrix_age_i").i-1,prm.at("SSP_matrix_z_i").i);
	for(int i=0; i<prm.at("SSP_matrix_z_i").i; ++i)
		StarsFlowSpeed.col(i)=prm.at("Time_step_Myr_i").i/conv_to<Col<float> >::from(prm.at("SSP_ages_usv").usv.rows(1,prm.at("SSP_ages_usv").usv.n_rows-1)-
	prm.at("SSP_ages_usv").usv.rows(0,prm.at("SSP_ages_usv").usv.n_rows-2));
	//AgeFlow=mat(prm.at("SSP_matrix_age_i").i-1,prm.at("SSP_matrix_z_i").i);
	//==================================================================
	// GLX RINGS
	//------------------------------------------------------------------
	GlxRings = fmat(prm.at("Rings_i").i, gr.size(), fill::zeros);
	GlxRings.col(gr.at("RadPc"))=linspace<fvec>(0,(prm.at("Rings_i").i-1)*prm.at("Cell_size_pc_f").f,prm.at("Rings_i").i);
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxRings(i,gr.at("AziCor"))=2*datum::pi*gsl_rng_uniform (random_number_gsl);
	}
	//to ensure proper behaviour of neighbor finding algorythm
	GlxRings(0,gr.at("AziCor"))=GlxRings(0,gr.at("AziCor"));
	//GlxRings.col(gr["AziCor"]).fill(0);
	read_rotation_curve(); //RotRad
	GlxRings.col(gr.at("Mgas")) = exp(-GlxRings.col(gr.at("RadPc"))/2000.);
	GlxRings.col(gr.at("Zgas")).fill(-1);
	GlxRings.col(gr.at("CellW"))=2.0*datum::pi/(linspace<fvec>(0,(prm.at("Rings_i").i-1),prm.at("Rings_i").i)*6);
	GlxRings(0,gr.at("CellW"))=2*datum::pi;
	//GlxRings.col(gr.at("TMgas"))=prm.at("ACC_Amp_0_f").f*exp(-GlxRings.col(gr.at("RadPc"))/prm.at("ACC_Re_0_f").f)+
	//(1.-prm.at("ACC_Amp_0_f").f)*exp(-GlxRings.col(gr.at("RadPc"))/prm.at("ACC_Re_1_f").f);
	//==================================================================
	// NEIGHBOURS
	//------------------------------------------------------------------
	//initialise neighbours coordinates matrix
	GlxNeigh = new neighbours* [prm.at("Rings_i").i];
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxNeigh[i] = new neighbours(i);
	}
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		if(i==0)
			GlxNeigh[i]->init_Up_Lo_adr(GlxNeigh[i+1], NULL);
		else if(i==prm.at("Rings_i").i-1)
			GlxNeigh[i]->init_Up_Lo_adr(NULL, GlxNeigh[i-1]);
		else
			GlxNeigh[i]->init_Up_Lo_adr(GlxNeigh[i+1], GlxNeigh[i-1]);
	}
	//==================================================================
	
	TOTAL_METALS=0;
	
	InitInfallScenario();

	if(prm.at("OUTPUT_TYPES_n").dict.at("0d"))
	{
		string igal_str = iname+"_igal.dat";
		igal.open(igal_str.c_str(), ios::out|ios::app);
	}
	if(prm.at("OUTPUT_TYPES_n").dict.at("pegase"))
	{
		string ipeg_str = iname+"_ipeg.dat";
		ipeg.open(ipeg_str.c_str(), ios::out|ios::app);
	}
	if(prm.at("OUTPUT_TYPES_n").dict.at("ring"))
	{
		stringstream iring_str;
		iring_str <<iname <<"_i"<<prm.at("RING_i").i<<+"ring.dat";
		iring.open(iring_str.str().c_str(), ios::out|ios::app);
	}

	if(prm["INIT_SFR_b"].initialized )
	{
		cerr<<"INTI_SFR"<<endl;
		map<unsigned short,unsigned short>::iterator it;
		for(it=prm.at("INIT_SFR_b").idict.begin(); it!=prm.at("INIT_SFR_b").idict.end(); ++it)
		{
			GlxStars(it->first)(0,0,it->second)=1e9;
			SF_event(it->first, it->second);
		}
		cerr<<"INTI_SFR finished"<<endl;
	}
}
galaxy::~galaxy()
{
	GlxRings.clear();
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		delete GlxNeigh[i];
		GlxCells(i).clear();
		GlxStars(i).clear();

	}
	delete GlxNeigh;
	prm.erase(prm.begin(),prm.end());

	if(iring.is_open())
		iring.close();
	if(igal.is_open())
		igal.close();
	if(ipeg.is_open())
		ipeg.close();
}

void galaxy::flush_all_dict(string fname)
{
	fstream log;
	log.open(fname.c_str(),ios::out);
	// show content:
  for (std::map<string,parameter>::iterator it=prm.begin(); it!=prm.end(); ++it)
  {
    if (it->first[it->first.size()-1]=='f')
    {
		 log << it->first <<  " => " << it->second.f << endl;
	}
	else if (it->first[it->first.size()-1]=='i')
	{
		log << it->first <<  " => " << it->second.i << endl;
	}
	else if (it->first[it->first.size()-1]=='s')
	{
		log << it->first <<  " => " << it->second.str << endl;
	}
	else if (it->first[it->first.size()-1]=='v')
	{
		log << it->first <<  " => ";
		if (it->first[it->first.size()-2]=='f')
		{
			for(unsigned int i=0; (i<it->second.fv.size()) && (i<5); i++)
				log<< it->second.fv(i) << " ";
			if(it->second.fv.size() >= 5)
				log <<"... " << it->second.fv(it->second.fv.size()-1);
			log << endl;
		}
		if ((it->first[it->first.size()-2]=='s') && (it->first[it->first.size()-3]=='u'))
		{
			for(unsigned int i=0; (i<it->second.usv.size()) && (i<5); i++)
				log<< it->second.usv(i) << " ";
			if(it->second.usv.size() >= 5)
				log <<"... " << it->second.usv(it->second.usv.size()-1);
			log << endl;
		}	
	}
	else if (it->first[it->first.size()-1]=='n')
	{
		log << it->first <<  " => True" << endl;
	}
	else if (it->first[it->first.size()-1]=='c')
	{
		log << it->first <<  " => True" << endl;
	}
	else if (it->first[it->first.size()-1]=='a')
	{
		log << it->first <<  " => True" << endl;
	}
	else if (it->first[it->first.size()-1]=='b')
	{
		log << it->first <<  " => ";
		if(it->second.initialized)
			log << "True" << endl;
		else
			log << "False" << endl;
	}
	else
	{
		log << it->first <<  " => MISSING" << endl;
	}
  }
  log.close();
}

void galaxy::Read_params(string fname)
{
	fstream fin;
	fin.open(fname.c_str(),ios::in);
	if(!fin.is_open())
	{
		cerr<< "Input file:"<<fname<<" failed to open"<<endl;
		exit(EXIT_FAILURE);
	}
	
	string in_str;
	for(;;)
	{
		fin >> in_str;

		if(in_str.size() == 0 || fin.eof())
		{
			break;
		}
		// check if line contents are comments, otherwise interpret
		else if(in_str.find_first_of("#") < in_str.npos )
            getline(fin, in_str);
        else if(in_str.find_first_of("#") == in_str.npos )
		{
			// ----------------------------------------
			// SFR 'directly' controlling parameters
			if(in_str==string("SFE_POW"))
            {
				string check_input;
				fin >> check_input;
				if( isdigit(check_input[0]) )
				{
					cerr << "SFE_POW set to constant value" << endl;
					prm["SFE_POW_f"].f = atof(check_input.c_str());
					prm["SFE_POW_f"].initialized=true;
				}
				else
				{
					cerr << "SFE_POW set as time variable" << endl;
					fstream input_file;
					input_file.open(check_input.c_str(),ios::in);
					strip_header(input_file);
					prm["SFE_POW_f"].fv.load(input_file,raw_ascii);
					cerr<< "SFE_POW entries:" << prm["SFE_POW_f"].fv.n_rows << endl;
					input_file.close();
				}
			}
			else if(in_str==string("SFE"))
            {
				string check_input;
				fin >> check_input;
				if( isdigit(check_input[0]) )
				{
					cerr << "SFE set to constant value" << endl;
					prm["SFE_f"].f = atof(check_input.c_str());
					prm["SFE_f"].initialized=true;
				}
				else
				{
					cerr << "SFE set as time variable" << endl;
					fstream input_file;
					input_file.open(check_input.c_str(),ios::in);
					strip_header(input_file);
					prm["SFE_f"].fv.load(input_file,raw_ascii);
					cerr<< "SFE entries:" << prm["SFE_f"].fv.n_rows << endl;
					input_file.close();
				}
			}
			else if(in_str==string("TRIGGERED"))
			{
				string check_input;
				fin >> check_input;
				if( isdigit(check_input[0]) )
				{
					cerr << "TRIGGERED_PROBABILITY set to constant value" << endl;
					prm["TRIG_PROB_f"].f = atof(check_input.c_str());
					prm["TRIG_PROB_f"].initialized=true; 
				}
				else
				{
					cerr << "TRIGGERED_PROBABILITY set as time variable" << endl;
					fstream input_file;
					input_file.open(check_input.c_str(),ios::in);
					strip_header(input_file);
					prm["TRIG_PROB_f"].fv.load(input_file,raw_ascii);
					cerr<< "TRIGG entries:" << prm["TRIG_PROB_f"].fv.n_rows << endl;
					input_file.close();
				}
				fin >> prm["TRIG_TIME_i"].i;
			}
			//-----------------------------------------
			else if(in_str==string("SSP"))
			{
				string fage, fz;
				fin >> fage >> fz;
				prm["SSP_ages_usv"].usv.load(fage,raw_ascii);
				prm["SSP_metalls_fv"].fv.load(fz,raw_ascii);

			   	prm["SSP_matrix_age_i"].i = prm.at("SSP_ages_usv").usv.n_rows;
			   	prm["SSP_matrix_z_i"].i = prm.at("SSP_metalls_fv").fv.n_rows;
            }
            else if(in_str==string("ACC_METALLICITY"))
            {
				fin >> ACCRETION_METTALICITY;
			}
			else if(in_str==string("NRM_SFE"))
            {
				fin >> prm["NRM_SFE_f"].f;
				prm.at("NRM_SFE_f").f=prm.at("NRM_SFE_f").f*prm.at("CELL_AREA_f").f;
			}
			else if(in_str==string("MINIMUM_SFE"))
            {
				fin >> prm["MINIMUM_SFE_f"].f;
			}
			else if(in_str=="MAXIMUM_SFE")
			{
				fin >> prm["MAXIMUM_SFE_f"].f;
			}
			else if(in_str==string("TRIGG_MASS"))
            {
				fin >> prm["MINIMUM_MASS_f"].f >> prm["TRIGG_MASS_f"].f ;
			}
			else if(in_str==string("GALAXY_AGE"))
			{
				fin >> prm["Galaxy_age_Myr_i"].i;
			}
			else if(in_str==string("TIMESTEP"))
			{
				fin >> prm["Time_step_Myr_i"].i;
			}
			else if(in_str==string("SEED"))
			{
				fin >> prm["RND_SEED_i"].i;
			}
			else if(in_str==string("OUTPUT"))
			{
				fin >> prm["OUTPUT_TYPES_n"].n;
				//set initial values to false
				prm.at("OUTPUT_TYPES_n").dict["0d"]=false;
				prm.at("OUTPUT_TYPES_n").dict["1d"]=false;
				prm.at("OUTPUT_TYPES_n").dict["2d"]=false;
				prm.at("OUTPUT_TYPES_n").dict["cmd"]=false;
				prm.at("OUTPUT_TYPES_n").dict["pegase"]=false;
				prm.at("OUTPUT_TYPES_n").dict["ring"]=false;

				string tmp_string;
				for(int i=0; i<prm.at("OUTPUT_TYPES_n").n; ++i)
				{
					fin>>tmp_string;
					prm.at("OUTPUT_TYPES_n").dict.at(tmp_string)=true;
				}
				
				string check_input;
				fin >> check_input;
				std::vector<unsigned short> t;
				if( isdigit(check_input[0]) ) {
					prm["Output_Times_ia"].an=atof(check_input.c_str());
					unsigned short tin;
					for(int i=0; i<prm.at("Output_Times_ia").an; ++i)
					{
						fin >>tin;
						cerr << tin << endl;
						t.push_back(tin);
					}
					prm["Output_Times_ia"].usv = Col<unsigned short>(t);
					prm.at("Output_Times_ia").cur=0;
					t.clear();
				}
				else{
					unsigned short tin0, tin1, t_step;
					fin >>tin0 >> tin1 >> t_step;
					prm["Output_Times_ia"].an=(tin1-tin0)/t_step+1;
					for(int i=0; i<prm.at("Output_Times_ia").an; ++i)
					{
						cerr << tin0+t_step*i << endl;
						t.push_back(tin0+t_step*i);
					}
					prm["Output_Times_ia"].usv = Col<unsigned short>(t);
					prm.at("Output_Times_ia").cur=0;
					t.clear();
				}
			}
			// Used for ouput detailed information of selected ring
			else if(in_str==string("RING"))
			{
				fin >> prm["RING_i"].i;
				cerr<<"RING for output:" << prm["RING_i"].i << endl;
			}

			else if(in_str==string("GAS_DIFFUSION"))
			{
				fin >> prm["Diffusion_const_f"].f;
				if (prm.at("Diffusion_const_f").f != 0)
				{
					prm.at("Diffusion_const_f").f = 1./(prm.at("Diffusion_const_f").f/prm.at("Time_step_Myr_i").i);
				}
			}
			else if(in_str==string("ROTATION"))
			{
				fin >> prm["Rotation_curve_s"].str;
			}
			else if(in_str==string("GRID_SIZE"))
			{
				fin >> prm["Rings_i"].i >> prm["Rings_BUFF_i"].i;
			}
			else if(in_str==string("SN_TIMESCALE"))
			{
				fin >> prm["SN_Timescale_10Myr_i"].i;
				prm["SN_Timescale_10Myr_i"].i /= prm.at("Time_step_Myr_i").i;
			}
			else if(in_str==string("CELL_SIZE"))
			{
				fin >>prm["Cell_size_pc_f"].f;
				prm["CELL_AREA_f"].f=pow(prm.at("Cell_size_pc_f").f,2)*datum::pi/3;
			}
			else if(in_str==string("PEGASE"))
			{
				fin >>prm["MGAS_PEGASE_s"].str >> prm["ZGAS_PEGASE_s"].str;
				MGAS_PEGASE.load(prm.at("MGAS_PEGASE_s").str, raw_ascii);
				ZGAS_PEGASE.load(prm.at("ZGAS_PEGASE_s").str, raw_ascii);
			}
			else if(in_str==string("PHOTOMETRY"))
			{
				prm["PHOTOMETRY_sc"].initialized=true;
				unsigned short idx;
				fin >>idx;
				string loc, passb;
				prm["PHOTOMETRY_sc"].as=field<std::string>(idx);
				for(unsigned short int i=0; i<idx; ++i)
				{
					fin >> loc >> prm.at("PHOTOMETRY_sc").as(i);

					PHOT_PEGASE[prm.at("PHOTOMETRY_sc").as(i)].load(loc,raw_ascii);
					gp[prm.at("PHOTOMETRY_sc").as(i)]=i;
				}
			}
			//==========================================================
			// SPONTANEOUS SF PARAMETERS
			else if(in_str==string("GAS2_SPONT"))
			{
				SPONT_TYPE=0;
				fin >> prm["GAS2_SPONT_NRM_f"].f ;
			}
			else if(in_str==string("COUNT_ACTIVE"))
			{
				prm["COUNT_ACTIVE"].initialized=true ;
			}
			//==========================================================
			else if(in_str==string("GAS_SFR_TRESHOLD"))
			{
				fin >>prm["GAS_SFR_TRESHOLD_f"].f;
			}
			else if(in_str==string("SHARP_SFTRESHOLD"))
			{
				prm["SHARP_SFTRESHOLD"].initialized=true;
			}
			else if(in_str==string("STELLAR_DIFFUSION"))
			{
				fin >>prm["STELLAR_DIFUSSION_f"].f >> prm["Minimum_Diff_Mass_f"].f;
			}
			else if(in_str==string("INIT_SFR"))
			{
				prm["INIT_SFR_b"].initialized=true;
				
				int idx, i_ring, i_cell;
				fin >> idx;
				for(int i=0; i<idx; ++i)
				{
					fin >> i_ring >> i_cell;
					cerr<< i_ring << " " << i_cell << endl;
					prm.at("INIT_SFR_b").idict[i_ring]=i_cell;
				}
			}
			else if(in_str==string("NUM_OF_THREADS"))
			{
				fin>>prm["NUM_OF_THREADS_i"].i;
				omp_set_num_threads(prm.at("NUM_OF_THREADS_i").i);
			}
			else if(in_str==string("CMD_LIMIT"))
			{
				//first one -- mass limit, second -- age
				fin>>prm["CMD_MASS_LIMIT_f"].f >> prm["CMD_AGE_LIMIT_f"].f;
				prm.at("CMD_MASS_LIMIT_f").initialized=true;
				prm.at("CMD_AGE_LIMIT_f").initialized=true;
			}
			else if(in_str==string("ACCRETION"))
			{
				prm["ACCRETION_b"].initialized=true;
				fin>>prm["ACCRETION_RAD_PROFILE_s"].str >> prm["ACCRETION_TIME_PROFILE_s"].str;
			}
			else if(in_str==string("DISTANCE"))
			{
				fin >> prm["DISTANCE_MPC_f"].f;
				PC2_TO_ARCSEC2 = pow(atan(1./(prm["DISTANCE_MPC_f"].f*1e6))*180./datum::pi*3600, 2);//0.060296 --> 840 KPC;
			}
			else if(in_str==string("OUTFLOW"))
			{
				fin >>prm["OUTFLOW_FRAC_f"].f >> prm["OUTFLOW_VEL_RATIO_f"].f;
				prm["OUTFLOW_b"].initialized=true;
			}
		}
	}
	cerr <<"Conversion PC^2 --> ARCSEC^2: " <<  PC2_TO_ARCSEC2 << endl;
	flush_all_dict(fname+".log");
	fin.close();
}


/*
 *  The function calculates gas and metals added to the cell by evolving SSPs within the cell.
 *  According to the outflow parameter calculated in SF_event() funtion, fraction of metals
 *  and gas are discarded here.
 * 
 *  The function operates on internal class variables:
 *  GlxCells
 *  GlxRings
 * 
 * 	Updates output variables:
 *  ROTFL_GAS(i,OTFL_INDEX)
 *  ROTFL_METALS(i,OTFL_INDEX)
 * 
 *  Updates consistency check variables:
 * 	GLOB_OUTFLOW_MGAS
 * 	TOTAL_METALS
 * 	GLOB_STELLAR_ACC
 *  
 */
void galaxy::EvolveMetals()
{
	float LEFT_AFTER_OUTFLOW=1.-prm.at("OUTFLOW_FRAC_f").f;
	double SUM=0;
	double METALS=0;
	double STELLAR_ACC=0;
	#pragma omp parallel for schedule(guided) reduction(+:SUM,METALS, STELLAR_ACC)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		float eMgas=0, eZgas=0, delta_Mgas=0, delta_Zgas=0;
		for(unsigned short j=0; j<NCells(i); ++j)
		{
			// returned gas from evolved SPPs
			eMgas = accu(MGAS_PEGASE % GlxStars(i).slice(j));
			eZgas = accu(ZGAS_PEGASE % GlxStars(i).slice(j));
			
			//  case of star formation event later then SN timescale
			if(GlxCells(i)(j,gc.at("RefTime")) > prm.at("SN_Timescale_10Myr_i").i)
			{
				GlxCells(i)(j, gc.at("Metals"))+=eZgas;
				GlxCells(i)(j, gc.at("Mgas"))+= eMgas;
				GlxCells(i)(j, gc.at("Zgas"))= GlxCells(i)(j, gc.at("Metals")) / GlxCells(i)(j, gc.at("Mgas"));
				METALS+=eZgas;
				STELLAR_ACC+=eMgas;
			}
			//  case of star formation event earlier then SN timescale
			else if(GlxCells(i)(j,gc.at("RefTime")) >= 0 && GlxCells(i)(j,gc.at("RefTime")) <= prm.at("SN_Timescale_10Myr_i").i)
			{
				if ( GlxCells(i)(j,gc.at("TVEL")) > 1 && prm["OUTFLOW_b"].initialized)
				{
					delta_Mgas=eMgas*prm.at("OUTFLOW_FRAC_f").f;
					delta_Zgas=eZgas*prm.at("OUTFLOW_FRAC_f").f;
	
					GlxRings(i, gr.at("O_metals"))+=delta_Zgas;
					GlxRings(i, gr.at("O_mgas"))+= delta_Mgas;
					
					SUM+=delta_Mgas;
					
					ROTFL_GAS(i,OTFL_INDEX)+=delta_Mgas;
					ROTFL_METALS(i,OTFL_INDEX)+=delta_Zgas;
					
					GlxCells(i)(j, gc.at("Metals"))+=LEFT_AFTER_OUTFLOW*eZgas;
					GlxCells(i)(j, gc.at("Mgas"))  +=LEFT_AFTER_OUTFLOW*eMgas;
					GlxCells(i)(j, gc.at("Zgas"))   = GlxCells(i)(j, gc.at("Metals")) / GlxCells(i)(j, gc.at("Mgas"));
					METALS+=LEFT_AFTER_OUTFLOW*eZgas;
					STELLAR_ACC+=LEFT_AFTER_OUTFLOW*eMgas;
				}
				else
				{
					GlxCells(i)(j, gc.at("Metals"))+=eZgas;
					GlxCells(i)(j, gc.at("Mgas"))+= eMgas;
					GlxCells(i)(j, gc.at("Zgas"))= GlxCells(i)(j, gc.at("Metals")) / GlxCells(i)(j, gc.at("Mgas"));
					METALS+=eZgas;
					STELLAR_ACC+=eMgas;
				}
			}
			// Consistency check
			if(GlxCells(i)(j, gc.at("Metals"))/GlxCells(i)(j, gc.at("Mgas"))<0.0001)
			{
				cerr<<"WHAT? i: "<<i << " j: " << j << " eZgas: "<<eZgas;
				cerr<<" eMgas: "<< eMgas << " Mgas: "<<GlxCells(i)(j, gc.at("Mgas")) << endl;
				cerr<<" delta_mgas: "<< delta_Mgas << " delta_Zgas: "<<delta_Zgas << endl;
				cerr<<" Zgas: "<<GlxCells(i)(j, gc.at("Zgas")) << endl;
				cerr << " outflow: "<<  GlxCells(i)(j,gc.at("TVEL")) << " ref time: " << GlxCells(i)(j,gc.at("RefTime")) <<endl;
			}
		}
	}
	GLOB_OUTFLOW_MGAS=SUM;
	TOTAL_METALS+=METALS;
	GLOB_STELLAR_ACC=STELLAR_ACC;
}

/*
 *  The function shifts matrix of stellar mases in age direction.
 *  This is executed for every cell in the model
 * 
 *  The function operates on internal class variables:
 *  GlxStars
 * 
 *  The function uses temporary variable AgeFlow
 * 
 * 	Updates output variables:
 *  GlxStars
 */
void galaxy::EvolveAge()
{
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		//allocate matrix which will be used for shifts
		fmat AgeFlow = fmat(prm.at("SSP_matrix_age_i").i,prm.at("SSP_matrix_z_i").i,fill::zeros);
		for(unsigned short j=0; j<NCells(i); ++j)
		{
			AgeFlow(span(1,prm.at("SSP_matrix_age_i").i-1),span::all) = GlxStars(i).slice(j)(span(0,prm.at("SSP_matrix_age_i").i-2),span::all);
			GlxStars(i).slice(j)=AgeFlow;
		}
	}
}

/*
 *  The function moves gas mass between cells
 * 
 *  The function operates on internal class variables:
 *  GlxCells
 * 
 * 	Updates consistency check variables:
 *  METALS_REMAIN
 *  MGAS_REMAIN
 */
void galaxy::gas_diffusion()
{
	#pragma omp parallel for
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxCells(i).col(gc.at("Diff_gas_rez")).zeros();
		GlxCells(i).col(gc.at("Diff_zgas_rez")).zeros();
		GlxCells(i).col(gc.at("DiffN")).zeros();
	}
	
	//to avoid concurent writing into memory location, gas echange is splitted
	violent_gas_exchange(0, 3);
	violent_gas_exchange(1, 3);
	violent_gas_exchange(2, 3);
	
	calm_gas_exchange(0, 2);
	calm_gas_exchange(1, 2);
	
	// safety checks in principle can be disabled
	
	//safety check 1
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i; i++)
		for(int j=0; j<NCells(i); ++j)
		{
			if(GlxCells(i)(j,gc.at("DiffN"))!=GlxNeigh[i]->Num(j))
			{
				#pragma omp critical
				{
				cerr<<endl;
				cerr<<"TIMESTEP: "<< TIMESTEP << " reftime: "<< GlxCells(i)(j,gc.at("RefTime"))<<endl;
				cerr<<"!Neighbors exchange ERROR! "<<GlxNeigh[i]->Num(j)<<" "<<GlxCells(i)(j,gc.at("DiffN"))<<" "<<i<<" "<<j<<endl;
				for(int n=0; n<GlxNeigh[i]->Num(j); ++n)
				{
					cerr<<GlxNeigh[i]->Rings(j,n)<<" "<<GlxNeigh[i]->Cells(j,n)<<" "<<GlxNeigh[i]->Areas(j,n)<<" ";
					cerr<<GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("RefTime"))<<" ";
					cerr<<GlxNeigh[GlxNeigh[i]->Rings(j,n)]->Num(GlxNeigh[i]->Cells(j,n))<<" ";
					cerr<<GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("DiffN"))<< endl;
				}
				cerr<<"Shortlist"<<endl;
				for(int n=0; n<GlxNeigh[i]->oN(j); ++n)
				{
					cerr<<GlxNeigh[i]->oR(j,n)<<" "<<GlxNeigh[i]->oC(j,n)<<" "<<GlxNeigh[i]->oA(j,n)<<" ";
					cerr<<GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("RefTime"))<<" ";
					cerr<<GlxNeigh[GlxNeigh[i]->oR(j,n)]->Num(GlxNeigh[i]->oC(j,n))<<" "<<GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("DiffN"))<< endl;
				}
				cerr<<"Details"<<endl;
				for(int n=0; n<GlxNeigh[i]->Num(j); ++n)
				{
					cerr<<GlxNeigh[i]->Rings(j,n)<<" ## "<<GlxNeigh[i]->Cells(j,n)<<endl;
					for(int nn=0; nn<GlxNeigh[GlxNeigh[i]->Rings(j,n)]->oN(GlxNeigh[i]->Cells(j,n)); ++nn)
					{
						cerr<<GlxNeigh[GlxNeigh[i]->Rings(j,n)]->oR(GlxNeigh[i]->Cells(j,n),nn)<<" "<<GlxNeigh[GlxNeigh[i]->Rings(j,n)]->oC(GlxNeigh[i]->Cells(j,n),nn)<<endl;
					}
				}
				}
			}
		}
	//safety check 2
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
		for(int j=0; j<NCells(i); ++j)
		{
			if( ( GlxCells(i)(j, gc.at("Metals"))+GlxCells(i)(j,gc.at("Diff_zgas_rez")))/
			(GlxCells(i)(j, gc.at("Mgas"))+GlxCells(i)(j,gc.at("Diff_gas_rez"))   )<0.0001)
			{
				#pragma omp critical
				{
				cerr<<"TIMESTEP" << TIMESTEP << endl;
				cerr<<"WHAT?  GasDiff Zgas BEFORE: "<<GlxCells(i)(j, gc.at("Metals"))/GlxCells(i)(j, gc.at("Mgas")) << endl;
				cerr<<" Zgas AFTER: "<<( GlxCells(i)(j, gc.at("Metals"))+GlxCells(i)(j,gc.at("Diff_zgas_rez")))/
			(GlxCells(i)(j, gc.at("Mgas"))+GlxCells(i)(j,gc.at("Diff_gas_rez"))) << endl;
				cerr<<" Mgas BEFORE: "<<GlxCells(i)(j, gc.at("Mgas")) << endl;
				cerr<<" Mgas after: "<<GlxCells(i)(j, gc.at("Mgas"))+GlxCells(i)(j,gc.at("Diff_gas_rez")) << endl;
				cerr<<" Metals BEFORE: "<<GlxCells(i)(j, gc.at("Metals")) << endl;
				cerr<<" Metals after: "<<GlxCells(i)(j, gc.at("Metals"))+GlxCells(i)(j,gc.at("Diff_zgas_rez")) << endl;
				cerr<<" ref time: "<<GlxCells(i)(j,gc.at("RefTime")) <<endl;
				}
			}
		}
			
	// results of gas exchange are writen into arrays
	double METALS_REMAIN=0, MGAS_REMAIN=0;
	#pragma omp parallel for schedule(guided) reduction(+:METALS_REMAIN,MGAS_REMAIN)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxCells(i).col(gc.at("Metals"))+=GlxCells(i).col(gc.at("Diff_zgas_rez"));
		GlxCells(i).col(gc.at("Mgas"))+= GlxCells(i).col(gc.at("Diff_gas_rez"));
		GlxCells(i).col(gc.at("Zgas"))=GlxCells(i).col(gc.at("Metals"))/GlxCells(i).col(gc.at("Mgas"));
		METALS_REMAIN+=sum(GlxCells(i).col(gc.at("Diff_zgas_rez")));
		MGAS_REMAIN+=sum(GlxCells(i).col(gc.at("Diff_gas_rez")));
	}
	// consistency check 3
	if(fabs(METALS_REMAIN)>1e-8 || fabs(MGAS_REMAIN) >1e-8)
	{
		cerr<<"TIMESTEP: "<<TIMESTEP<<" METALS: "<<METALS_REMAIN<<" MGAS :"<<MGAS_REMAIN<<endl;
	}
}

/*
 *  The function is subfunction of gas_diffusion()
 *  Takes 2 parameters:
 *  int start_idx -  index of the ring on which to execute
 *  int incr_idx  -  step in index, i.e. value 3 would mean
 *  to run function on every third ring - 0, 3, 6, 9
 * 
 *  The function operates on internal class variables (GlxCells),
 *  but writes only in temporary storage variables -- gc.at("Diff_gas_rez"), gc.at("Diff_zgas_rez").
 *  Actual writing occurs in the function gas_diffusion. 
 * 
 * 	Updates consistency check variables:
 *  GlxCells(i)(j,gc.at("DiffN"))
 */
void galaxy::violent_gas_exchange(int start_idx, int incr_idx)
{
	#pragma omp parallel for schedule(guided)
	for(int i=start_idx; i<prm.at("Rings_i").i; i+=incr_idx)
	{
		double diff_Mgas, diff_Zgas;		
		for(int j=0; j<NCells(i); ++j)
		{
			if(GlxCells(i)(j,gc.at("RefTime")) <= prm.at("SN_Timescale_10Myr_i").i)
			{
				for(int n=0; n<GlxNeigh[i]->Num(j); ++n)
				{
					// gas exchange in violent mode
					if(GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("RefTime")) > prm.at("SN_Timescale_10Myr_i").i)
					{
						diff_Mgas = GlxCells(i)(j,gc.at("Mgas"))*GlxNeigh[i]->Areas(j,n)*0.95;
						diff_Zgas = GlxCells(i)(j,gc.at("Metals"))*GlxNeigh[i]->Areas(j,n)*0.95;
						
						if(diff_Zgas/diff_Mgas < 0.0001)
						{
							cerr<<"How? mgas: "<<diff_Mgas<<" metals: "<<diff_Zgas<<" zgas: "<< diff_Zgas/diff_Mgas <<endl;
						}
						float dg=GlxCells(i)(j,gc.at("Diff_gas_rez"));
						float dz=GlxCells(i)(j,gc.at("Diff_zgas_rez"));

						GlxCells(i)(j,gc.at("Diff_gas_rez"))-= diff_Mgas;
						GlxCells(i)(j,gc.at("Diff_zgas_rez"))-= diff_Zgas;
						GlxCells(i)(j,gc.at("DiffN"))+=1;
						
						if(GlxCells(i)(j,gc.at("Diff_zgas_rez"))/GlxCells(i)(j,gc.at("Diff_gas_rez"))/0.0001 < 1)
						{
							#pragma omp critical(PRINT)
							{
								cerr<<"How rez 1? i,j: "<< i <<" "<<j<<" metals: "<<GlxCells(i)(j,gc.at("Diff_zgas_rez"))<<" mgas: "<<GlxCells(i)(j,gc.at("Diff_gas_rez"));
								cerr<<" zgas: "<< GlxCells(i)(j,gc.at("Diff_zgas_rez"))/GlxCells(i)(j,gc.at("Diff_gas_rez")) <<endl;
								cerr<<" dmgas: " << diff_Mgas << " dzgas: "<< diff_Zgas <<" "<< diff_Zgas/diff_Mgas << endl;
								cerr<<" mgas: "<<GlxCells(i)(j,gc.at("Mgas"))<< " metals: "<<GlxCells(i)(j,gc.at("Metals"))<<" Z: "<<GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas"))<<endl;
								cerr<<GlxCells(i)(j,gc.at("DiffN"))<<endl;
								cerr<<dg << " "<<dz<<endl;
							}
						}
						GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("Diff_gas_rez"))+= diff_Mgas;
						GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("Diff_zgas_rez")) += diff_Zgas;
						GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("DiffN"))+=1;
					}
					// if neighbour also in violent mode, gas exchange does not occur
					else
					{
						GlxCells(i)(j,gc.at("DiffN"))+=1;
					}
				}
			}
		}
	}
}

/*
 *  The function is subfunction of gas_diffusion()
 *  Takes 2 parameters:
 *  int start_idx -  index of the ring on which to execute
 *  int incr_idx  -  step in index, i.e. value 3 would mean
 *  to run function on every third ring - 0, 3, 6, 9
 * 
 *  The function operates on internal class variables (GlxCells),
 *  but writes only in temporary storage variables -- gc.at("Diff_gas_rez"), gc.at("Diff_zgas_rez").
 *  Actual writing occurs in the function gas_diffusion. 
 * 
 * 	Updates consistency check variables:
 *  GlxCells(i)(j,gc.at("DiffN"))
 */
void galaxy::calm_gas_exchange(int start_idx, int incr_idx)
{
	#pragma omp parallel for schedule(guided)
	for(int i=start_idx; i<prm.at("Rings_i").i; i+=incr_idx)
	{
		double diff_Mgas, diff_Zmas, diff_nZmas, diff_nMgas, nrm_Mgas, DGAS, DZMAS;
		for(int j=0; j<NCells(i); j++)
		{
			if(GlxCells(i)(j,gc.at("RefTime")) > prm.at("SN_Timescale_10Myr_i").i)
			{
				for(int n=0; n<GlxNeigh[i]->oN(j); n++)
				{
					// both cells in quiescent mode -- diffusion
					if(GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("RefTime")) > prm.at("SN_Timescale_10Myr_i").i)
					{			
						nrm_Mgas = (GlxRings(         i          ,gr.at("AccEf")) + GlxRings(GlxNeigh[i]->oR(j,n),gr.at("AccEf")));

						diff_Mgas  = GlxCells(         i          )(          j          , gc.at("Mgas"))*
						             GlxRings(GlxNeigh[i]->oR(j,n),gr.at("AccEf"));
						
						diff_nMgas = GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n) , gc.at("Mgas"))*
									 GlxRings(         i          ,gr.at("AccEf"));
						
						diff_Zmas = GlxCells(         i          )(          j          , gc.at("Metals"))*
						             GlxRings(GlxNeigh[i]->oR(j,n),gr.at("AccEf"));
						             
						diff_nZmas = GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n) , gc.at("Metals"))*
									 GlxRings(         i          ,gr.at("AccEf"));
						
						DGAS =( diff_nMgas - diff_Mgas )/nrm_Mgas*GlxNeigh[i]->oA(j,n)*prm.at("Diffusion_const_f").f;
						
						DZMAS=( diff_nZmas - diff_Zmas )/nrm_Mgas*GlxNeigh[i]->oA(j,n)*prm.at("Diffusion_const_f").f;
						
						GlxCells(i)(j,gc.at("Diff_gas_rez"))+= DGAS;
						GlxCells(i)(j,gc.at("Diff_zgas_rez"))+= DZMAS;						
						GlxCells(i)(j,gc.at("DiffN"))+=1;
						
						GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("Diff_gas_rez"))-= DGAS;
						GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("Diff_zgas_rez"))-= DZMAS;
						GlxCells(GlxNeigh[i]->oR(j,n))(GlxNeigh[i]->oC(j,n),gc.at("DiffN"))+=1;
					}
				}
			}
		}
	}
}

/*
 *  The funtion for the calculation of stellar diffusion
 * 
 *  At first results are stored in the temporary array GlxStarsBuff,
 *  and then written into GlxStars.
 * 
 *  The parameter  prm.at("Minimum_Diff_Mass_f").f controls the treshold mass for the
 *  diffusion to proceed.
 */
void galaxy::stars_diffusion()
{
	int max_age=prm.at("SSP_matrix_age_i").i-1;
	#pragma omp parallel for
	for(int i=0; i<prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i; ++i)
	{
		GlxStarsBuff(i).zeros();
	}
	//work on even rings, this will avoid memory collision
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i; i+=2)
	{
		for(int j=0; j<NCells(i); ++j)
		{
			for(int n=0; n<GlxNeigh[i]->oN(j); ++n)
			{
				fmat res = GlxNeigh[i]->oA(j,n)*(GlxStars(i).slice(j)(span(8,max_age),span::all)
				 - GlxStars(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all));
				
				if (accu( abs(res) ) > prm.at("Minimum_Diff_Mass_f").f)
				{
					GlxStarsBuff(i).slice(j)(span(8,max_age),span::all)= GlxStarsBuff(i).slice(j)(span(8,max_age),span::all)-res;
					GlxStarsBuff(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all)=
					GlxStarsBuff(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all)+res;
				}
			}
		}
	}
	//work on odd rings
	#pragma omp parallel for schedule(guided)
	for(int i=1; i<prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i; i+=2)
	{
		for(int j=0; j<NCells(i); ++j)
		{
			for(int n=0; n<GlxNeigh[i]->oN(j); ++n)
			{
				fmat res = GlxNeigh[i]->oA(j,n)*(GlxStars(i).slice(j)(span(8,max_age),span::all)
				 - GlxStars(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all));
				
				if (accu( abs(res) ) > prm.at("Minimum_Diff_Mass_f").f)
				{
					GlxStarsBuff(i).slice(j)(span(8,max_age),span::all)= GlxStarsBuff(i).slice(j)(span(8,max_age),span::all)-res;
					GlxStarsBuff(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all)=
					GlxStarsBuff(GlxNeigh[i]->oR(j,n)).slice(GlxNeigh[i]->oC(j,n))(span(8,max_age),span::all)+res;
				}
			}
		}
	}
	// write results
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i; ++i)
	{
		GlxStars(i)+=prm.at("STELLAR_DIFUSSION_f").f*GlxStarsBuff(i);
	}
}

/*
 *  The funtion for reading rotation curve data from file.
 *  
 *  First lines starting with # are skipped (can be more than one)
 *  First columns is skipped, assumed that user set rotation curve for every ring of the galaxy correctly
 *  Second columns is expected to be in km/s units
 */
void galaxy::read_rotation_curve()
{
	fstream rotation_curve;
	rotation_curve.open(prm.at("Rotation_curve_s").str.c_str(),ios::in);
	if(!rotation_curve.is_open())
	{
		cerr<<"Rotation curve named: "<<prm.at("Rotation_curve_s").str << " was not found." << endl;
		exit(EXIT_FAILURE);
	}
	strip_header(rotation_curve);

	fmat raw_rotation_curve;
	raw_rotation_curve.load(rotation_curve, raw_ascii);
	GlxRings.col(gr.at("RotRad")) = raw_rotation_curve.col(1);
	
	//convert km/s to rad/time step
	GlxRings.col(gr.at("RotRad")) = 2*datum::pi * KMS_PER_MYR_TO_PC_PER_MYR * prm.at("Time_step_Myr_i").i * GlxRings.col(gr.at("RotRad"))/(2*datum::pi*GlxRings.col(gr.at("RadPc")));
	GlxRings(0,gr.at("RotRad")) = GlxRings(1,gr.at("RotRad"));
	rotation_curve.close();
	raw_rotation_curve.clear();
}

/*
 *  The funtion updates cells neighbors information (called after each rotation)
 *  
 */
void galaxy::find_neighbours()
{
	#pragma omp parallel for schedule(guided)
	for(int ring_id=1; ring_id<prm.at("Rings_i").i; ++ring_id)
	{
		double upper_diff, lower_diff;
		lower_diff = fmod( (GlxRings(ring_id, gr.at("AziCor")) - GlxNeigh[ring_id]->HCellW)-
			(GlxRings(ring_id - 1, gr.at("AziCor")) - GlxNeigh[ring_id - 1]->HCellW),
			2*datum::pi)/GlxNeigh[ring_id - 1]->CellW + NCells(ring_id-1);

		upper_diff = 0;

		if( ring_id + 1 < prm.at("Rings_i").i )
		{
			upper_diff = fmod((GlxRings(ring_id, gr.at("AziCor")) - GlxNeigh[ring_id]->HCellW)-
				(GlxRings(ring_id + 1, gr.at("AziCor")) - GlxNeigh[ring_id+1]->HCellW),
				2*datum::pi)/GlxNeigh[ring_id+1]->CellW + NCells(ring_id+1) + NCells(ring_id+1);
		}
		GlxNeigh[ring_id]->find_neighbours(lower_diff, upper_diff);
	}
}

/*
 *  The funtion creates triggered star formation events in the cells
 *  
 */

void galaxy::IzotropicTriggeredStarFormation()
{
	unsigned int SUM=0;
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		for(int j=0; j<NCells(i); ++j)
		{
			if(GlxCells(i)(j,gc.at("TrigTime")) > -0.5 && GlxCells(i)(j,gc.at("TrigTime")) < 0.5)
			{
				for(int n=0; n<GlxNeigh[i]->Num(j); ++n)
				{
					if(GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("SFR_TRESH")) > 0 )
					{
						double RND;
						float MStars;
						RND=gsl_rng_uniform(random_number_gsl);
						if( prm.at("TRIG_PROB_f").initialized ) {
							if( RND < prm.at("TRIG_PROB_f").f*GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("SFR_TRESH"))) {
								++SUM;
								MStars=SF_event(GlxNeigh[i]->Rings(j,n),GlxNeigh[i]->Cells(j,n));
								GlxRings(i,gr.at("SF_event"))+=1;
								GLX.at("SFR")+=MStars;
							}
						}
						else {
							if( RND < prm.at("TRIG_PROB_f").fv[TIMESTEP]*GlxCells(GlxNeigh[i]->Rings(j,n))(GlxNeigh[i]->Cells(j,n),gc.at("SFR_TRESH"))) {
								++SUM;
								MStars=SF_event(GlxNeigh[i]->Rings(j,n),GlxNeigh[i]->Cells(j,n));
								GlxRings(i,gr.at("SF_event"))+=1;
								GLX.at("SFR")+=MStars;
							}
						}
					}
				}
			}
		}
	}
	GLX.at("TRIGGERED_EVENTS")=SUM;
}

/*
 *  The funtion for the creation of spontaneous star formation events in the cells
 *  
 */
 
void galaxy::SpontaneousGas2StarFormation()
{
	double SUM=0;
	#pragma omp parallel for schedule(static) reduction(+:SUM)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxCells(i).col(gc.at("SPONT_BUFFER"))=GlxCells(i).col(gc.at("Mgas"))%GlxCells(i).col(gc.at("Mgas"));
		SUM+=sum( GlxCells(i).col( gc.at("SPONT_BUFFER") ) );
	}
	SUM=prm.at("GAS2_SPONT_NRM_f").f/SUM;

	unsigned int SUM0=0;
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		double RND;
		float MStars;
		for(unsigned short j=0; j<NCells(i); ++j)
		{
			if(GlxCells(i)(j,gc.at("SFR_TRESH")) > 0)
			{

				RND=gsl_rng_uniform (random_number_gsl);

				if(RND < GlxCells(i)(j,gc.at("SPONT_BUFFER"))*SUM)
				{
					MStars=SF_event(i,j);
					++SUM0;
					GlxRings(i,gr.at("SF_event"))+=1;
					GlxRings(i,gr.at("SPONT_event"))+=1;
					
					GLX.at("SFR")+=MStars;
				}
			}
		}
	}
	GLX.at("SPONTANEOUS_EVENTS")=SUM0;
}

/*
 *  The funtion is called from the evolve loop. It is the place for
 *  the implementation of custom spontaneous SF funtion
 *  
 */

void galaxy::SpontaneousGeneric()
{
	switch(SPONT_TYPE)
	{
		/* Custom types of spontaneous star formation can be added hereOther */
		case 0:
			SpontaneousGas2StarFormation();
			break;
		default:
			cerr<<"Unsupported spontaneous type:"<< SPONT_TYPE <<", [0-5]"<<endl;
			exit(EXIT_FAILURE);
	};
}

/*
 *  The funtion is called for every SF event in the cell, it creates SSP according to the SFE law defined in the parameters file
 *  ON/OFF switch of the metal enriched outflow is also calculated here using description of
 *  Baumgartner V., Breitschwerdt D. 2013, A&A, 557, A140
 *  
 */
double galaxy::SF_event(unsigned short i, unsigned short j)
{
	//const double  SQRT_INVERSE_AREA_POW_1_5 = pow(prm.at("CELL_AREA_f").f,-1.5);
	//const double INVERSE_AREA = 1./prm.at("CELL_AREA_f").f;
	const double NORMALIZE_SFE_MASS = prm.at("NRM_SFE_f").f;

	if(GlxCells(i)(j,gc.at("Mgas")) < 0 || GlxCells(i)(j,gc.at("Metals")) < 0)
	{
		cerr<<i<<" "<<j<<" "<<GlxCells(i)(j,gc.at("Mgas"))<< " "<< GlxCells(i)(j,gc.at("Zgas"))<<endl;
	}
	double SFE=-99999;//, mstr, MOL_RATIO, Pressure;
	
	// Check if SFR controlling variables are set as time variables or constants
	if( prm.at("SFE_f").initialized && prm.at("SFE_POW_f").initialized) {
		SFE = prm.at("SFE_f").f*pow(GlxCells(i)(j,gc.at("Mgas"))/NORMALIZE_SFE_MASS, prm.at("SFE_POW_f").f);
	}
	else if( prm.at("SFE_f").initialized && !prm.at("SFE_POW_f").initialized) {
		SFE = prm.at("SFE_f").f*pow(GlxCells(i)(j,gc.at("Mgas"))/NORMALIZE_SFE_MASS, prm.at("SFE_POW_f").fv[TIMESTEP]);
	}
	else if( !prm.at("SFE_f").initialized && prm.at("SFE_POW_f").initialized) {
		SFE = prm.at("SFE_f").fv[TIMESTEP]*pow(GlxCells(i)(j,gc.at("Mgas"))/NORMALIZE_SFE_MASS, prm.at("SFE_POW_f").f);
	}
	else if( !prm.at("SFE_f").initialized && !prm.at("SFE_POW_f").initialized) {
		SFE = prm.at("SFE_f").fv[TIMESTEP]*pow(GlxCells(i)(j,gc.at("Mgas"))/NORMALIZE_SFE_MASS, prm.at("SFE_POW_f").fv[TIMESTEP]);
	}
	//------------------------------------
	if (SFE < prm.at("MINIMUM_SFE_f").f)
		SFE = prm.at("MINIMUM_SFE_f").f;

	if(SFE > prm.at("MAXIMUM_SFE_f").f)
		SFE = prm.at("MAXIMUM_SFE_f").f;

	double Mstars = GlxCells(i)(j,gc.at("Mgas"))*SFE;

	// if mass of formed stars is bellow limit for ~1 supernovue, do not initiate triggering
	if(Mstars < prm.at("MINIMUM_MASS_f").f)
	{
		return 0;
	}
	else if (Mstars > prm.at("TRIGG_MASS_f").f)
	{
		GlxCells(i)(j,gc.at("RefTime"))=-1;
		GlxCells(i)(j,gc.at("TrigTime"))=-prm.at("TRIG_TIME_i").i;
	}

	GlxCells(i)(j,gc.at("SFE"))=SFE;
	GlxCells(i)(j,gc.at("last_Mstr"))=Mstars;
	GlxRings(i,gr.at("SFR"))+=Mstars;
	RSFH(i,SFR_INDEX)+=Mstars;
	
	//double L_0 = Mstars*1e44/0.29*0.0014*(5./11)/40e6/3.156e7;
	double N_SN = Mstars/0.29*0.0014;
	// 0.29   -- average mass of stars
	// 0.0014 -- fraction fo SN
	// 2/11   -- fraction of energy converted to kinetic energy
	// 40e6   -- period of SN activity in years
	// 3.156e7-- number of seconds per year
	// 1e44 energy of SN, equvivavlent to 1e51 erg
	double Rho_0 = GlxCells(i)(j,gc.at("Mgas"))*2e30/prm.at("CELL_AREA_f").f/(2.*100.)/pow(3.085e18,3)/1.66e-27;
	GlxCells(i)(j,gc.at("rho")) = Rho_0;
	// 2e30    -- Sun mass
	// 10472   -- area of the cell
	// 2*100   -- 2x times scale height of the gas density
	// 3.085e18^-3 -- convert from kg/pc^3 -> kg/cm^3
	// 1.66e-27    -- atomical mass unit
	
	double H0 = 100.;
	//double H0 = 3.085e16*100;
	// scale height in meters
	
	const double delta=0.82;
	if (prm.at("OUTFLOW_VEL_RATIO_f").f > 0)
	{
		if (N_SN <= 1)
		{
			GlxCells(i)(j,gc.at("TVEL"))=-9.99;
		}
		else if (N_SN >= 500)
		{
			GlxCells(i)(j,gc.at("TVEL"))=9999;
		}
		else
		{
			double Vacc=prm.at("OUTFLOW_VEL_RATIO_f").f;
			
			double D0=1854.5*(Rho_0/0.5)*pow(H0/500.,2.-delta)*pow(Vacc,delta+3);
			double N_SN_REQ = 0.011*pow(D0,1.041)-3.772;
			if (N_SN > N_SN_REQ)
			{
				GlxCells(i)(j,gc.at("TVEL"))=2;
			}
			else
			{
				GlxCells(i)(j,gc.at("TVEL"))=-1;
			}
			GlxCells(i)(j,gc.at("OflN0"))=D0;
		}
	}
	else if (prm.at("OUTFLOW_VEL_RATIO_f").f == 0)
	{
		GlxCells(i)(j,gc.at("TVEL"))=2.;
	}
	else if (prm.at("OUTFLOW_VEL_RATIO_f").f < 0)
	{
		GlxCells(i)(j,gc.at("TVEL"))=-99.99;
	}

	//GlxCells(i)(j,gc.at("TVEL"))=pow(L_0/(Rho_0*H0*H0),1./3.)*0.7*0.001; // km/s
	float Zg_before=GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas"));
	unsigned short zidx = assign_stellar_metallicity(GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas")), GlxCells(i)(j,gc.at("Mgas")), i, j);
	GlxStars(i)(0,zidx,j)+=Mstars;
	GlxCells(i)(j,gc.at("Mgas"))-=Mstars;
	GlxCells(i)(j,gc.at("Metals"))-=Mstars*prm.at("SSP_metalls_fv").fv(zidx);
	
	#pragma omp atomic
	TOTAL_METALS-=Mstars*prm.at("SSP_metalls_fv").fv(zidx);
	
	if(GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas"))<0.0001)
	{
		cerr<<"What? SF_event: "<<GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas"));
		cerr<<" Zs: "<< prm.at("SSP_metalls_fv").fv(zidx)<<" Za"<< Zg_before <<" TIMESTEP: "<<TIMESTEP<<endl;
	}
	return Mstars;
}

/*
 *  The funtion is called every time step in the evolve loop.
 *  It marks cells which are currently capable to form stars
 *  
 */

void galaxy::MarkActiveCells()
{
	static const float  INVERSE_GAS_TRESH = 1./(prm.at("GAS_SFR_TRESHOLD_f").f*prm.at("CELL_AREA_f").f);
	static const float COMPARE_VAL = prm.at("GAS_SFR_TRESHOLD_f").f*prm.at("CELL_AREA_f").f;

	if(!prm["SHARP_SFTRESHOLD"].initialized)
	{
		#pragma omp parallel for schedule(guided)
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			for(int j=0; j<NCells(i); ++j)
			{
				if( GlxCells(i)(j,gc.at("RefTime"))>prm.at("SN_Timescale_10Myr_i").i)
				{
					if ( GlxCells(i)(j,gc.at("Mgas")) >=  COMPARE_VAL )
						GlxCells(i)(j, gc.at("SFR_TRESH")) = 1;
					else
						GlxCells(i)(j, gc.at("SFR_TRESH")) = GlxCells(i)(j,gc.at("Mgas"))*INVERSE_GAS_TRESH;
				}
				else
					GlxCells(i)(j, gc.at("SFR_TRESH"))=0;
			}
			if( prm["COUNT_ACTIVE"].initialized )
			{
				GlxRings(i,gr.at("ACTIVE"))+=sum(GlxCells(i).col(gc.at("SFR_TRESH")));
			}
		}
	}
	else
	{
		#pragma omp parallel for schedule(guided)
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			for(int j=0; j<NCells(i); ++j)
			{
				if(GlxCells(i)(j,gc.at("RefTime"))>prm.at("SN_Timescale_10Myr_i").i    &&
				  GlxCells(i)(j,gc.at("Mgas")) >=  COMPARE_VAL 									)
				{
					GlxCells(i)(j, gc.at("SFR_TRESH"))=1;
				}
				else
				{
					GlxCells(i)(j, gc.at("SFR_TRESH"))=0;
				}
			}
			if( prm["COUNT_ACTIVE"].initialized )
			{
				GlxRings(i,gr.at("ACTIVE"))+=sum(GlxCells(i).col(gc.at("SFR_TRESH")));
			}
		}
	}
}

/*
 *  The funtion time integrates galaxy parameters
 *  
 */
void galaxy::Evolve()
{
	SFR_INDEX=0; // sFR index is used to get longer then 10 Myr SFR profile
	OTFL_INDEX=0;
	TIMESTEP=0;
	for(int t=0; t<prm.at("Galaxy_age_Myr_i").i; t+=prm.at("Time_step_Myr_i").i)
	{		
		// rotation of the cells
		GlxRings.col(gr.at("AziCor"))=GlxRings.col(gr.at("AziCor"))+GlxRings.col(gr.at("RotRad"));
		// after the rotation, recalculate neighbors of the cells
		find_neighbours();
		
		if (SFR_INDEX >= MAX_SFR_INDEX)
		{
			SFR_INDEX=0;
		}
		if (OTFL_INDEX >= MAX_OTFL_INDEX)
		{
			OTFL_INDEX=0;
		}
		RSFH.col(SFR_INDEX).fill(0);
		ROTFL_GAS.col(OTFL_INDEX).fill(0);
		ROTFL_METALS.col(OTFL_INDEX).fill(0);

		Accretion();

		//! evolve metals before redistributing
		EvolveMetals();
		if(prm.at("Diffusion_const_f").f >= 0)
		{
			gas_diffusion();
		}

		if(prm.at("STELLAR_DIFUSSION_f").f > 0)
		{
			stars_diffusion();
		}
		EvolveAge();

		MarkActiveCells();

		SpontaneousGeneric();

		IzotropicTriggeredStarFormation();

		++TIMESTEP;
		#pragma omp parallel for
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			GlxCells(i).col(gc.at("RefTime"))=GlxCells(i).col(gc.at("RefTime"))+1;
			GlxCells(i).col(gc.at("TrigTime"))=GlxCells(i).col(gc.at("TrigTime"))+1;
		}

		Output();
		ResetVariables();
		
		SFR_INDEX++;
		OTFL_INDEX++;	
	}
	check_metals();
}

/*
 *  The funtion called ant the end of simulation prints diagnostic
 *  information about metallicity content in galaxy
 *  
 */
void galaxy::check_metals()
{
	double CUR_METALS=0;
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		for(int j=0; j<NCells(i); ++j)
		{
			CUR_METALS+=GlxCells(i)(j,gc.at("Metals"));
		}
	}
	cerr<<"CUR_METALS:"<<CUR_METALS<< " " <<" TOTAL_METALS: "<<TOTAL_METALS << " "<<CUR_METALS/TOTAL_METALS<<endl;
}

/*
 *  The funtion process the model's ouputs
 *  
 */
void galaxy::Output()
{
	if(TIMESTEP*prm.at("Time_step_Myr_i").i == prm.at("Output_Times_ia").usv(prm.at("Output_Times_ia").cur) || prm.at("OUTPUT_TYPES_n").dict.at("pegase") )
	{
		if(prm.at("OUTPUT_TYPES_n").dict.at("pegase"))
		{
			GenRadialPhotometry();
		}
	}
	if(TIMESTEP*prm.at("Time_step_Myr_i").i == prm.at("Output_Times_ia").usv(prm.at("Output_Times_ia").cur) && prm.at("OUTPUT_TYPES_n").dict.at("2d"))
	{
		stringstream s, sdump;
		s<<setprecision(3);
		s<<"#r-kpc x-kpc y-kpc ref_t-tstep sfr_t sp_buff mgas-msol-pc2 zgas mstr-msol-pc2 sfe-ratio last_mstr-msol";
		if( prm.at("OUTFLOW_b").initialized )
		{
			s<<" TVEL rho N0" ;
		}
		if( prm.at("PHOTOMETRY_sc").initialized )
		{
			for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
				s<<" "<< it->first << "-mag-arcsec2 ";
			for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
				s<<" ave_"<< it->first << " ";
		}
		s<<endl;
		#pragma omp parallel for  schedule(dynamic)
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			float mstr;
			stringstream ss;
			for(int j=0; j<NCells(i); ++j)
			{
				ss<<GlxRings(i,gr.at("RadPc"))/1000<<" ";
				ss<<GlxRings(i,gr.at("RadPc"))*cos(GlxRings(i,gr.at("AziCor"))+GlxRings(i,gr.at("CellW"))*j)/1000<<" ";
				ss<<GlxRings(i,gr.at("RadPc"))*sin(GlxRings(i,gr.at("AziCor"))+GlxRings(i,gr.at("CellW"))*j)/1000<<" ";
				ss<<GlxCells(i)(j,gc.at("RefTime"))<<" ";
				ss<<GlxCells(i)(j, gc.at("SFR_TRESH"))<<" ";
				ss<<GlxCells(i)(j,gc.at("SPONT_BUFFER"))<<" ";
				ss<<GlxCells(i)(j,gc.at("Mgas"))/prm.at("CELL_AREA_f").f<<" ";
				ss<<GlxCells(i)(j,gc.at("Metals"))/GlxCells(i)(j,gc.at("Mgas"))<<" ";
				mstr=accu(GlxStars(i).slice(j))/prm.at("CELL_AREA_f").f;
				ss<<mstr<<" ";
				ss<<GlxCells(i)(j,gc.at("SFE"))<<" ";
				ss<<GlxCells(i)(j,gc.at("last_Mstr"))<<" ";
				if( prm.at("OUTFLOW_b").initialized )
				{
					ss << GlxCells(i)(j,gc.at("TVEL")) << " ";
					ss << GlxCells(i)(j,gc.at("rho")) << " ";
					ss << GlxCells(i)(j,gc.at("OflN0")) << " ";
				}
				if( prm.at("PHOTOMETRY_sc").initialized )
				{
					if( mstr != 0 )
					{
						for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
						{
							if (it->first == string("LHa"))
							{
								if(GlxPhot( i )(j, gp.at( it->first )) > 0)
									ss<< log10(GlxPhot( i )(j, gp.at( it->first ))/(prm.at("CELL_AREA_f").f))<<" ";
								else
									ss<<" \"\" ";
								
							}
							else
								ss<< -2.5*log10(GlxPhot( i )(j, gp.at( it->first ))/(prm.at("CELL_AREA_f").f*PC2_TO_ARCSEC2))<<" ";
						}
						//average flux
						for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
						{
							if(i != 0 && is_finite( log10(AveFluxBlock( i, j, gp.at( it->first )))) )
							{
								ss<< -2.5*log10( AveFluxBlock( i, j, gp.at( it->first )) )<<" ";
							}
							else
							{
								ss<<" \"\" ";
							}
						}
					}
					else
					{
						for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
							ss<<" \"\" "<< " \"\" " ;
					}
				}
				ss<<endl;
			}
			#pragma omp critical
			{
				s<<ss.str();
			}
		}
		fstream f;
		stringstream fs, fdump;
		fs << iname<<"_cells_"<<TIMESTEP*prm.at("Time_step_Myr_i").i<<".dat";
		f.open(fs.str().c_str(),ios::out);
		f << s.str();
		f.close();
	}
	if(TIMESTEP*prm.at("Time_step_Myr_i").i == prm.at("Output_Times_ia").usv(prm.at("Output_Times_ia").cur) && prm.at("OUTPUT_TYPES_n").dict.at("1d"))
	{
		stringstream s;
		float mstr;
		s<<setprecision(3);
		s<<"#r-kpc mgas-msol-pc2 zgas mstr-msol-pc2 SF_events-num SP_events-num SFR-msol-pc2-tstep SFR100-msol-pc2-10tsteps ACTIVE-num Tgas-msol-pc2";
		if( prm.at("OUTFLOW_b").initialized )
		{
			s<<" Ogas_tot Ogas_cur Ometals_tot Ometals_cur";
		}
		if( prm.at("PHOTOMETRY_sc").initialized )
			for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
				s<<" "<< it->first;
		s<<endl;
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			s<<GlxRings(i,gr.at("RadPc"))/1000<<" ";
			s<<mean(GlxCells(i).col(gc.at("Mgas")))/prm.at("CELL_AREA_f").f <<" ";
			s<<median(GlxCells(i).col(gc.at("Metals"))/GlxCells(i).col(gc.at("Mgas"))) <<" ";
			mstr = accu(GlxStars(i))/(prm.at("CELL_AREA_f").f*NCells(i));
			s<< mstr <<" ";
			s<<GlxRings(i,gr.at("SF_event") )<<" ";
			s<<GlxRings(i,gr.at("SPONT_event"))<<" ";
			s<<GlxRings(i,gr.at("SFR"))/(prm.at("CELL_AREA_f").f*NCells(i)*prm.at("Time_step_Myr_i").i)<<" ";
			s<<sum(RSFH.row(i))/(prm.at("CELL_AREA_f").f*NCells(i)*prm.at("Time_step_Myr_i").i*MAX_SFR_INDEX)<<" ";
			s<<GlxRings(i,gr.at("ACTIVE") )<<" ";
			s<<GlxRings(i,gr.at("TMgas"))/prm.at("CELL_AREA_f").f <<" ";
			if( prm["OUTFLOW_b"].initialized )
			{
				s << GlxRings(i,gr.at("O_mgas"))/NCells(i) << " ";
				s << sum(ROTFL_GAS.row(i))/(prm.at("CELL_AREA_f").f*NCells(i)*prm.at("Time_step_Myr_i").i*MAX_OTFL_INDEX) << " ";
				s << GlxRings(i,gr.at("O_metals"))/NCells(i) << " ";
				s << sum(ROTFL_METALS.row(i))/(prm.at("CELL_AREA_f").f*NCells(i)*prm.at("Time_step_Myr_i").i*MAX_OTFL_INDEX) << " ";
			}
			
			if( prm.at("PHOTOMETRY_sc").initialized )
			{
				if( mstr != 0 )
				{
					for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
					{
						if (it->first == string("LHa"))
						{
							s<< log10( accu( GlxPhot(i).col(gp.at( it->first )) )/(prm.at("CELL_AREA_f").f*NCells(i)) ) <<" ";
						}
						else if (it->first == string("M")  || (it->first == string("MWD")) || (it->first == string("MBHNS")))
						{
							s<< accu( GlxPhot(i).col(gp.at( it->first )) )/(prm.at("CELL_AREA_f").f*NCells(i))  <<" ";
						}
						else
							s<< -2.5*log10( accu( GlxPhot(i).col(gp.at( it->first )) )/(prm.at("CELL_AREA_f").f*NCells(i)*PC2_TO_ARCSEC2) ) <<" ";
					}
				}
				else
				{
					for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
						s<<" \"\" ";
				}
				s << endl;
			}
		}
		fstream f;
		stringstream fs;
		fs << iname<<"_rings_"<<TIMESTEP*prm.at("Time_step_Myr_i").i<<".dat";
		f.open(fs.str().c_str(),ios::out);
		f << s.str();
		f.close();

		//==============================================================
		// Reset some varibles
		//--------------------------------------------------------------
		GlxRings.col(gr.at("SPONT_event")).zeros();
	}
	if(TIMESTEP*prm.at("Time_step_Myr_i").i == prm.at("Output_Times_ia").usv(prm.at("Output_Times_ia").cur) && prm.at("OUTPUT_TYPES_n").dict.at("cmd"))
	{
		static const float CMD_AGE_LIMIT=prm.at("CMD_AGE_LIMIT_f").f;
		static const float CMD_MASS_LIMIT=prm.at("CMD_MASS_LIMIT_f").f;
		cerr<<"cmd"<<endl;
		fstream f;
		stringstream fs;
		fs << iname<<"_cmd_"<<TIMESTEP*prm.at("Time_step_Myr_i").i<<".dat";
		f.open(fs.str().c_str(),ios::out);
		stringstream s;

		int NAges = prm.at("SSP_ages_usv").usv.n_rows;
		int Nz =  prm.at("SSP_metalls_fv").fv.n_rows;
		cerr<<NAges <<" "<< Nz << endl;
		f<<"#r-kpc phi-radians logM-msol T-myr Z" <<endl;
		f<<setprecision(3);
		#pragma omp parallel for  schedule(dynamic)
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			stringstream ss;
			for(int j=0; j<NCells(i); ++j)
			{
				for(int a=0; a<NAges; ++a)
				{
					if(prm.at("SSP_ages_usv").usv(a) < CMD_AGE_LIMIT )
					for(int z=0; z<Nz; ++z)
					{
						if( GlxStars(i)(a,z,j)> CMD_MASS_LIMIT )
						{
							ss<< GlxRings(i,gr.at("RadPc"))/1000 <<" ";
							ss<< setprecision(5);
							ss<< fmodf(GlxRings(i,gr.at("AziCor")),(float)2*datum::pi)+GlxRings(i,gr.at("CellW"))*j << " ";
							ss<< log10(GlxStars(i)(a,z,j)) <<" ";
							ss<< prm.at("SSP_ages_usv").usv(a) <<" ";
							ss<< prm.at("SSP_metalls_fv").fv(z) <<endl;
						}
					}
				}
			}
			#pragma omp critical
			{
				s<<ss.str();
			}
		}
		f << s.str();
		f.close();
	}
	if(TIMESTEP*prm.at("Time_step_Myr_i").i == prm.at("Output_Times_ia").usv(prm.at("Output_Times_ia").cur))
	{
		if(prm.at("Output_Times_ia").cur + 1 < prm.at("Output_Times_ia").an)
		{
			++prm.at("Output_Times_ia").cur;
		}
	}
	//output integrated galaxy properties
	if(prm.at("OUTPUT_TYPES_n").dict.at("0d"))
	{
		if(TIMESTEP == 1)
		{
			igal <<"#t-myr SP_E-num TR_E-num ACC-msol-yr ST_GAS_ACC-msol-yr OTFL-msol-yr STARS-msol GAS-msol ZGAS TSFR-msol-yr"<<endl;
		}
		igal<< TIMESTEP*prm.at("Time_step_Myr_i").i <<" ";
		igal<< GLX.at("SPONTANEOUS_EVENTS")<<" ";
		igal<< GLX.at("TRIGGERED_EVENTS")  <<" ";
		igal<< GLX.at("ACCRETION")/(prm.at("Time_step_Myr_i").i*1e6)<<" ";
		igal<< GLOB_STELLAR_ACC/(prm.at("Time_step_Myr_i").i*1e6)<<" ";
		igal<< GLOB_OUTFLOW_MGAS/(prm.at("Time_step_Myr_i").i*1e6)<<" ";
		
		double Galaxy_Mgas=0, Galaxy_Mstars=0;
		#pragma omp parallel for schedule(static) reduction(+:Galaxy_Mgas,Galaxy_Mstars)
		for(int i=0; i<prm.at("Rings_i").i; ++i)
		{
			Galaxy_Mgas+=accu(GlxCells(i).col(gc.at("Mgas")));
			Galaxy_Mstars+=accu(GlxStars(i));
		}
		igal<< Galaxy_Mstars<<" ";
		igal<< Galaxy_Mgas<<" ";
		igal<< TOTAL_METALS/Galaxy_Mgas<<" ";
		
		igal<<" "<<GLX.at("SFR")/(prm.at("Time_step_Myr_i").i*1e6)<<endl;
		for(map<string,float>::iterator it=GLX.begin(); it!= GLX.end(); ++it)
		{
			it->second=0;
		}
	}
	if(prm.at("OUTPUT_TYPES_n").dict.at("pegase"))
	{	
		if(TIMESTEP == 1)
		{
			ipeg <<"#t";
			for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
				for(int i=0; i<prm.at("PEG_RINGS_i").i;++i)
					ipeg<<" "<< it->first<<i;
			ipeg<<endl;
		}
		ipeg<< TIMESTEP*prm.at("Time_step_Myr_i").i;
		for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
			for(int i=0; i<prm.at("PEG_RINGS_i").i;++i)
			{
				if( GlxRadPhot(i,gp.at( it->first )) > 0 )
					ipeg<<" "<<-2.5*log10( GlxRadPhot(i,gp.at( it->first )) );
				else
					ipeg<<" "<<"\"\"";
			}
		ipeg<<endl;
	}
	if(prm.at("OUTPUT_TYPES_n").dict.at("ring"))
	{
		fvec mstr = fvec(NCells(prm.at("RING_i").i));
		fvec fuv = fvec(NCells(prm.at("RING_i").i));
		fvec nuv = fvec(NCells(prm.at("RING_i").i));
		if(TIMESTEP == 1)
		{
			iring <<"#t";
			for(int j=0; j<NCells(prm.at("RING_i").i);++j)
			{
				iring<<" "<< "mstr"<<j<<" "<<"fuv"<<j<<" "<<"nuv"<<j<<" reft"<<j;
			}
			iring<<endl;
		}
		iring<< TIMESTEP;

		#pragma omp parallel for
		for(int j=0; j<NCells(prm.at("RING_i").i);++j)
		{
			mstr(j) = accu(GlxStars(prm.at("RING_i").i).slice(j))/prm.at("CELL_AREA_f").f;
			if( mstr(j) > 0 )
			{
				fuv(j)=-2.5*log10( accu(GlxStars(prm.at("RING_i").i).slice(j) % PHOT_PEGASE.at("FUV") ) );
				nuv(j)=-2.5*log10( accu(GlxStars(prm.at("RING_i").i).slice(j) % PHOT_PEGASE.at("NUV") ) );
			}
		}
		for(int j=0; j<NCells(prm.at("RING_i").i);++j)
		{
			if( mstr(j) > 0 )
				iring<<" "<<mstr(j) << " "<< fuv(j) << " " << nuv(j) <<" "<<GlxCells(prm.at("RING_i").i)(j,gc.at("RefTime"));
			else
			{
				for(int foo=0; foo<3; ++foo)
					iring<<" "<<"\"\"";
				iring<<" "<<GlxCells(prm.at("RING_i").i)(j,gc.at("RefTime"));
			}
		}
		iring<<endl;
	}

}

/*
 *  The funtion sums flux from all ssp present in the cells
 *  to produce 2-D surface brightness image
 *  
 */

void galaxy::GenPhotometry()
{
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		for(unsigned short j=0; j<NCells(i); ++j)
		{
			for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
			{
				GlxPhot(i)(j,it->second) = accu(GlxStars(i).slice(j) % PHOT_PEGASE.at(it->first) );
			}
		}
	}
}

/*
 *  The funtion sums flux from all cells in the ring to produce 1-D 
 *  surface brightness profiles
 *  
 */

void galaxy::GenRadialPhotometry()
{
	GlxRadPhot.zeros();
	#pragma omp parallel for schedule(guided)
	for(int i=0; i<prm.at("PEG_RINGS_i").i; ++i)
	{	
		for(map<string,unsigned short>::iterator it=gp.begin(); it!= gp.end(); ++it)
			{
				GlxRadPhot(i,it->second) = sum( GlxPhot(i).col(it->second)  );
			}
	}
}

/*
 *  The funtion averages flux between nearby cells.
 *  Basically, the averaging occurs between nearby four cells
 *  
 */
 
float galaxy::AveFluxBlock(int i, int j, int filter)
{
	float CELL2 = prm.at("Cell_size_pc_f").f*prm.at("Cell_size_pc_f").f;
	float same_ring_area = ((i+0.5)*(i+0.5)-(i*i))*CELL2*datum::pi;
	float AveFlux = GlxPhot( i )(j, filter)*same_ring_area*
	GlxNeigh[i]->HCellW /(2*datum::pi)/prm.at("CELL_AREA_f").f;

	//neighbour in same ring
	AveFlux += GlxPhot( GlxNeigh[i]->oR(j,0) )(GlxNeigh[i]->oC(j,0), filter)*same_ring_area*
	GlxNeigh[i]->HCellW /(2*datum::pi)/prm.at("CELL_AREA_f").f;

	float AveFluxArea = 2*same_ring_area*GlxNeigh[i]->HCellW /(2*datum::pi);

	//upper ring this neighbour
	float AreaSum=0;
	float up_ring_area = ((i+1)*(i+1)-(i+0.5)*(i+0.5))*CELL2*datum::pi;
	for(int n=1; n<GlxNeigh[i]->oN(j); ++n)
	{
		AreaSum+=GlxNeigh[i]->oA(j,n)*GlxNeigh[i]->PTotal/GlxNeigh[i]->Pdown;
		if(AreaSum > 0.5)
		{
			AveFlux += GlxPhot( GlxNeigh[i]->oR(j,n) )(GlxNeigh[i]->oC(j,n), filter)*up_ring_area*(AreaSum-0.5)*GlxNeigh[ GlxNeigh[i]->oR(j,n) ]->CellW
			/(2*datum::pi)/prm.at("CELL_AREA_f").f;
			AveFluxArea += up_ring_area*(AreaSum-0.5)*GlxNeigh[ GlxNeigh[i]->oR(j,n) ]->CellW/(2*datum::pi);
		}
	}
	//upper ring next neighbour
	AreaSum=0;
	float OldAreaSum;
	int jj = GlxNeigh[i]->oC(j,0);
	for(int n=1; n<GlxNeigh[i]->oN( jj ); ++n)
	{
		OldAreaSum=AreaSum;
		AreaSum+=GlxNeigh[i]->oA( jj, n)*GlxNeigh[i]->PTotal/GlxNeigh[i]->Pdown;

		if(AreaSum < 0.5)
		{
			AveFlux += GlxPhot( GlxNeigh[i]->oR(jj,n) )(GlxNeigh[i]->oC(jj,n), filter)*up_ring_area*(AreaSum)*GlxNeigh[ GlxNeigh[i]->oR(jj, n) ]->CellW
			/(2*datum::pi)/prm.at("CELL_AREA_f").f;
			AveFluxArea += up_ring_area*(AreaSum)*GlxNeigh[ GlxNeigh[i]->oR(jj, n) ]->CellW/(2*datum::pi);
		}
		else
		{
			AveFlux += GlxPhot( GlxNeigh[i]->oR(jj,n) )(GlxNeigh[i]->oC(jj,n), filter)*up_ring_area*(0.5-OldAreaSum)*GlxNeigh[ GlxNeigh[i]->oR(jj, n) ]->CellW
			/(2*datum::pi)/prm.at("CELL_AREA_f").f;
			AveFluxArea += up_ring_area*(0.5-OldAreaSum)*GlxNeigh[ GlxNeigh[i]->oR(jj, n) ]->CellW/(2*datum::pi);
			break;
		}
	}
	return AveFlux/(AveFluxArea*PC2_TO_ARCSEC2);
}

/*
 *  The funtion assigns discrete metallicity for the formed ssp
 *  
 */
unsigned short galaxy::assign_stellar_metallicity(float z, float mgas, int ii, int jj)
{
	short Rval=-1;
	for(int i=0; i<prm.at("SSP_matrix_z_i").i-1; ++i)
	{
		if(z>=prm.at("SSP_metalls_fv").fv(0) && z<prm.at("SSP_metalls_fv").fv(1))
		{
			Rval=0;
			break;
		}
		else if(z>=prm.at("SSP_metalls_fv").fv(i) && z<prm.at("SSP_metalls_fv").fv(i+1))
		{
			//Rval=i;
			Rval=(z - prm.at("SSP_metalls_fv").fv(i))/prm.at("SSP_metalls_fv").fv(i) <
			(prm.at("SSP_metalls_fv").fv(i+1)-z)/prm.at("SSP_metalls_fv").fv(i+1) ?
			i : i+1;
			//chosing higher metallicity in rare cases results in lower than minimum metalicity
			//in the model
			//Rval=i;
			break;
		}
		else if(z>=prm.at("SSP_metalls_fv").fv(prm.at("SSP_matrix_z_i").i-1))
		{
			Rval=prm.at("SSP_matrix_z_i").i-1;
			break;
		}
	}
	if(Rval == -1 || Rval > prm.at("SSP_matrix_z_i").i-1)
	{
		cerr<<"Metallicity assigning error: t: "<< TIMESTEP<<" z: "<< z<<"  m: "
		<< mgas/prm.at("CELL_AREA_f").f<<" "<<prm.at("SSP_metalls_fv").fv(0)<<"  ";
		cerr<<prm.at("SSP_metalls_fv").fv(prm.at("SSP_matrix_z_i").i-1);
		cerr <<" | " << GlxCells(ii)(jj, gc.at("Diff_gas_rez")) << " " <<
		GlxCells(ii)(jj, gc.at("Diff_zgas_rez"))<< " " << GlxCells(ii)(jj, gc.at("Diff_zgas_rez")) / GlxCells(ii)(jj, gc.at("Diff_gas_rez"))<<endl;
		Rval=0;
	}
	return Rval;
}


/*
 *  The funtion reads accretion data files:
 *  - radial accretion efficiency;
 *  - time variable accretion into galaxy
 * 
 *  Here thearray with accretion values is created and later used
 *  in the Accretion() function each time step
 *  
 */
 
void galaxy::InitInfallScenario()
{
	int NSTEPS = prm.at("Galaxy_age_Myr_i").i/prm.at("Time_step_Myr_i").i;
	LCDM_ACC = fvec(NSTEPS,fill::zeros);
	
	fmat ACC_RAD;
	fmat ACC_TIME;
	//if accretion file is present, ignore analytic accretion efficiency
	if(prm["ACCRETION_b"].initialized)
	{
		fstream fACC_RAD;

		fACC_RAD.open(prm.at("ACCRETION_RAD_PROFILE_s").str.c_str());
		if(!fACC_RAD.is_open())
		{
			cerr<< "File :" << prm.at("ACCRETION_RAD_PROFILE_s").str.c_str() << " not found." << endl;
			exit(EXIT_FAILURE);
		}
		strip_header(fACC_RAD);
		ACC_RAD.load(fACC_RAD);
		GlxRings(span(0,prm.at("Rings_i").i-1), gr.at("AccEf"))=ACC_RAD.col(1);
		
		
		fstream fACC_TIME;
		
		fACC_TIME.open(prm.at("ACCRETION_TIME_PROFILE_s").str.c_str());
		if(!fACC_TIME.is_open())
		{
			cerr<< "File :" << prm.at("ACCRETION_TIME_PROFILE_s").str.c_str() << " not found." << endl;
			exit(EXIT_FAILURE);
		}
		strip_header(fACC_TIME);
		ACC_TIME.load(fACC_TIME);
				
		LCDM_ACC = ACC_TIME.col(1);
	}
	else
	{
		cerr<< "Accretion is not defined!" << endl;
		exit(EXIT_FAILURE);
	}
	
	cout.precision(4);
	cerr << "Galaxy mass : " << accu(LCDM_ACC)/1e9 << " 10^9 Msol" <<endl;
	
	//it should be normalized to the sum of unity
	GlxRings.col(gr.at("AccEf"))=GlxRings.col(gr.at("AccEf"))/
		accu(NCells(span(0,(prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i-1))) % GlxRings(span(0,(prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i-1)),gr.at("AccEf")));

	cerr << "After RADIAL_ACC_PROFILE normalization central and last value: " << GlxRings(0,gr.at("AccEf")) << " "
	 << GlxRings((prm.at("Rings_i").i-prm.at("Rings_BUFF_i").i-1),gr.at("AccEf")) <<endl;

}

/*
 *  The funtion reads every time step is called from the Evolve loop.
 *  It adds accreted gas mass to the cells in the galaxy disk
 *  
 */
void galaxy::Accretion()
{		
	double SUM=0;
	#pragma omp parallel for schedule(guided) reduction(+:SUM)
	for(int i=0; i<prm.at("Rings_i").i; ++i)
	{
		GlxRings(i,gr.at("TMgas"))+=LCDM_ACC(TIMESTEP)*GlxRings(i,gr.at("AccEf"));
		GlxCells(i).col(gc.at("Mgas"))+=LCDM_ACC(TIMESTEP)*GlxRings(i,gr.at("AccEf"));
		GlxCells(i).col(gc.at("Metals"))+=LCDM_ACC(TIMESTEP)*GlxRings(i,gr.at("AccEf"))*ACCRETION_METTALICITY;
		GlxCells(i).col(gc.at("Zgas"))=GlxCells(i).col(gc.at("Metals"))/GlxCells(i).col(gc.at("Mgas"));
		
		SUM+=LCDM_ACC(TIMESTEP)*GlxRings(i,gr.at("AccEf"))*NCells(i);
	}
	GLX.at("ACCRETION")=SUM;
	TOTAL_METALS+=LCDM_ACC(TIMESTEP)*ACCRETION_METTALICITY;
}


/*
 *  The function is called verey time to reset variables in the Evolve loop
 *  
 */
void galaxy::ResetVariables()
{
	// Put here variables which need to be setted to zero each timestep
	GlxRings.col(gr.at("SFR")).zeros(); // storage of SFR intensity in the galaxy ring
}

