#include <string>
#include <armadillo>
#include <map>
#include <gsl/gsl_rng.h>

using namespace arma;

struct parameter
{
	parameter(): initialized( false )
	{}
    int i;
    std::string str;
    float f;
    int* ai;
    field<std::string> as;
    Col<unsigned short> usv;
    fvec fv;
    short int an;
    short int cur;
    std::string type;
    int n;
    std::map<std::string,bool> dict;
    std::map<unsigned short, unsigned short> idict;
    bool initialized;
};

class galaxy
{
	public:
	galaxy(std::string pfile);
	~galaxy();


	//void AnizotropicTriggeredStarFormation();
	void IzotropicTriggeredStarFormation();
	//void SpontaneousStarFormation();
	void SpontaneousGas2StarFormation();
	//void SpontaneousManualStarFormation();
	//void SpontaneousNotFixedStarFormation();
	//void SpontaneousFlatStarFormation();
	//void SpontaneousRadFlatStarFormation();
	void SpontaneousGeneric();
	void MarkActiveCells();
	void Evolve();
	void EvolveMetals();
	void EvolveAge();
	void ResetVariables();
	void Output();
	void Debug();
	void Debug1();
	void check_metals();
	void Read_params(std::string fname);
	void flush_all_dict(std::string fname);
	void violent_gas_exchange(int start_idx, int incr_idx);
	void calm_gas_exchange(int start_idx, int incr_idx);


	private:
	void read_rotation_curve();
	void find_neighbours();
	void rotate();
	void Fill();
	void gas_diffusion();
	void stars_diffusion();
	unsigned short assign_stellar_metallicity(float z, float mgas, int i, int j);
	double SF_event(unsigned short i, unsigned short j);
	void RefractoryTime();
	void InitInfallScenario();
	void Accretion();
	void GenPhotometry();
	float AveFluxBlock(int i, int j, int filter);
	void GenRadialPhotometry();

	unsigned int TIMESTEP;
	unsigned int SFR_INDEX, MAX_SFR_INDEX, OTFL_INDEX,MAX_OTFL_INDEX;
	fmat GlxRings;
	Col<unsigned short> NCells; //number of cells in ring
	field<mat> GlxCells;		//variable holding rings of galaxy as matrices varying in size
	field<fcube> GlxStars;		//variable holding SSP popualtions of galaxy
	field<fcube> GlxStarsBuff;	//variable holding SSP popualtions of galaxy
	fmat StarsFlowSpeed;		//rate of stellar aging in non-linear SSP grid
	//mat AgeFlow;				//temporary matrix for SSP ageing
	class neighbours** GlxNeigh;
	fmat MGAS_PEGASE;			// python produced pegase matrix
	fmat ZGAS_PEGASE;			// python produced pegase matrix
	fmat RSFH;                  // to store radial SFH for later producing averaged SFR profiles
	fmat ROTFL_GAS;             // to store radial outflow of gas
	fmat ROTFL_METALS;          // to store radial outflow of metals
	//field< Cube<unsigned short> > UpGlxNeigh;   //coordinates of neighbours
	//field< Cube<unsigned short> > LoGlxNeigh;   //coordinates of neighbours
	std::map<std::string,parameter> prm; // model parameters
	std::map<std::string,unsigned short> gr; // map of GlxRings columns and parameter name
	std::map<std::string,unsigned short> gc; // map of GlxRings columns and parameter name

	fvec LCDM_ACC;				// accretion variables

	gsl_rng*	random_number_gsl;			// gsl random number generator variables
	const gsl_rng_type * gsl_type;

	std::string iname;				    //parameter file name
	std::map<std::string,float> GLX;	//integrated galaxy parameters
	std::map<std::string,float> GLX_NRM;//normalisation of galaxy parameters
	std::map<int,std::string> i2s;
	std::map<int,std::string> i2m;
	std::map<int,std::string> i2z;
	std::map<int,std::string> i2l;

	std::map<std::string,fmat> PHOT_PEGASE;
	field<fmat> GlxPhot;
	fmat GlxRadPhot;
	fmat ESCAPE_VEL;
	std::map<std::string,unsigned short> gp;

	std::fstream igal, ipeg, iring;

	fmat MSPONT;
	unsigned short SPONT_TYPE;
	
	double PC2_TO_ARCSEC2;
	double GLOB_OUTFLOW_MGAS,GLOB_OUTFLOW_METALS;
	
	//for debugging
	double TOTAL_METALS;
	double GLOB_STELLAR_ACC;
};


