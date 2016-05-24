#include <string>
#include "galaxy.h"
#include <armadillo>

int main(int argc, char* argv[])
{
	arma::wall_clock timer;
	if(argc < 2)
	{
		cerr<<"missing input parameter file and/or num of threads" << endl;
		exit(EXIT_FAILURE);
	}
	// creation of galaxy class; during startup it executes:
	// - allocation of memory
	// - initialization of a random number generator
	// - reading parameter file
	// - reading rotation curve file
	// - reading gas accretion file
	// - opening output files of type 0d, pegase and ring
	galaxy* glx = new galaxy(std::string(argv[1]));
	cerr<<"Evolve"<<endl;
	timer.tic();
	
	// Integration of the galaxy parameters takes place here
	glx->Evolve();

	// prints time spent on calculations
	cerr << "took " << timer.toc() << " seconds" << endl;

	// Galaxy class allocated memory clean-up
	delete glx;
	return 0;
}
