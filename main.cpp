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
	galaxy* glx = new galaxy(std::string(argv[1]));
	cerr<<"Evolve"<<endl;
	timer.tic();
	glx->Evolve();

	cerr << "took " << timer.toc() << " seconds" << endl;

	//glx->Output();
	delete glx;
	return 0;
}
