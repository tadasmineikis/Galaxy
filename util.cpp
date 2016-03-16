#include <armadillo>
#include <iostream>

using namespace std;

/*
 * Linear 1D interpolation of monothonically increasing values of x
 */
arma::fvec interp1d(arma::fvec data_x, arma::fvec data_y, arma::fvec xx)
{
	/*;
	data_x.save("data_x.dat", arma::raw_ascii);
	data_y.save("data_y.dat", arma::raw_ascii);
	xx.save("xx.dat",arma::raw_ascii);*/

	arma::fvec yy= arma::fvec(xx.n_rows);
	arma::fvec diff0=data_x.rows(1, data_x.n_rows-1)-data_x.rows(0,data_x.n_rows-2);
	if ( arma::any(diff0<=0) )
	{
		cerr<<"data_x values are not monothonically increasing!" << endl;
		return yy.fill(arma::datum::nan);
	}
	diff0.clear();

	arma::fvec diff1 = xx.rows(1, xx.n_rows-1)-xx.rows(0,xx.n_rows-2);
	if ( arma::any(diff1)<=0)
	{
		cerr<<"xx values are not monothonically increasing!" << endl;
		return yy.fill(arma::datum::nan);
	}
	diff1.clear();

	float diff=-1e9;
	for(unsigned int i=0, j=1; i<xx.n_rows; ++i)
	{
		if(xx(i)<data_x(0) || xx(i)>data_x(data_x.n_rows-1))
		{
			//if xx out of range, return nan
				yy(i)=arma::datum::nan;
		}
		else
			for(; j<data_x.n_rows; )
			{
				if(data_x(j-1) < xx(i) && xx(i) <= data_x(j))
				{
					// if xx between data, interpolate linearly
					diff = (xx(i)-data_x(j-1)) / ( data_x(j)-data_x(j-1) );
					yy(i)=data_y(j-1)*(1.-diff)+data_y(j)*diff;
					break;
				}
				else if (xx(i) > data_x(j))
				{
					// if xx falls out of interval, jump to next data interval
					++j;
				}
			}
	}
	return yy;
}
/*
 * Strip header of input file stream. Header is identified by symbol -> # <- at the begining of each line
*/
void strip_header(fstream& in_stream)
{
	char peek;
	in_stream >> peek;
	//35 is equvivalent to symbol #
	while(peek == 35)
	{
		//throws away input of file length or util finds character \n
		in_stream.ignore(numeric_limits<streamsize>::max(), '\n');
		// check if header continues on next line
		in_stream >> peek;
	}
	// if header ended, putback peek'ed value
	in_stream.putback(peek);
}
