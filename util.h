/*
 * Linear 1D interpolation of monothonically increasing values of x
 */
arma::fvec interp1d(arma::fvec data_x, arma::fvec data_y, arma::fvec xx);

/*
 * Strip header of input file stream. Header is identified by symbol -> # <- at the begining of each line
*/
void strip_header(std::fstream& in_stream);
