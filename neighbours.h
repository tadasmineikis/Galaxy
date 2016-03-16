#include <armadillo>

class neighbours
{
	public:
	neighbours(int ring_id);
	~neighbours();

	void find_neighbours(double lower_diff, double upper_diff);
	void init_Up_Lo_adr(neighbours* rUpper, neighbours* rLower);

	unsigned int Ring_id;
	double CellW;					//Cell width in radian units
	double HCellW;					//Half cell width in radian units
	double PTotal;					//total perimeter in pc
	double Pup, Pdown, Pside;		// upper, lower and side length of cell in relative units (side=1)
	double LowCorr;  				// correct lower mass flows
	unsigned int NCells;			// number of cell in the ring
	Mat<unsigned short> Cells, oC;		//coordinates (cell id)
	Mat<unsigned short> Rings, oR;		//coordinates (ring id)
	Col<unsigned short> Num, oN;		//ocuppied places in matrix
	mat Areas, oA;						//Areas of interaction
	class neighbours* Upper;		//reference to the upper ring
	class neighbours* Lower;		//reference to the lower ring
};
