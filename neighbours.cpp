#include "galaxy.h"
#include "neighbours.h"
#include <float.h>

neighbours::neighbours(int ring_id)
{
	Ring_id = ring_id;
	NCells = (ring_id == 0)?
	1 : ring_id*6;

	Cells = Mat <unsigned short>(NCells, 7);
	Rings = Mat <unsigned short>(NCells, 7);
	Num = Col <unsigned short>(NCells);
	Areas = mat(NCells, 7);

	Cells.fill(-1);
	Rings.fill(-1);
	Areas.zeros();

	oC = Mat <unsigned short>(NCells, 6);
	oR = Mat <unsigned short>(NCells, 6);
	oN = Col <unsigned short>(NCells);
	oA = mat(NCells, 6);

	oC.fill(-1);
	oR.fill(-1);
	oA.fill(-1);


	if(ring_id == 0)
	{
		Pside=0;
		Pup=2.0*datum::pi*((float)ring_id+0.5)/NCells;
		Pdown=0;
		PTotal=Pside*2+Pup+Pdown;

		Cells.row(0)=linspace<Col<unsigned short> >(0,6,7).t();
		Rings.row(0).fill(1);
		Areas.fill(0.12788);
		Num(0)=6;

		oC.row(0)=linspace<Col<unsigned short> >(0,5,6).t();
		oR.row(0).fill(1);
		//because there is discontinuity here
		// 1/6 is to large for next ring of cell
		oA.fill(0.12788);
		oN(0)=6;
		LowCorr = 0.25/0.12788;
	}
	else
	{
		Pside=1.0;
		//Pup=2.0*datum::pi*((float)ring_id+0.5)/NCells;
		//Pdown=2.0*datum::pi*((float)ring_id-0.5)/NCells;
		Pup=1.0;
		Pdown=1.0;
		PTotal=Pside*2+Pup+Pdown;

		Cells.col(0)=linspace<Col<unsigned short> >(0,NCells-1,NCells) - 1;
		Cells(0,0)=NCells - 1;
		Rings.col(0).fill(ring_id);

		Cells.col(1)=linspace<Col<unsigned short> >(0,NCells-1,NCells) + 1;
		Cells(NCells - 1,1)= 0;
		Rings.col(1).fill(ring_id);
		Num.fill(2);

		Areas.submat(0,0,Areas.n_rows-1,1).fill(1./PTotal);

		oC.col(0)=linspace<Col<unsigned short> >(0,NCells-1,NCells) + 1;
		oC(NCells-1,0)=0;
		oR.col(0).fill(ring_id);
		oN.fill(1);
		oA.col(0).fill(1./PTotal);
		LowCorr = (NCells+6.)/(NCells);
	}
	CellW=2*datum::pi/NCells;
	HCellW = CellW*0.5;

	Upper = 0;
	Lower = 0;
}

neighbours::~neighbours()
{
	Cells.clear();
	Rings.clear();
	Areas.clear();
	Num.clear();
	Upper = NULL;
	Lower = NULL;
}

void neighbours::init_Up_Lo_adr(neighbours* rUpper, neighbours* rLower)
{
	Upper = rUpper;
	Lower = rLower;
}
/*
 * Class method which upon call fills in arrays with neighboouring cells coordinates
 */
void neighbours::find_neighbours(double lower_diff, double upper_diff)
{
	double cell_id0, cell_id1, fract0, fract1;

	if(Ring_id == 1)
	{
		Cells.col(2).fill(0);
		Rings.col(2).fill(0);
		Areas.col(2).fill(Pdown/PTotal);
		Num.fill(3);
	}
	else
	{
		Num.fill(2);// reset neighbours
		for(unsigned int id=0; id<NCells; ++id)
		{
			while(lower_diff >= Lower->NCells )
			{
				lower_diff = lower_diff - Lower->NCells;
			}
			fract0 = modf(lower_diff, &cell_id0);
			fract1 = modf(lower_diff+CellW/Lower->CellW, &cell_id1);

			if( fract1 > fract0 )
			{
				Cells(id, Num(id)) = (int)cell_id0;
				Rings(id, Num(id)) = Ring_id - 1;
				Areas(id, Num(id)) = Pdown/PTotal;
				Num(id)+=1;
			}
			else
			{
				Cells(id, Num(id)) = (int)cell_id0;
				Rings(id, Num(id)) = Ring_id - 1;
				Areas(id, Num(id)) = (1.-fract0)*Lower->CellW/CellW*Pdown/PTotal;
				Num(id)+=1;

				if(cell_id1 == Lower->NCells)
					Cells(id, Num(id))=0;
				else
					Cells(id, Num(id))= (int)cell_id1;
				Rings(id, Num(id))= Ring_id - 1;
				Areas(id, Num(id))=fract1*Lower->CellW/CellW*Pdown/PTotal;
				Num(id)+=1;
			}
			lower_diff += CellW/Lower->CellW;
		}
	}

	oN.fill(1);
	for(unsigned int id=0; id<NCells && Upper; ++id)
	{
		while(upper_diff >= Upper->NCells)
		{
			upper_diff -= Upper->NCells;
		}
		fract0 = modf(upper_diff, &cell_id0);
		fract1 = modf(upper_diff+CellW/Upper->CellW, &cell_id1);

		Cells(id,Num(id)) = (int)cell_id0;
		Rings(id,Num(id)) = Ring_id + 1;
		Areas(id,Num(id)) = (1.-fract0)*Upper->CellW/CellW*Pup/PTotal;

		oC(id,oN(id)) = Cells(id,Num(id));
		oR(id,oN(id)) = Rings(id,Num(id));
		oA(id,oN(id)) = Areas(id,Num(id));

		oN(id) +=1;
		Num(id)+=1;

		if((int)cell_id0+2 == (int)cell_id1)
		{
			if(cell_id0+1 == Upper->NCells)
				Cells(id,Num(id))=0;
			else
				Cells(id,Num(id)) = (int)cell_id0+1;
			Rings(id,Num(id)) = Ring_id + 1;
			Areas(id,Num(id)) = Upper->CellW/CellW*Pup/PTotal;

			oC(id,oN(id)) = Cells(id,Num(id));
			oR(id,oN(id)) = Rings(id,Num(id));
			oA(id,oN(id)) = Areas(id,Num(id));

			oN(id) +=1;
			Num(id)+=1;
		}
		if (cell_id1 >= Upper->NCells)
			Cells(id,Num(id)) = (int)cell_id1 - Upper->NCells;
		else
			Cells(id,Num(id)) = (int)cell_id1;
		Rings(id,Num(id)) = Ring_id + 1;
		Areas(id,Num(id)) = fract1*Upper->CellW/CellW*Pup/PTotal;

		oC(id,oN(id)) = Cells(id,Num(id));
		oR(id,oN(id)) = Rings(id,Num(id));
		oA(id,oN(id)) = Areas(id,Num(id));

		oN(id) +=1;
		Num(id)+=1;
		upper_diff += CellW/Upper->CellW;
	}
	/*Cells.save("neigh_cell_dump.dat",raw_ascii);
	Areas.save("neigh_cell_aras_dump.dat",raw_ascii);
	exit(0);*/
}
