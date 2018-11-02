#include "Grid.h"
#include <iostream>
#include <petsc.h>

Grid::Grid()
{
	n1 = 6480;
	materialarray = new PetscScalar [n1];
	n2 = 560;
	corearray = new PetscScalar [n2];
	n3 = 18;//16+2
	dr = new PetscScalar [n3];
	n4 = 37;//35+2
	dz = new PetscScalar [n4];
	n5 = 18;
	rcor = new PetscScalar [n5];
	n6 = 36;
	zcor = new PetscScalar [n6];

}

Grid::~Grid()
{
	delete []materialarray;
	delete []corearray;
	delete []dr;
	delete []dz;
	delete []rcor;
	delete []zcor;
}

void Grid::setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6)
{
	for (int i = 0; i < 6480; i++)
	{
		materialarray[i] = array1[i];
	}

	for (int i = 0; i < 560; i++)
	{
		corearray[i] = array2[i];
	}

	for (int i = 0; i < 18; i++)
	{
		dr[i] = array3[i];
	}

	for (int i = 0; i < 37; i++)
	{
		dz[i] = array4[i];
	}

	for (int i = 0; i < 18; i++)
	{
		rcor[i] = array5[i];
	}

	for (int i = 0; i < 36; i++)
	{
		zcor[i] = array6[i];
	}
}

PetscScalar Grid::getmaterialarray(PetscInt a)
{
	return materialarray[a];
}

PetscScalar Grid::getae(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dxright,dx,dx_plus,dx_minus;
	PetscReal D_E,D_P;
	PetscReal ae,aright;
	PetscInt  id,index1,idright,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dx        = dr[j+1];
	dxright   = dr[j+2];
	dx_minus  = dr[j+1] / 2.0 ;
	dx_plus   = dr[j+2] / 2.0 ;

	if (j < (16-1) )
	{
		idright = corearray[II+1];
	    index2 = (idright-1)*36 + (group -1)*9 ;
	    D_E = 1.0 / (3.0 * materialarray[index2]);
		ae = - 1.0 / ( dx*dx_plus / D_E + dx*dx_minus / D_P ) ;
		return ae;
	}
	else
	{
		aright = - 1.0 / ( dx*dx_plus / D_P + dx*dx_minus / D_P ) ;
		return aright;
	}
}

PetscScalar Grid::getaw(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dxleft,dx,dx_plus,dx_minus;
	PetscReal D_W,D_P;
	PetscReal aw,aleft;
	PetscInt  id,index1,idleft,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dx        = dr[j+1];
	dxleft    = dr[j];
	dx_plus   = dr[j+1] / 2.0 ;
	dx_minus  = dr[j] / 2.0 ;
	if (j > 0 )
	{
		idleft = corearray[II-1];
	    index2 = (idleft-1)*36 + (group -1)*9 ;
	    D_W = 1.0 / (3.0 * materialarray[index2]);
		aw = - 1.0 / ( dx*dx_minus / D_W + dx*dx_plus / D_P ) ;
		return aw;
	}
	else
	{
		aleft = - 1.0 / ( dx*dx_minus / D_P + dx*dx_plus / D_P ) ;
		return aleft;
	}
}

PetscScalar Grid::getan(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dyup,dy,dy_plus,dy_minus;
	PetscReal D_N,D_P;
	PetscReal an,aup;
	PetscInt  id,index1,idup,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dy        = dz[i+1];
	dyup      = dz[i+2];
	dy_plus   = dz[i+2] / 2.0 ;
	dy_minus  = dz[i+1] / 2.0 ;
	if (i < (35-1) )
	{
		idup   = corearray[II+16];
	    index2 = (idup-1)*36 + (group -1)*9 ;
	    D_N = 1.0 / (3.0 * materialarray[index2]);
		an = - 1.0 / ( dy*dy_plus / D_N + dy*dy_minus / D_P ) ;
		return an;
	}
	else
	{
		aup = - 1.0 / ( dy*dy_plus / D_P + dy*dy_minus / D_P ) ;
		return aup;
	}
}

PetscScalar Grid::getas(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dydown,dy,dy_plus,dy_minus;
	PetscReal D_S,D_P;
	PetscReal as,adown;
	PetscInt  id,index1,iddown,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dy          = dz[i+1];
	dydown      = dz[i];
	dy_plus     = dz[i+1] / 2.0 ;
	dy_minus    = dz[i] /2.0 ;
	if (i > 0 )
	{
		iddown   = corearray[II-16];
	    index2   = (iddown-1)*36 + (group -1)*9 ;
	    D_S = 1.0 / (3.0 * materialarray[index2]);
		as = - 1.0 / ( dy*dy_minus / D_S + dy*dy_plus / D_P ) ;
		return as;
	}
	else
	{
		adown = - 1.0 / ( dy*dy_minus / D_P + dy*dy_plus / D_P ) ;
		return adown;
	}
}

PetscScalar Grid::getchisigmaf(PetscInt II, PetscInt chi,PetscInt group)
{
	PetscReal *chi1 = new PetscReal [3] ;
	PetscInt  id,index;
	PetscReal chisigmaf;

	chi1[0] = 0.98439;
    chi1[1] = 0.015616;
    chi1[2] = 0.67690e-7;
    id = corearray[II];
    index = (id - 1) * 36 + 2 + (group-1)*9;
    chisigmaf = chi1[chi-1] * materialarray[index] ;
    delete [] chi1;
    return chisigmaf;
}

PetscScalar Grid::getsigmat(PetscInt II,PetscInt group)
{
	PetscInt  id,index;
	PetscReal sigmat;

	id = corearray[II];
	index = (id-1)*36 + 3 + (group-1)*9 ;
	sigmat = materialarray[index];
	return sigmat;
}

PetscScalar Grid::getsigmas(PetscInt II,PetscInt group11,PetscInt group22)
{
	PetscInt  id,index;
	PetscReal sigmas;

	id = corearray[II];
	switch (group11)
	{
		case 1 : index = (id-1)*36 + 5 + (group22-1)*10 ;
		         break;
		case 2 : index = (id-1)*36 + 14 + (group22-2)*10 ;
		         break;
		case 3 : index = (id-1)*36 + 23 + (group22-3)*10 ;
		         break;
		case 4 : index = (id-1)*36 + 32;
		         break;
	}
	sigmas = materialarray[index];
	return sigmas;
}

PetscScalar Grid::getrzae(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal re,drr,dre_plus,dre_minus;
	PetscReal D_E,D_P;
	PetscReal ae,aright;
	PetscInt  id,index1,idright,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	re         = rcor[j+1];
	drr        = dr[j+1];
	dre_minus  = dr[j+1] / 2.0 ;
	dre_plus   = dr[j+2] / 2.0 ;

	if (j < (16-1) )
	{
		idright = corearray[II+1];
	    index2 = (idright-1)*36 + (group -1)*9 ;
	    D_E = 1.0 / (3.0 * materialarray[index2]);
		ae = - re / ( drr * dre_plus / D_E + drr * dre_minus / D_P ) ;
		return ae;
	}
	else
	{
		aright = - re / ( drr * dre_plus / D_P + drr * dre_minus / D_P ) ;
		return aright;
	}
}

PetscScalar Grid::getrzaw(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal rw,drr,drw_plus,drw_minus;
	PetscReal D_W,D_P;
	PetscReal aw,aleft;
	PetscInt  id,index1,idleft,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	rw         = rcor[j];
	drr        = dr[j+1];
	drw_plus   = dr[j+1] / 2.0 ;
	drw_minus  = dr[j] / 2.0 ;

	if (j > 0 )
	{
		idleft = corearray[II-1];
	    index2 = (idleft-1)*36 + (group -1)*9 ;
	    D_W = 1.0 / (3.0 * materialarray[index2]);
		aw = - rw / ( drr * drw_minus / D_W + drr * drw_plus / D_P ) ;
		return aw;
	}
	else
	{
		aleft = 0.0 ;/* cylindrical-coordinate system */
		return aleft;
	}
}

PetscScalar Grid::getrzan(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
    PetscReal rp,dzz,dzn_plus,dzn_minus;
	PetscReal D_N,D_P;
	PetscReal an,aup;
	PetscInt  id,index1,idup,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
    rp         = (rcor[j+1] + rcor[j]) / 2.0;
	dzz        = dz[i+1];
	dzn_plus   = dz[i+2] / 2.0 ;
	dzn_minus  = dz[i+1] / 2.0 ;

	if (i < (35-1) )
	{
		idup   = corearray[II+16];
	    index2 = (idup-1)*36 + (group -1)*9 ;
	    D_N = 1.0 / (3.0 * materialarray[index2]);
		an = - rp / ( dzz * dzn_plus / D_N + dzz * dzn_minus / D_P ) ;
		return an;
	}
	else
	{
		aup = - rp / ( dzz * dzn_plus / D_P + dzz * dzn_minus / D_P ) ;
		return aup;
	}	
}

PetscScalar Grid::getrzas(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal rp,dzz,dzs_plus,dzs_minus;
	PetscReal D_S,D_P;
	PetscReal as,adown;
	PetscInt  id,index1,iddown,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	rp           = (rcor[j+1] + rcor[j]) / 2.0 ;
	dzz          = dz[i+1];
	dzs_plus     = dz[i+1] / 2.0 ;
	dzs_minus    = dz[i] /2.0 ;

	if (i > 0 )
	{
		iddown   = corearray[II-16];
	    index2   = (iddown-1)*36 + (group -1)*9 ;
	    D_S = 1.0 / (3.0 * materialarray[index2]);
		as = - rp / ( dzz * dzs_minus / D_S + dzz * dzs_plus / D_P ) ;
		return as;
	}
	else
	{
		adown = - rp / ( dzz * dzs_minus / D_P + dzz * dzs_plus / D_P ) ;
		return adown;
	}
}

PetscScalar Grid::getrzchisigmaf(PetscInt j,PetscInt II, PetscInt chi,PetscInt group)
{
	PetscReal *chi1 = new PetscReal [3] ;
	PetscInt  id,index;
	PetscReal chisigmaf;
	PetscReal rp;

	rp      = (rcor[j+1] + rcor[j]) / 2.0 ;
	chi1[0] = 0.98439;
    chi1[1] = 0.015616;
    chi1[2] = 0.67690e-7;
    id = corearray[II];
    index = (id - 1) * 36 + 2 + (group-1)*9;
    chisigmaf = chi1[chi-1] * materialarray[index] ;
    delete [] chi1;
    return chisigmaf*rp ;
}

PetscScalar Grid::getrzsigmat(PetscInt j,PetscInt II,PetscInt group)
{
	PetscInt  id,index;
	PetscReal sigmat;
	PetscReal rp;

	rp     = (rcor[j+1] + rcor[j]) / 2.0 ;
	id     = corearray[II];
	index  = (id-1)*36 + 3 + (group-1)*9 ;
	sigmat = materialarray[index];
	return sigmat*rp ;
}

PetscScalar Grid::getrzsigmas(PetscInt j,PetscInt II,PetscInt group11,PetscInt group22)
{
	PetscInt  id,index;
	PetscReal sigmas;
	PetscReal rp;

	rp = (rcor[j+1] + rcor[j]) / 2.0 ;
	id = corearray[II];
	switch (group11)
	{
		case 1 : index = (id-1)*36 + 5 + (group22-1)*10 ;
		         break;
		case 2 : index = (id-1)*36 + 14 + (group22-2)*10 ;
		         break;
		case 3 : index = (id-1)*36 + 23 + (group22-3)*10 ;
		         break;
		case 4 : index = (id-1)*36 + 32;
		         break;
	}
	sigmas = materialarray[index];
	return sigmas*rp ;
}