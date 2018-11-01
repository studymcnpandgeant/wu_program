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
	PetscReal dxright,dx;
	PetscReal D_E,D_P;
	PetscReal ae,aright;
	PetscInt  id,index1,idright,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dx        = dr[j+1];
	dxright   = dr[j+2];
	if (j < (16-1) )
	{
		idright = corearray[II+1];
	    index2 = (idright-1)*36 + (group -1)*9 ;
	    D_E = 1.0 / (3.0 * materialarray[index2]);
		ae = - 1.0 / ( dx*dxright / (2.0 * D_E) + dx*dx / ( 2.0 * D_P )) ;
		return ae;
	}
	else
	{
		aright = - 1.0 / ( dx*dxright / (2.0 * D_P) + dx*dx / ( 2.0 * D_P )) ;
		return aright;
	}
}

PetscScalar Grid::getaw(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dxleft,dx;
	PetscReal D_W,D_P;
	PetscReal aw,aleft;
	PetscInt  id,index1,idleft,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dx        = dr[j+1];
	dxleft    = dr[j];
	if (j > 0 )
	{
		idleft = corearray[II-1];
	    index2 = (idleft-1)*36 + (group -1)*9 ;
	    D_W = 1.0 / (3.0 * materialarray[index2]);
		aw = - 1.0 / ( dx*dxleft / (2.0 * D_W) + dx*dx / ( 2.0 * D_P )) ;
		return aw;
	}
	else
	{
		aleft = - 1.0 / ( dx*dxleft / (2.0 * D_P) + dx*dx / ( 2.0 * D_P )) ;
		return aleft;
	}
}

PetscScalar Grid::getan(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dyup,dy;
	PetscReal D_N,D_P;
	PetscReal an,aup;
	PetscInt  id,index1,idup,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dy        = dz[i+1];
	dyup      = dz[i+2];
	if (i < (35-1) )
	{
		idup   = corearray[II+16];
	    index2 = (idup-1)*36 + (group -1)*9 ;
	    D_N = 1.0 / (3.0 * materialarray[index2]);
		an = - 1.0 / ( dy*dyup / (2.0 * D_N) + dy*dy / ( 2.0 * D_P )) ;
		return an;
	}
	else
	{
		aup = - 1.0 / ( dy*dyup / (2.0 * D_P) + dy*dy / ( 2.0 * D_P )) ;
		return aup;
	}
}

PetscScalar Grid::getas(PetscInt i,PetscInt j,PetscInt II,PetscInt group)
{
	PetscReal dydown,dy;
	PetscReal D_S,D_P;
	PetscReal as,adown;
	PetscInt  id,index1,iddown,index2;

	id = corearray[II];
	index1 = (id-1)*36 + (group -1)*9 ;
	D_P = 1.0 / (3.0 * materialarray[index1]);
	dy          = dz[i+1];
	dydown      = dz[i];
	if (i > 0 )
	{
		iddown   = corearray[II-16];
	    index2   = (iddown-1)*36 + (group -1)*9 ;
	    D_S = 1.0 / (3.0 * materialarray[index2]);
		as = - 1.0 / ( dy*dydown / (2.0 * D_S) + dy*dy / ( 2.0 * D_P )) ;
		return as;
	}
	else
	{
		adown = - 1.0 / ( dy*dydown / (2.0 * D_P) + dy*dy / ( 2.0 * D_P )) ;
		return adown;
	}
}