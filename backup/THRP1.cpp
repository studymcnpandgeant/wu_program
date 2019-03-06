#include "TH1.h"
#include "RPower.h"
#include "THRP1.h"
#include <iostream>
#include <petsc.h>
#include <cmath>


void THRP1::settotalinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10,
		PetscScalar *array11,PetscScalar *array12,PetscScalar *array13,PetscScalar *array14,PetscScalar *array15,PetscScalar *array16,PetscScalar *array17,PetscScalar *array18,PetscScalar *array19)
		{
			TH1::setTHinitial(array1,array2,array3,array4,array5,array6,array7,array8,array9,array10);
			RPower::setinitial(array11,array12,array13,array14,array15,array16,array17,array18,array19);
		}

PetscScalar THRP1::RPgeta1(PetscInt i, PetscInt j, PetscInt group)
{
	RPower::geta1(i, j, group);
}

PetscScalar THRP1::RPgeta2(PetscInt i, PetscInt j, PetscInt group)
{
	RPower::geta2(i, j, group);
}

PetscScalar THRP1::RPgeta3(PetscInt i, PetscInt j, PetscInt group)
{
	RPower::geta3(i, j, group);
}

PetscScalar THRP1::RPgeta4(PetscInt i, PetscInt j, PetscInt group)
{
	RPower::geta4(i, j, group);
}

PetscScalar THRP1::RPgeta0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4)
{
	RPower::geta0(i, j, group, a1, a2, a3, a4);
}

PetscScalar THRP1::RPgetrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff)
{
	RPower::getrs(i, j, group, flux1, flux2, flux3, flux4, keff);
}


PetscScalar THRP1::RPgetsigf(PetscInt i, PetscInt j, PetscInt group)
{
	RPower::getsigf(i, j, group);
}

PetscScalar THRP1::RPgetcvv(PetscInt i, PetscInt j)
{
	RPower::getcvv(i, j);
}

void THRP1::RPiniQf()
{
	RPower::iniQf();
}

void THRP1::RPaddQQQ()
{
	RPower::addQQQ();
}