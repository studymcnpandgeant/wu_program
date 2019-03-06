#ifndef THRP1_H
#define THRP1_H
#include <petsc.h>
#include <cmath>
#include "TH1.h"
#include "RPower.h"

class THRP1 : public RPower, public TH1
{
public:
	THRP1(){};
	~THRP1(){};
	void settotalinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10,
		PetscScalar *array11,PetscScalar *array12,PetscScalar *array13,PetscScalar *array14,PetscScalar *array15,PetscScalar *array16,PetscScalar *array17,PetscScalar *array18,PetscScalar *array19);
	PetscScalar RPgeta1(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar RPgeta2(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar RPgeta3(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar RPgeta4(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar RPgeta0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4);
	PetscScalar RPgetrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff);
	PetscScalar RPgetsigf(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar RPgetcvv(PetscInt i, PetscInt j);
	void RPiniQf();
	void RPaddQQQ();
};



#endif