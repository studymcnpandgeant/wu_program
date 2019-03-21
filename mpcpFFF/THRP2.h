#ifndef THRP2_H
#define THRP2_H
#include <petsc.h>
#include <cmath>
#include "TH2.h"
#include "RPower2.h"

class THRP2 : public RPower2, public TH2
{
private:
	int nRP,nTH;

public:
	PetscScalar *uRP,*uTH;
	THRP2();
	~THRP2();
	void settotalinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10,
		PetscScalar *array11,PetscScalar *array12,PetscScalar *array13,PetscScalar *array14,PetscScalar *array15,PetscScalar *array16,PetscScalar *array17,PetscScalar *array18,PetscScalar *array19,PetscScalar *array20);
	void vecGetTHRPArray(const PetscScalar *u);
};



#endif