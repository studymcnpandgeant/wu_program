#ifndef THRP3_H
#define THRP3_H
#include <petsc.h>
#include <cmath>
#include "TH3.h"
#include "RPower3.h"

class THRP3 : public RPower3, public TH3
{
private:
	int nRP,nTH;

public:
	PetscScalar *uRP,*uTH;
	THRP3();
	~THRP3();
	void settotalinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10,
		PetscScalar *array11,PetscScalar *array12,PetscScalar *array13,PetscScalar *array14,PetscScalar *array15,PetscScalar *array16,PetscScalar *array17,PetscScalar *array18,PetscScalar *array19,PetscScalar *array20);
	void vecGetTHRPArray(const PetscScalar *u);
};



#endif