#ifndef RPTS_H
#define RPTS_H
#include <petsc.h>
#include "RPowertsa.h"


void computeInitialCdflux(RPowertsa *user,const PetscReal *array1);
void computeLastStepCdflux(RPowertsa *user,const PetscReal *array1);
void storageLastStepFlux(RPowertsa *user,const PetscReal *array1);
void computePower(RPowertsa *user,const PetscReal *array1,PetscScalar *power);
void addTimeDerivativeCoefficient(RPowertsa *user,PetscScalar *f,const PetscScalar *adotu);



#endif