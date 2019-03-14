#ifndef MAINFUNC_H
#define MAINFUNC_H
#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmredundant.h>
#include <petscdmcomposite.h>
#include <petscpf.h>
#include <petscsnes.h>
#include "RPowertsa.h"
#include "RPts.h"
#include <cmath>

PetscErrorCode FormInitialGuess(Vec,Vec,Vec,Vec,Vec,void*);
PetscErrorCode FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);
PetscErrorCode FormIJacobian(TS,PetscReal,Vec,Vec,PetscReal,Mat,Mat,void*);
PetscErrorCode Monitor(TS,PetscInt,PetscReal,Vec,void*);


#endif