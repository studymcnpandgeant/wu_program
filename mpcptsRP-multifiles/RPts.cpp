#include <iostream>
#include <petsc.h>
#include <cmath>
#include "RPts.h"
#include "RPowertsa.h"

void computeInitialCdflux(RPowertsa *user,const PetscReal *array1)
{
	PetscInt       i,j,igd,Interval1,Interval2,Interval3,II,indexdflux,n;
	PetscScalar    Qfd,phigroup1,phigroup2,phigroup3,phigroup4;

	Interval1 = 4680;
    Interval2 = 9360;
    Interval3 = 14040;
    n         = 45;
	for (j = 1; j < 104+1; j++)
    {
      for (i = 1; i < 45+1; i++)
      {
        for (igd = 1; igd < 7 ; igd++)
        {
          II = (j-1)*n + i - 1 ;
          phigroup1  = array1[II];
          phigroup2  = array1[II+Interval1];
          phigroup3  = array1[II+Interval2];
          phigroup4  = array1[II+Interval3];
          Qfd = phigroup1*user->getsigf(i,j,1)+phigroup2*user->getsigf(i,j,2)+phigroup3*user->getsigf(i,j,3)+phigroup4*user->getsigf(i,j,4);
          indexdflux = (j-1)*n + i - 1 + (igd-1)*45*104;
          user->Cdfluxold[indexdflux] = user->beta[igd-1] * Qfd / user->lamda[igd-1] ;
          user->Cdflux[indexdflux]    = user->Cdfluxold[indexdflux] ;  
        }
      }
    }


}

void computeLastStepCdflux(RPowertsa *user,const PetscReal *array1)
{
	PetscInt       i,j,igd,Interval1,Interval2,Interval3,II,indexdflux,n;
	PetscScalar    c31,c41;


	Interval1 = 4680;
    Interval2 = 9360;
    Interval3 = 14040;
    n         = 45;
	for (j = 1; j < 104+1; j++)
    {
      for (i = 1; i < 45+1; i++)
      {
        for (igd = 1; igd < 7 ; igd++)
        {
          II = (j-1)*n + i - 1 ;
          c31 = array1[II]*user->getsigf(i,j,1)+array1[II+Interval1]*user->getsigf(i,j,2)+array1[II+Interval2]*user->getsigf(i,j,3)+array1[II+Interval3]*user->getsigf(i,j,4);
          c41 = user->fluxold[II]*user->getsigf(i,j,1)+user->fluxold[II+Interval1]*user->getsigf(i,j,2)+user->fluxold[II+Interval2]*user->getsigf(i,j,3)+user->fluxold[II+Interval3]*user->getsigf(i,j,4); 
          indexdflux = (j-1)*n + i - 1 + (igd-1)*45*104;
          user->Cdflux[indexdflux] = user->beta[igd-1]/user->lamda[igd-1]*(
            ((1-(1-exp(-user->lamda[igd-1]*user->dt))/(user->lamda[igd-1]*user->dt))*c31)
              +((1-exp(-user->lamda[igd-1]*user->dt))/(user->lamda[igd-1]*user->dt)-exp(-user->lamda[igd-1]*user->dt))*c41)
            +user->Cdfluxold[indexdflux]*exp(-user->lamda[igd-1]*user->dt);
        }
      }
    }

    for (j = 1; j < 104+1; j++)
    {
      for (i = 1; i < 45+1; i++)
      {
        for (igd = 1; igd < 7 ; igd++)
        {
          II = (j-1)*n + i - 1 ; 
          indexdflux = (j-1)*n + i - 1 + (igd-1)*45*104;
          user->Cdfluxold[indexdflux] = user->Cdflux[indexdflux] ;
        }
      }
    }

}


void storageLastStepFlux(RPowertsa *user,const PetscReal *array1)
{
  PetscInt   i;
  for (i = 0; i < 45*104*4; i++)
  {
    user->fluxold[i] = array1[i];
  }

}


void computePower(RPowertsa *user,const PetscReal *array1,PetscScalar *power)
{
	PetscScalar e,sigf,cvv;
	PetscInt    i,j,group,indexpower,n;


    (*power)  = 0.0;
    e         = 3.2e-14;
    n         = 45;
    for (j = 1; j < 104+1; j++)
    {
      for (i = 1; i < 45+1; i++)
      {
        for (group = 1; group < 5 ; group++)
        {
          indexpower = (j-1)*n + i -1 + (group-1)*45*104;
          sigf = user->getsigf(i,j,group);
          cvv   = user->getcvv(i,j);
          (*power) = (*power) + array1[indexpower]* sigf *e * cvv / 1.0e6 ;
        }
      }
    }
}

void addTimeDerivativeCoefficient(RPowertsa *user,PetscScalar *f,const PetscScalar *adotu)
{
	PetscInt i,j,group,II,n,indexc;

	n = 45;
	for (j=1; j<104+1; j++){
      for (i=1; i<45+1; i++){
        for ( group = 1; group < 5; group++)
        {
          II = (j-1)*n + i - 1 ;
          indexc = (group-1)*104*45 + II;
          f[indexc] = 1.0 / user->v[group-1] *adotu[indexc] * user->getcvv(i,j) + f[indexc] ;
        }  
      }
    }
}