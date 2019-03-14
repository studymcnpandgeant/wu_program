#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmredundant.h>
#include <petscdmcomposite.h>
#include <petscpf.h>
#include <petscsnes.h>
#include "RPowertsa.h"
#include "RPts.h"
#include "Mainfunc.h"
#include <cmath>

PetscErrorCode FormInitialGuess(Vec U,Vec iniphi1,Vec iniphi2,Vec iniphi3,Vec iniphi4,void* dummy)
{
  RPowertsa      *user = (RPowertsa*)dummy;
  PetscErrorCode ierr;
  PetscScalar    *u;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
  PetscInt       row,Ilast,i,j,II,n;
  PetscScalar    *initial_phi1,*initial_phi2,*initial_phi3,*initial_phi4;

  PetscFunctionBeginUser;
  ierr =VecGetArray(iniphi1,&initial_phi1);
  ierr =VecGetArray(iniphi2,&initial_phi2);
  ierr =VecGetArray(iniphi3,&initial_phi3);
  ierr =VecGetArray(iniphi4,&initial_phi4);

  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  Iend = 4680;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  n         = 45;

  for (j=1; j<104+1; j++){
    for (i=1; i<45+1; i++){
        II = (j-1)*n + i - 1 ;
        u[II] =  initial_phi1[II];
        u[II+Interval1] =  initial_phi2[II];
        u[II+Interval2] =  initial_phi3[II];
        u[II+Interval3] =  initial_phi4[II];
    }
  }

  for (j=1; j<104+1; j++){
    for (i=1; i<45+1; i++){
        II = (j-1)*n + i - 1 ;
        user->fluxold[II] =  initial_phi1[II];
        user->fluxold[II+Interval1] =  initial_phi2[II];
        user->fluxold[II+Interval2] =  initial_phi3[II];
        user->fluxold[II+Interval3] =  initial_phi4[II];
    }
  }


  user->keff = initial_phi4[4680];
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi1,&initial_phi1);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi2,&initial_phi2);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi3,&initial_phi3);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi4,&initial_phi4);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec u,Vec udot,Vec FU,void *dummy)
{
  RPowertsa      *user = (RPowertsa*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N,II,n;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       group1=1,group2=2,group3=3,group4=4,group;
  PetscInt       indexc,indexl,indexr,indexu,indexd;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  const PetscScalar    *arrayu,*arraydotu;
  PetscScalar    *fu;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,cc1,cc2,cc3,cc4,ac1,ac2,ac3,ac4,ac0;
  PetscScalar    phigroup1,phigroup2,phigroup3,phigroup4;
  PetscScalar    keff;
  PetscScalar    Tpower0,Tpower,c,e,sigf,cvv;//power related
  PetscInt       indexpower,indexpower2;
  PetscReal      sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;
  PetscReal      sigmas12,sigmas13,sigmas14,sigmas23,sigmas24,sigmas34;
  PetscReal      chisigmaf1,chisigmaf2,chisigmaf3,chisigmaf4;
  PetscReal      phi1boundary,phi2boundary,phi3boundary,phi4boundary;
  PetscInt       indexdflux,igd;//   dflux area

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(u,&arrayu);CHKERRQ(ierr);
  ierr = VecGetArrayRead(udot,&arraydotu);CHKERRQ(ierr);
  ierr = VecGetArray(FU,&fu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(u, &start, &end);CHKERRQ(ierr);
  n     = 45;
  phi1boundary = 0.0;
  phi2boundary = 0.0;
  phi3boundary = 0.0;
  phi4boundary = 0.0;
  xs = 0; ys = 0; xm = 45; ym = 104;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;

  /* residual f_lambda */
  //if (xs == 0 && ys == 0 ) { /* only first processor computes this */
  //  dotphi = 0.0;
  //  for ( j=1 ; j < 104+1 ; j++ ){
  //     for ( i=1 ; i < 45+1 ; i++ ){
  //        II = (j-1)*n + i-1;
  //        dotphi += u[II]*u[II] + u[II+Interval1]*u[II+Interval1] + u[II+Interval2]*u[II+Interval2] + u[II+Interval3]*u[II+Interval3] ;
  //      }
  //  }
  //  fu[Ilast] =  dotphi - 1.0 ;
  //}

  /* inner points */
  for (j=1; j<104+1; j++){
    for (i=1; i<45+1; i++){
      for ( group = 1; group < 5; group++)
      {
        II = (j-1)*n + i - 1 ;
      /*phi1 field*/
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      phigroup1  = arrayu[II];
      phigroup2  = arrayu[II+Interval1];
      phigroup3  = arrayu[II+Interval2];
      phigroup4  = arrayu[II+Interval3];
      keff       = user->keff;
      rightside = user->gettsrs(i,j,group,phigroup1,phigroup2,phigroup3,phigroup4,keff,t) ;
      if ( i > 1 )
      {
        Ileft = II - 1 ;
        indexl = Ileft + (group-1)*45*104;
        ac1 = user->gettsa1(i,j,group);
        cc1 = ac1 * arrayu[indexl];
      }
      if ( i < 45 )
      {
        Iright = II + 1 ;
        indexr = Iright + (group-1)*45*104;
        ac2 = user->gettsa2(i,j,group);
        cc2 = ac2 * arrayu[indexr];
      }
      if ( j > 1 )
      {
        Idown = II - 45 ;
        indexd = Idown + (group-1)*45*104 ;
        ac3 = user->gettsa3(i,j,group);
        cc3 = ac3 * arrayu[indexd];
      }
      if ( j < 104 )
      {
        Iup = II + 45 ;
        indexu = Iup + (group-1)*45*104 ;
        ac4 = user->gettsa4(i,j,group);
        cc4 = ac4 * arrayu[indexu];
      }
      indexc = (group-1)*104*45 + II;
      ac1 = user->gettsa1(i,j,group);
      ac2 = user->gettsa2(i,j,group);
      ac3 = user->gettsa3(i,j,group);
      ac4 = user->gettsa4(i,j,group);
      ac0 = user->gettsa0(i,j,group,ac1,ac2,ac3,ac4);
      fu[indexc] = ac0 * arrayu[indexc] +cc1+cc2+cc3+cc4-rightside;
      }  
    }
  }
  /* add the time coefficient */
  addTimeDerivativeCoefficient(user,fu,arraydotu);

  ierr = VecRestoreArrayRead(u,&arrayu);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(udot,&arraydotu);CHKERRQ(ierr);
  ierr = VecRestoreArray(FU,&fu);CHKERRQ(ierr);
 
  PetscLogFlops(xm*ym*4.0*27.0);

  PetscFunctionReturn(0);
}

PetscErrorCode FormIJacobian(TS ts,PetscReal t,Vec u,Vec udot,PetscReal s,Mat J,Mat B,void *dummy)
{
  RPowertsa      *user = (RPowertsa*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N,II,n;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       group,group1=1,group2=2,group3=3,group4=4;
  PetscInt       chigroup1=1,chigroup2=2,chigroup3=3;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  PetscInt       idcenter,idright,idleft,idup,iddown;
  const PetscScalar    *arrayu,*arraydotu;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,e,b,c,dd,a;
  PetscReal      sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;

  PetscFunctionBeginUser;
  //if (user->oshift == s) return 0;//This statement should require more attention
  ierr = VecGetArrayRead(u,&arrayu);CHKERRQ(ierr);
  ierr = VecGetArrayRead(udot,&arraydotu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(u, &start, &end);CHKERRQ(ierr);
  n     = 45;
  xs = 0; ys = 0; xm = 45; ym = 104;
  Iend = 4680;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;

  /* inner points */
  for (j=1; j<105; j++){
    for (i=1; i<46; i++){
      for (group = 1; group < 5; group++)
      {
        II = (j-1)*n + i -1 + (group-1)*45*104 ;
        idcenter = II;
        /*phi1 field*/
        if ( i > 1 )
      {
        idleft   = II -  1;
        b  = user->gettsa1(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idleft,&b,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 45 )
      {
        idright  = II +  1;
        c  = user->gettsa2(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idright,&c,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 1 )
      {
        iddown   = II - 45;
        dd = user->gettsa3(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&iddown,&dd,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 104 )
      {
        idup     = II + 45;
        e  = user->gettsa4(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idup,&e,INSERT_VALUES);CHKERRQ(ierr);
      }
        b = user->gettsa1(i,j,group);
        c = user->gettsa2(i,j,group);
        dd = user->gettsa3(i,j,group);
        e = user->gettsa4(i,j,group);        
        a  = 1.0 / user->v[group-1] * user->getcvv(i,j) * s + user->gettsa0(i,j,group,b,c,dd,e) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idcenter,&a,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = VecRestoreArrayRead(u,&arrayu);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(udot,&arraydotu);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (J != B) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscLogFlops(xm * ym * 4.0 * 5.0 );
  user->oshift = s;
  PetscFunctionReturn(0);

}


PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal time,Vec u,void *dummy)
{
  RPowertsa      *user = (RPowertsa*)dummy;
  PetscErrorCode ierr;
  const PetscReal      *arrayu;
  PetscInt       i,j,igd,Interval1,Interval2,Interval3,II,indexdflux;
  PetscScalar    Qfd,phigroup1,phigroup2,phigroup3,phigroup4,c31,c41;
  PetscScalar    Tpower,e,sigf,cvv;
  PetscInt       indexpower,n,group;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Timestep %D:\n",step);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"CurrentTime %g:\n",time);CHKERRQ(ierr);

  ierr = VecGetArrayRead(u,&arrayu);CHKERRQ(ierr);
  
  /* Compute the Cdflux and Cdfluxold */
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  n     = 45;

  if (time == 0.1)
  {
    user->perturbation(time);
  }

  if ( time == 0.0 )
  {
    computeInitialCdflux(user,arrayu);
  }
  else
  {
    computeLastStepCdflux(user,arrayu);
  }

  /* storage the last step flux */
  storageLastStepFlux(user,arrayu);

  /* Compute the power */
  Tpower = 0.0;
  computePower(user,arrayu,&Tpower);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Power ======= %g\n",(double)Tpower);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(u,&arrayu);CHKERRQ(ierr);

  PetscFunctionReturn(0);
  return 0;
}