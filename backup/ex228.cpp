static const char help[] = "Solves PDE optimization problem using full-space method, treats state and adjoint variables separately.\n\n";

/*T
   requires: !single
T*/

#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmredundant.h>
#include <petscdmcomposite.h>
#include <petscpf.h>
#include <petscsnes.h>
#include "RPowerplus.h"

/*

       

*/

typedef struct
{
  PetscScalar *vv1,*vv2,*vv3;
  PetscScalar *vv4,*vv5,*vv6;
  PetscScalar *vv7,*vv8,*vv9;
} AppCtx;

extern PetscErrorCode FormInitialGuess(Vec,Vec,Vec,Vec,Vec);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  AppCtx         data;
  PetscErrorCode ierr;
  PetscInt       its;
  PetscInt       nx=45,ny=104,nglobal,nq,ndflux;
  Vec            U,FU,vlambda,vphi1,vphi2,vphi3,vphi4;
  Vec            u1,u2,u3,u4,u5,u6,u7,u8,u9;
  Vec            iniphi1,iniphi2,iniphi3,iniphi4;
  Vec            Q,vQQQ;
  Vec            DFLUX;
  PetscScalar    *arraydflux;
  PetscScalar    *qfth,*lQQQ;
  Mat            A1,A2,A3,A4;/*Matrix for storage initial guess, four groups*/
  Mat            B;//Precondition Matrix
  SNES           snes;
  RPowerplus     user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Reading vectors in binary format */
  /* Read four vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-mpcp2",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u4);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u5);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u6);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u7);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u8);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u9);CHKERRQ(ierr);
  ierr = VecLoad(u1,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = VecLoad(u2,viewer);CHKERRQ(ierr);/* grid's materials infomation (INT format) */
  ierr = VecLoad(u3,viewer);CHKERRQ(ierr);/* dx-infomation */
  ierr = VecLoad(u4,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u5,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u6,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = VecLoad(u7,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u8,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u9,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the arrays in user's appctx */
  ierr =VecGetArray(u1,&data.vv1);
  ierr =VecGetArray(u2,&data.vv2);
  ierr =VecGetArray(u3,&data.vv3);
  ierr =VecGetArray(u4,&data.vv4);
  ierr =VecGetArray(u5,&data.vv5);
  ierr =VecGetArray(u6,&data.vv6);
  ierr =VecGetArray(u7,&data.vv7);
  ierr =VecGetArray(u8,&data.vv8);
  ierr =VecGetArray(u9,&data.vv9);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n reading vector in binary from QQQ.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-mpcpQQQ",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&vQQQ);CHKERRQ(ierr);
  ierr = VecLoad(vQQQ,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr =VecGetArray(vQQQ,&lQQQ);
    for (int i = 0; i < 4680; i++)
  {
    user.QQQ[i] = lQQQ[i];
  }
  ierr = VecRestoreArray(vQQQ,&lQQQ);CHKERRQ(ierr);
  ierr = VecDestroy(&vQQQ);CHKERRQ(ierr);

  user.setinitial(data.vv1,data.vv2,data.vv3,data.vv4,data.vv5,data.vv6,data.vv7,data.vv8,data.vv9);

  /*Read the file and store them in the matrix*/
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f option");
  }
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&iniphi1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&iniphi2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&iniphi3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&iniphi4);CHKERRQ(ierr);
  ierr = VecLoad(iniphi1,viewer);CHKERRQ(ierr);
  ierr = VecLoad(iniphi2,viewer);CHKERRQ(ierr);
  ierr = VecLoad(iniphi3,viewer);CHKERRQ(ierr);
  ierr = VecLoad(iniphi4,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);



  nglobal = 4 * nx * ny + 1 ;//set the global size of the solution vectors: four groups and a egivenvalue.
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetSizes(U, PETSC_DECIDE, nglobal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(U,iniphi1,iniphi2,iniphi3,iniphi4);CHKERRQ(ierr);

  /* create nonlinear solver */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* Form Precondition Matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nglobal,nglobal);CHKERRQ(ierr);
  ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(B,5,NULL);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,B,B,FormJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes,Monitor,&user,0);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  /* Output the data */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpout20190301",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(U,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  nq   = nx * ny ;
  ierr = VecCreate(PETSC_COMM_WORLD,&Q);CHKERRQ(ierr);
  ierr = VecSetSizes(Q, PETSC_DECIDE, nq);CHKERRQ(ierr);
  ierr = VecSetFromOptions(Q);CHKERRQ(ierr);
  ierr = VecGetArray(Q,&qfth);
  for (int i = 0; i < 4680; i++)
  {
    qfth[i] =  user.Qf[i];
  }
  ierr = VecRestoreArray(Q,&qfth);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpQf",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(Q,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ndflux   = nx * ny * 6 ;
  ierr = VecCreate(PETSC_COMM_WORLD,&DFLUX);CHKERRQ(ierr);
  ierr = VecSetSizes(DFLUX, PETSC_DECIDE, ndflux);CHKERRQ(ierr);
  ierr = VecSetFromOptions(DFLUX);CHKERRQ(ierr);
  ierr = VecGetArray(DFLUX,&arraydflux);
  for (int i = 0; i < 28080; i++)
  {
    arraydflux[i] =  user.Cdflux[i];
  }
  ierr = VecRestoreArray(DFLUX,&arraydflux);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpdflux",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(DFLUX,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);


  ierr = SNESDestroy(&snes);CHKERRQ(ierr);


  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&FU);CHKERRQ(ierr);
  ierr = VecDestroy(&u1);CHKERRQ(ierr);
  ierr = VecDestroy(&u2);CHKERRQ(ierr);
  ierr = VecDestroy(&u3);CHKERRQ(ierr);
  ierr = VecDestroy(&u4);CHKERRQ(ierr);
  ierr = VecDestroy(&u5);CHKERRQ(ierr);
  ierr = VecDestroy(&u6);CHKERRQ(ierr);
  ierr = VecDestroy(&u7);CHKERRQ(ierr);
  ierr = VecDestroy(&u8);CHKERRQ(ierr);
  ierr = VecDestroy(&u9);CHKERRQ(ierr);
  ierr = VecDestroy(&Q);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi1);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi2);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi3);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi4);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


PetscErrorCode FormInitialGuess(Vec U,Vec iniphi1,Vec iniphi2,Vec iniphi3,Vec iniphi4)
{
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
  Ilast     = 18720;
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
  u[Ilast] = 0.999884;
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
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  RPowerplus     *user = (RPowerplus*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N,II,n;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       group1=1,group2=2,group3=3,group4=4,group;
  PetscInt       indexc,indexl,indexr,indexu,indexd;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  const PetscScalar    *u;
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
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(FU,&fu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
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
  if (xs == 0 && ys == 0 ) { /* only first processor computes this */
    dotphi = 0.0;
    for ( j=1 ; j < 104+1 ; j++ ){
       for ( i=1 ; i < 45+1 ; i++ ){
          II = (j-1)*n + i-1;
          dotphi += u[II]*u[II] + u[II+Interval1]*u[II+Interval1] + u[II+Interval2]*u[II+Interval2] + u[II+Interval3]*u[II+Interval3] ;
       }
    }
    fu[Ilast] =  dotphi - 1.0 ;
  }

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
      phigroup1  = u[II];
      phigroup2  = u[II+Interval1];
      phigroup3  = u[II+Interval2];
      phigroup4  = u[II+Interval3];
      keff       = u[Ilast];
      rightside = user->getrs(i,j,group,phigroup1,phigroup2,phigroup3,phigroup4,keff) ;
      if ( i > 1 )
      {
        Ileft = II - 1 ;
        indexl = Ileft + (group-1)*45*104;
        ac1 = user->geta1(i,j,group);
        cc1 = ac1 * u[indexl];
      }
      if ( i < 45 )
      {
        Iright = II + 1 ;
        indexr = Iright + (group-1)*45*104;
        ac2 = user->geta2(i,j,group);
        cc2 = ac2 * u[indexr];
      }
      if ( j > 1 )
      {
        Idown = II - 45 ;
        indexd = Idown + (group-1)*45*104 ;
        ac3 = user->geta3(i,j,group);
        cc3 = ac3 * u[indexd];
      }
      if ( j < 104 )
      {
        Iup = II + 45 ;
        indexu = Iup + (group-1)*45*104 ;
        ac4 = user->geta4(i,j,group);
        cc4 = ac4 * u[indexu];
      }
      indexc = (group-1)*104*45 + II;
      ac1 = user->geta1(i,j,group);
      ac2 = user->geta2(i,j,group);
      ac3 = user->geta3(i,j,group);
      ac4 = user->geta4(i,j,group);
      ac0 = user->geta0(i,j,group,ac1,ac2,ac3,ac4);
      fu[indexc] = ac0 * u[indexc] +cc1+cc2+cc3+cc4-rightside;
      }  
    }
  }

  Tpower = 0.0;
  e      = 3.2e-14;
  Tpower0= 2.5e8;
  for (j = 1; j < 104+1; j++)
  {
    for (i = 1; i < 45+1; i++)
    {
      for (group = 1; group < 5 ; group++)
      {
        indexpower = (j-1)*n + i -1 + (group-1)*45*104;
        sigf = user->getsigf(i,j,group);
        cvv   = user->getcvv(i,j);
        Tpower = Tpower + u[indexpower]* sigf *e * cvv / 1.0e6 ;
      }
    }
  }

  c = Tpower0 / Tpower;
  user->iniQf();

  for (j = 1; j < 104+1; j++)
  {
    for (i = 1; i < 45+1; i++)
    {
      for (group = 1; group < 5 ; group++)
      {
        indexpower = (j-1)*n + i -1 + (group-1)*45*104;
        indexpower2= (j-1)*n + i -1;
        sigf = user->getsigf(i,j,group);
        user->Qf[indexpower2] = user->Qf[indexpower2] + u[indexpower]* sigf *e * c ;
      }
    }
  }
  user->addQQQ();

  for (j = 1; j < 104+1; j++)
  {
    for (i = 1; i < 45+1; i++)
    {
      for (igd = 1; igd < 7 ; igd++)
      {
        II = (j-1)*n + i - 1 ;
        phigroup1  = u[II];
        phigroup2  = u[II+Interval1];
        phigroup3  = u[II+Interval2];
        phigroup4  = u[II+Interval3];
        indexdflux = (j-1)*n + i - 1 + (igd-1)*45*104;
        user->Cdflux[indexdflux] = user->getdflux(i,j,igd,phigroup1,phigroup2,phigroup3,phigroup4) ;
      }
    }
  }

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(FU,&fu);CHKERRQ(ierr);
 
  PetscLogFlops(xm*ym*4.0*27.0);

  PetscFunctionReturn(0);
}

PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat B,void *dummy)
{
  RPowerplus     *user = (RPowerplus*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N,II,n;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       group,group1=1,group2=2,group3=3,group4=4;
  PetscInt       chigroup1=1,chigroup2=2,chigroup3=3;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  PetscInt       idcenter,idright,idleft,idup,iddown;
  const PetscScalar    *u;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,e,b,c,dd,a;
  PetscReal      sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  n     = 45;
  xs = 0; ys = 0; xm = 45; ym = 104;
  Iend = 4680;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;

  /* last diagonal element in Precondition Matrix */
  if (xs == 0 && ys == 0 ) {
    a = 1.0;
    ierr  = MatSetValues(B,1,&Ilast,1,&Ilast,&a,INSERT_VALUES);CHKERRQ(ierr);
  }


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
        b  = user->geta1(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idleft,&b,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 45 )
      {
        idright  = II +  1;
        c  = user->geta2(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idright,&c,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 1 )
      {
        iddown   = II - 45;
        dd = user->geta3(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&iddown,&dd,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 104 )
      {
        idup     = II + 45;
        e  = user->geta4(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idup,&e,INSERT_VALUES);CHKERRQ(ierr);
      }        
        a  = user->geta0(i,j,group,b,c,dd,e) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idcenter,&a,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (J != B) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscLogFlops(xm * ym * 4.0 * 5.0 );
  PetscFunctionReturn(0);

}


PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *dummy)
{
  RPowerplus     *user = (RPowerplus*)dummy;
  PetscErrorCode ierr;
  Vec            U;
  const PetscReal      *u;
  PetscReal      keff;
  PetscInt       Ilast;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);

  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  Ilast = 18720;

  keff = u[Ilast] ;

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, Keff ======= %g\n",its,(double)keff);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
