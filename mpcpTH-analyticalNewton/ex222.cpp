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
#include "TH.h"

/*

       

*/

typedef struct
{
  PetscScalar *vv1,*vv2,*vv3;
  PetscScalar *vv4,*vv5,*vv6;
  PetscScalar *vv7,*vv8,*vv9;
  PetscScalar *vv10;
  PetscScalar *vv21,*vv22,*vv23,*vv24,*vv25,*vv26;
} AppCtx;

extern PetscErrorCode FormInitialGuess(Vec,Vec,Vec,Vec,Vec,Vec);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
// extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
// extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  AppCtx         data;
  PetscErrorCode ierr;
  PetscInt       its;
  PetscInt       nx=22,ny=76,nglobal,nq;
  Vec            U,FU,vlambda,vphi1,vphi2,vphi3,vphi4;
  Vec            u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;
  Vec            u21,u22,u23,u24,u25,u26;
  // Vec            Q,vQQQ;
  // PetscScalar    *qfth,*lQQQ;
  // Mat            A1,A2,A3,A4;Matrix for storage initial guess, four groups
  Mat            B;//Precondition Matrix
  SNES           snes;
  TH             user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Reading vectors in binary format */
  /* Read four vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THdata1",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u4);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u5);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u6);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u7);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u8);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u9);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u10);CHKERRQ(ierr);
  ierr = VecLoad(u1,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = VecLoad(u2,viewer);CHKERRQ(ierr);/* grid's materials infomation (INT format) */
  ierr = VecLoad(u3,viewer);CHKERRQ(ierr);/* dx-infomation */
  ierr = VecLoad(u4,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u5,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u6,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = VecLoad(u7,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u8,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u9,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = VecLoad(u10,viewer);CHKERRQ(ierr);/* y-infomation */
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
  ierr =VecGetArray(u10,&data.vv10);

  user.setinitial(data.vv1,data.vv2,data.vv3,data.vv4,data.vv5,data.vv6,data.vv7,data.vv8,data.vv9,data.vv10);

  /*Read the file and store them in the matrix*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading THini vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THini",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u21);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u22);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u23);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u24);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u25);CHKERRQ(ierr);
  ierr = VecLoad(u21,viewer);CHKERRQ(ierr);//  THP
  ierr = VecLoad(u22,viewer);CHKERRQ(ierr);//  THTf
  ierr = VecLoad(u23,viewer);CHKERRQ(ierr);//  THTs
  ierr = VecLoad(u24,viewer);CHKERRQ(ierr);//  THU
  ierr = VecLoad(u25,viewer);CHKERRQ(ierr);//  THV
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading THQfth vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THQfth",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u26);CHKERRQ(ierr);
  ierr = VecLoad(u26,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the two array */
  ierr =VecGetArray(u21,&data.vv21);
  ierr =VecGetArray(u22,&data.vv22);
  ierr =VecGetArray(u23,&data.vv23);
  ierr =VecGetArray(u24,&data.vv24);
  ierr =VecGetArray(u25,&data.vv25);
  ierr =VecGetArray(u26,&data.vv26);

  user.setTHinitial(data.vv21,data.vv22,data.vv23,data.vv24,data.vv25,data.vv26);

  nglobal = 5 * nx * ny ;//set the global size of the solution vectors: five field.
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetSizes(U, PETSC_DECIDE, nglobal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(U,u21,u22,u23,u24,u25);CHKERRQ(ierr);

  /* create nonlinear solver */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  /* Form Precondition Matrix */
  // ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  // ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nglobal,nglobal);CHKERRQ(ierr);
  // ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
  // ierr = MatSeqAIJSetPreallocation(B,5,NULL);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  // ierr = SNESSetJacobian(snes,B,B,FormJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  // ierr = SNESMonitorSet(snes,Monitor,&user,0);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  /* Output the data */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_THout",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(U,viewer);CHKERRQ(ierr);
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
  ierr = VecDestroy(&u10);CHKERRQ(ierr);
  ierr = VecDestroy(&u21);CHKERRQ(ierr);
  ierr = VecDestroy(&u22);CHKERRQ(ierr);
  ierr = VecDestroy(&u23);CHKERRQ(ierr);
  ierr = VecDestroy(&u24);CHKERRQ(ierr);
  ierr = VecDestroy(&u25);CHKERRQ(ierr);
  ierr = VecDestroy(&u26);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


PetscErrorCode FormInitialGuess(Vec U,Vec u21,Vec u22,Vec u23,Vec u24,Vec u25)
{
  PetscErrorCode ierr;
  PetscScalar    *u;
  PetscInt       start,end,Iend,Interval1,Interval2,Interval3,Interval4;
  PetscInt       row,i,j,II,n;
  PetscScalar    *initial_THP,*initial_THTf,*initial_THTs,*initial_THU,*initial_THV;

  PetscFunctionBeginUser;
  ierr =VecGetArray(u21,&initial_THP);
  ierr =VecGetArray(u22,&initial_THTf);
  ierr =VecGetArray(u23,&initial_THTs);
  ierr =VecGetArray(u24,&initial_THU);
  ierr =VecGetArray(u25,&initial_THV);

  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  Iend = 1672;
  Interval1 = 1672;
  Interval2 = 3344;
  Interval3 = 5016;
  Interval4 = 6688;
  n         = 22;

  for (j=1; j<76+1; j++){
    for (i=1; i<22+1; i++){
        II = (j-1)*n + i - 1 ;
        u[II+Interval4] =  initial_THP[II];
        u[II+Interval1] =  initial_THTf[II];
        u[II]           =  initial_THTs[II];
        u[II+Interval2] =  initial_THU[II];
        u[II+Interval3] =  initial_THV[II];
    }
  }
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(u21,&initial_THP);CHKERRQ(ierr);
  ierr = VecRestoreArray(u22,&initial_THTf);CHKERRQ(ierr);
  ierr = VecRestoreArray(u23,&initial_THTs);CHKERRQ(ierr);
  ierr = VecRestoreArray(u24,&initial_THU);CHKERRQ(ierr);
  ierr = VecRestoreArray(u25,&initial_THV);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  TH             *user = (TH*)dummy;
  PetscErrorCode ierr;
  PetscInt       i,N,II,n;
  PetscInt       j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       indexc,indexl,indexr,indexu,indexd,indexTH;
  PetscInt       start,end,Interval1,Interval2,Interval3,Interval4;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  const PetscScalar    *u;
  PetscScalar    *fu;
  PetscScalar    rightside,cc1,cc2,cc3,cc4,ac1,ac2,ac3,ac4,ac0;
  PetscScalar    P,Tf,Ts,UU,V;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(FU,&fu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  n     = 22;
  Interval1 = 1672;
  Interval2 = 3344;
  Interval3 = 5016;
  Interval4 = 6688;

  /* write the current variables */
  for (i = 0; i < 1672; i++)
  {
    user->THTs[i] = u[i];
    user->THTf[i] = u[i+Interval1];
    user->THU[i]  = u[i+Interval2];
    user->THV[i]  = u[i+Interval3];
    user->THP[i]  = u[i+Interval4];

  }

  /*set parameter*/
  for (i = 1; i < 23; i++)
  {
    for (j = 12; j < 88; j++)
    {
      indexTH = (j-12)*22+i-1;
      P       = user->THP[indexTH];
      Tf      = user->THTf[indexTH];
      Ts      = user->THTs[indexTH];
      UU      = user->THU[indexTH];
      V       = user->THV[indexTH];
      user->setparameter(i,j,Ts,Tf,P,UU,V);
    }
  }

  /* set boundary coefficient */
  for (i = 1; i < 23; i++)
  {
    for (j = 12; j < 88; j++)
    {
      user->setboundarycoefficient(i,j);
    }
  }

  /* inner points */
  for (j=12; j<87+1; j++){
    for (i=1; i<22+1; i++){
      
      /*Ts field*/
      II = (j-12)*n + i - 1 ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTsa1(i,j);
        cc1 = ac1 * u[Ileft];
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getTsa2(i,j);
        cc2 = ac2 * u[Iright];
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getTsa3(i,j);
        cc3 = ac3 * u[Idown];
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getTsa4(i,j);
        cc4 = ac4 * u[Iup];
      }
      ac1 = user->getTsa1(i,j);
      ac2 = user->getTsa2(i,j);
      ac3 = user->getTsa3(i,j);
      ac4 = user->getTsa4(i,j);
      ac0 = user->getTsa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getTsrs(i,j,ac2,ac3,ac4);
      fu[II] = ac0 * u[II] +cc1+cc2+cc3+cc4-rightside;

      /* Tf field */ 
      II = (j-12)*n + i - 1 ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTfa1(i,j);
        cc1 = ac1 * u[Ileft+Interval1];
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getTfa2(i,j);
        cc2 = ac2 * u[Iright+Interval1];
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getTfa3(i,j);
        cc3 = ac3 * u[Idown+Interval1];
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getTfa4(i,j);
        cc4 = ac4 * u[Iup+Interval1];
      }
      ac1 = user->getTfa1(i,j);
      ac2 = user->getTfa2(i,j);
      ac3 = user->getTfa3(i,j);
      ac4 = user->getTfa4(i,j);
      ac0 = user->getTfa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getTfrs(i,j,ac3);
      fu[II+Interval1] = ac0 * u[II+Interval1] +cc1+cc2+cc3+cc4-rightside;

      /* U field */
      II = (j-12)*n + i - 1 ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      ac0 = 1.0;
      rightside = user->getUrs(i,j);
      fu[II+Interval2] = ac0 * u[II+Interval2] +cc1+cc2+cc3+cc4-rightside;

      /* V field */
      II = (j-12)*n + i - 1 ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      ac0 = 1.0;
      rightside = user->getVrs(i,j);
      fu[II+Interval3] = ac0 * u[II+Interval3] +cc1+cc2+cc3+cc4-rightside;

      /* P field */
      II = (j-12)*n + i - 1 ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getPa1(i,j);
        cc1 = ac1 * u[Ileft+Interval4];
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getPa2(i,j);
        cc2 = ac2 * u[Iright+Interval4];
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getPa3(i,j);
        cc3 = ac3 * u[Idown+Interval4];
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getPa4(i,j);
        cc4 = ac4 * u[Iup+Interval4];
      }
      ac1 = user->getPa1(i,j);
      ac2 = user->getPa2(i,j);
      ac3 = user->getPa3(i,j);
      ac4 = user->getPa4(i,j);
      ac0 = user->getPa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getPrs(i,j,ac4);
      fu[II+Interval4] = ac0 * u[II+Interval4] +cc1+cc2+cc3+cc4-rightside;
    }
  }

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(FU,&fu);CHKERRQ(ierr);
 
  // PetscLogFlops(xm*ym*4.0*27.0);

  PetscFunctionReturn(0);
}

// PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat B,void *dummy)
// {
//   RPower         *user = (RPower*)dummy;
//   PetscErrorCode ierr;
//   PetscInt       xs,xm,i,N,II,n;
//   PetscInt       ys,ym,j,M;/*y direction index*/
//   PetscInt       xints,xinte,yints,yinte;
//   PetscInt       group,group1=1,group2=2,group3=3,group4=4;
//   PetscInt       chigroup1=1,chigroup2=2,chigroup3=3;
//   PetscInt       start,end,Iend,Interval1,Interval2,Interval3;
//   PetscInt       Ileft,Iright,Idown,Iup,Ilast;
//   PetscInt       idcenter,idright,idleft,idup,iddown;
//   const PetscScalar    *u;
//   PetscScalar    dotphi;//wu-define
//   PetscScalar    rightside,e,b,c,dd,a;
//   PetscReal      sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;

//   PetscFunctionBeginUser;
//   ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
//   ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
//   n     = 45;
//   xs = 0; ys = 0; xm = 45; ym = 104;
//   Iend = 4680;
//   Interval1 = 4680;
//   Interval2 = 9360;
//   Interval3 = 14040;
//   Ilast     = 18720;

//   /* last diagonal element in Precondition Matrix */
//   if (xs == 0 && ys == 0 ) {
//     a = 1.0;
//     ierr  = MatSetValues(B,1,&Ilast,1,&Ilast,&a,INSERT_VALUES);CHKERRQ(ierr);
//   }


//   /* inner points */
//   for (j=1; j<105; j++){
//     for (i=1; i<46; i++){
//       for (group = 1; group < 5; group++)
//       {
//         II = (j-1)*n + i -1 + (group-1)*45*104 ;
//         idcenter = II;
//         /*phi1 field*/
//         if ( i > 1 )
//       {
//         idleft   = II -  1;
//         b  = user->geta1(i,j,group) ;
//         ierr  = MatSetValues(B,1,&idcenter,1,&idleft,&b,INSERT_VALUES);CHKERRQ(ierr);
//       }
//       if ( i < 45 )
//       {
//         idright  = II +  1;
//         c  = user->geta2(i,j,group) ;
//         ierr  = MatSetValues(B,1,&idcenter,1,&idright,&c,INSERT_VALUES);CHKERRQ(ierr);
//       }
//       if ( j > 1 )
//       {
//         iddown   = II - 45;
//         dd = user->geta3(i,j,group) ;
//         ierr  = MatSetValues(B,1,&idcenter,1,&iddown,&dd,INSERT_VALUES);CHKERRQ(ierr);
//       }
//       if ( j < 104 )
//       {
//         idup     = II + 45;
//         e  = user->geta4(i,j,group) ;
//         ierr  = MatSetValues(B,1,&idcenter,1,&idup,&e,INSERT_VALUES);CHKERRQ(ierr);
//       }        
//         a  = user->geta0(i,j,group,b,c,dd,e) ;
//         ierr  = MatSetValues(B,1,&idcenter,1,&idcenter,&a,INSERT_VALUES);CHKERRQ(ierr);
//       }
//     }
//   }

//   ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
//   ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   if (J != B) {
//     ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//     ierr = MatAssemblyEnd  (J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//   }
//   PetscLogFlops(xm * ym * 4.0 * 5.0 );
//   PetscFunctionReturn(0);

// }


// PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *dummy)
// {
//   RPower         *user = (RPower*)dummy;
//   PetscErrorCode ierr;
//   Vec            U;
//   const PetscReal      *u;
//   PetscReal      keff;
//   PetscInt       Ilast;

//   PetscFunctionBeginUser;
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
//   ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);

//   ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
//   Ilast = 18720;

//   keff = u[Ilast] ;

//   ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, Keff ======= %g\n",its,(double)keff);CHKERRQ(ierr);

//   PetscFunctionReturn(0);
// }
