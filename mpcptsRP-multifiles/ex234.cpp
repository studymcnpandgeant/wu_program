static const char help[] = "Solves PDE optimization problem using full-space method, treats state and adjoint variables separately.\n\n";

/*T
   requires: !single
T*/

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

/*

       

*/

typedef struct
{
  PetscScalar *vv1,*vv2,*vv3;
  PetscScalar *vv4,*vv5,*vv6;
  PetscScalar *vv7,*vv8,*vv9;
} AppCtx;


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
  TS             ts;
  RPowertsa      user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;
  PetscReal      dt,ftime;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  dt = 0.1;

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



  nglobal = 4 * nx * ny ;//set the global size of the solution vectors: four groups and a egivenvalue.
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetSizes(U, PETSC_DECIDE, nglobal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  

  /* create nonlinear solver */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,U);CHKERRQ(ierr);

  /* Form Precondition Matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nglobal,nglobal);CHKERRQ(ierr);
  ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(B,5,NULL);CHKERRQ(ierr);

  ierr = TSSetIFunction(ts,NULL,FormIFunction,&user);CHKERRQ(ierr);

  ierr = TSMonitorSet(ts,Monitor,&user,NULL);CHKERRQ(ierr);
  //ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  ierr = TSSetIJacobian(ts,B,B,FormIJacobian,&user);CHKERRQ(ierr);
  user.oshift = PETSC_MIN_REAL;
  //ierr = SNESSetJacobian(snes,B,B,FormJacobian,&user);CHKERRQ(ierr);
  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(U,iniphi1,iniphi2,iniphi3,iniphi4,&user);CHKERRQ(ierr);

  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);

  ierr = TSSetTimeStep(ts,dt);CHKERRQ(ierr);

  ierr = TSSetMaxSteps(ts,30);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts,1e12);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);


  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = TSSetUp(ts);CHKERRQ(ierr);

  ierr = TSSolve(ts,U);CHKERRQ(ierr);
  ierr = TSGetSolveTime(ts,&ftime);CHKERRQ(ierr);

  ierr = TSGetStepNumber(ts,&its);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of pseudo timesteps = %D final time %4.2e\n",its,(double)ftime);CHKERRQ(ierr);

  /* Output the data */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpout20190305",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
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
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpdflux20190305",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(DFLUX,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);


  ierr = TSDestroy(&ts);CHKERRQ(ierr);


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