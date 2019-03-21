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
#include "RPower2.h"
#include "TH2.h"
#include "THRP2.h"

/*

       

*/

typedef struct
{
  PetscScalar *vv1,*vv2,*vv3;
  PetscScalar *vv4,*vv5,*vv6;
  PetscScalar *vv7,*vv8,*vv9;
  PetscScalar *vv10;
  PetscScalar *vv21,*vv22,*vv23,*vv24,*vv25,*vv26,*vv27;
  PetscScalar *vv101,*vv102,*vv103;
  PetscScalar *vv104,*vv105,*vv106;
  PetscScalar *vv107,*vv108,*vv109;
  PetscScalar *vv110;
} AppCtx;

extern PetscErrorCode FormInitialGuess(Vec,Vec,Vec,Vec,Vec,Vec,Vec,Vec,Vec,Vec);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  AppCtx         data;
  PetscErrorCode ierr;
  PetscInt       its;
  PetscInt       nx=45,ny=104,THnx=22,THny=76,nglobal,nq;
  Vec            U,FU;
  Vec            u1,u2,u3,u4,u5,u6,u7,u8,u9,u10;//TH parameters
  Vec            u21,u22,u23,u24,u25,u26,u27;//TH initial guess and Q 
  Vec            u101,u102,u103,u104,u105,u106,u107,u108,u109,u110;//RP parameters
  Vec            iniphi1,iniphi2,iniphi3,iniphi4;
  Vec            Q,vQQQ;
  PetscScalar    *qfth,*lQQQ;
  Mat            B;//Precondition Matrix
  SNES           snes;
  THRP2          user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Reading vectors in binary format */
  /* Read four vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n This is the New version of mpcpower\n \n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-mpcp2",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u101);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u102);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u103);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u104);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u105);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u106);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u107);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u108);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u109);CHKERRQ(ierr);
  ierr = VecLoad(u101,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = VecLoad(u102,viewer);CHKERRQ(ierr);/* grid's materials infomation (INT format) */
  ierr = VecLoad(u103,viewer);CHKERRQ(ierr);/* dx-infomation */
  ierr = VecLoad(u104,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u105,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u106,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = VecLoad(u107,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u108,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u109,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the arrays in user's appctx */
  ierr =VecGetArray(u101,&data.vv101);
  ierr =VecGetArray(u102,&data.vv102);
  ierr =VecGetArray(u103,&data.vv103);
  ierr =VecGetArray(u104,&data.vv104);
  ierr =VecGetArray(u105,&data.vv105);
  ierr =VecGetArray(u106,&data.vv106);
  ierr =VecGetArray(u107,&data.vv107);
  ierr =VecGetArray(u108,&data.vv108);
  ierr =VecGetArray(u109,&data.vv109);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"BCO1",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u110);CHKERRQ(ierr);
  ierr = VecLoad(u110,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr =VecGetArray(u110,&data.vv110);

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

  user.settotalinitial(data.vv1,data.vv2,data.vv3,data.vv4,data.vv5,data.vv6,data.vv7,data.vv8,data.vv9,data.vv10,data.vv101,data.vv102,data.vv103,data.vv104,data.vv105,data.vv106,data.vv107,data.vv108,data.vv109,data.vv110);

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

  /*Read the file and store them in the matrix*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading THini vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THini1",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u21);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u22);CHKERRQ(ierr);
  //ierr = VecCreate(PETSC_COMM_WORLD,&u23);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u24);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u25);CHKERRQ(ierr);
  ierr = VecLoad(u21,viewer);CHKERRQ(ierr);//  THP
  ierr = VecLoad(u22,viewer);CHKERRQ(ierr);//  THTf
  //ierr = VecLoad(u23,viewer);CHKERRQ(ierr);//  THTs
  ierr = VecLoad(u24,viewer);CHKERRQ(ierr);//  THU
  ierr = VecLoad(u25,viewer);CHKERRQ(ierr);//  THV
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading THQfth vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THQfth",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u26);CHKERRQ(ierr);
  ierr = VecLoad(u26,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading THTsbig vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-THini2",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u27);CHKERRQ(ierr);
  ierr = VecLoad(u27,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the two array */
  ierr =VecGetArray(u21,&data.vv21);
  ierr =VecGetArray(u22,&data.vv22);
  //ierr =VecGetArray(u23,&data.vv23);
  ierr =VecGetArray(u24,&data.vv24);
  ierr =VecGetArray(u25,&data.vv25);
  ierr =VecGetArray(u26,&data.vv26);
  ierr =VecGetArray(u27,&data.vv27);

  user.setTHinitial(data.vv21,data.vv22,data.vv27,data.vv24,data.vv25,data.vv26);



  nglobal = 4 * nx * ny + 1 + 4 * THnx * THny + nx * ny ;//set the global size of the solution vectors: four groups and a egivenvalue.
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetSizes(U, PETSC_DECIDE, nglobal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess first 4 neutron field, last 5 thermal field */
  ierr = FormInitialGuess(U,iniphi1,iniphi2,iniphi3,iniphi4,u21,u22,u27,u24,u25);CHKERRQ(ierr);

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
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpFFF_0320_noadd",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
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
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_mpcpQf_0320",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(Q,viewer);CHKERRQ(ierr);
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
  ierr = VecDestroy(&u101);CHKERRQ(ierr);
  ierr = VecDestroy(&u102);CHKERRQ(ierr);
  ierr = VecDestroy(&u103);CHKERRQ(ierr);
  ierr = VecDestroy(&u104);CHKERRQ(ierr);
  ierr = VecDestroy(&u105);CHKERRQ(ierr);
  ierr = VecDestroy(&u106);CHKERRQ(ierr);
  ierr = VecDestroy(&u107);CHKERRQ(ierr);
  ierr = VecDestroy(&u108);CHKERRQ(ierr);
  ierr = VecDestroy(&u109);CHKERRQ(ierr);
  ierr = VecDestroy(&u21);CHKERRQ(ierr);
  ierr = VecDestroy(&u22);CHKERRQ(ierr);
  //ierr = VecDestroy(&u23);CHKERRQ(ierr);
  ierr = VecDestroy(&u24);CHKERRQ(ierr);
  ierr = VecDestroy(&u25);CHKERRQ(ierr);
  ierr = VecDestroy(&u26);CHKERRQ(ierr);
  ierr = VecDestroy(&u27);CHKERRQ(ierr);
  ierr = VecDestroy(&Q);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi1);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi2);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi3);CHKERRQ(ierr);
  ierr = VecDestroy(&iniphi4);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


PetscErrorCode FormInitialGuess(Vec U,Vec iniphi1,Vec iniphi2,Vec iniphi3,Vec iniphi4,Vec u21,Vec u22,Vec u27,Vec u24,Vec u25)
{
  PetscErrorCode ierr;
  PetscScalar    *u;
  PetscInt       start,end,Interval1,Interval2,Interval3;
  PetscInt       THInterval1,THInterval2,THInterval3,THInterval4;
  PetscInt       Ilast,i,j,II,n;
  PetscScalar    *initial_phi1,*initial_phi2,*initial_phi3,*initial_phi4;
  PetscScalar    *initial_THP,*initial_THTf,*initial_THTs,*initial_THU,*initial_THV;

  PetscFunctionBeginUser;
  ierr =VecGetArray(iniphi1,&initial_phi1);
  ierr =VecGetArray(iniphi2,&initial_phi2);
  ierr =VecGetArray(iniphi3,&initial_phi3);
  ierr =VecGetArray(iniphi4,&initial_phi4);
  ierr =VecGetArray(u21,&initial_THP);
  ierr =VecGetArray(u22,&initial_THTf);
  ierr =VecGetArray(u27,&initial_THTs);
  ierr =VecGetArray(u24,&initial_THU);
  ierr =VecGetArray(u25,&initial_THV);

  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;
  THInterval1 = 1672;
  THInterval2 = 3344;
  THInterval3 = 5016;
  THInterval4 = 6688;
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
  for (j=1; j<76+1; j++){
    for (i=1; i<22+1; i++){
        II = (j-1)*22 + i - 1  ;
        u[II + 18721]             =  initial_THP[II];
        u[II+THInterval1 + 18721] =  initial_THTf[II];
        u[II+THInterval2 + 18721] =  initial_THU[II];
        u[II+THInterval3 + 18721] =  initial_THV[II];
        
    }
  }
  
  for (j=1; j<104+1; j++){
    for (i=1; i<45+1; i++){
        II = (j-1)*45 + i - 1 ;
        u[II+THInterval4 + 18721] =  initial_THTs[II];    
    }
  }


  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi1,&initial_phi1);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi2,&initial_phi2);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi3,&initial_phi3);CHKERRQ(ierr);
  ierr = VecRestoreArray(iniphi4,&initial_phi4);CHKERRQ(ierr);
  ierr = VecRestoreArray(u21,&initial_THP);CHKERRQ(ierr);
  ierr = VecRestoreArray(u22,&initial_THTf);CHKERRQ(ierr);
  ierr = VecRestoreArray(u27,&initial_THTs);CHKERRQ(ierr);
  ierr = VecRestoreArray(u24,&initial_THU);CHKERRQ(ierr);
  ierr = VecRestoreArray(u25,&initial_THV);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  THRP2          *user = (THRP2*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,II,n,nth,ith;
  PetscInt       ys,ym,j;/*y direction index*/
  PetscInt       group;
  PetscInt       indexc,indexl,indexr,indexu,indexd;
  PetscInt       start,end,Interval1,Interval2,Interval3;
  PetscInt       THInterval1,THInterval2,THInterval3,THInterval4;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  const PetscScalar    *u;
  PetscScalar    *fu,*Ts,*Qf0;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,cc1,cc2,cc3,cc4,ac1,ac2,ac3,ac4,ac0;
  PetscScalar    phigroup1,phigroup2,phigroup3,phigroup4;
  PetscScalar    keff;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(FU,&fu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  n     = 45;
  xs = 0; ys = 0; xm = 45; ym = 104;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;
  nth     = 22;
  ith     = 18721;
  THInterval1 = 1672;
  THInterval2 = 3344;
  THInterval3 = 5016;
  THInterval4 = 6688;

  /* abstract RP and TH variables */
  //user->vecGetTHRPArray(u);
  /*abstract Ts */
  PetscMalloc1(Interval1,&Ts);
  for (i = 0; i < Interval1; i++)
  {
    Ts[i] = 0.0;
  }


  user->updateVoidCrossSection(Ts);


  /* Compute c and Qfth */
  user->computeQfth(u);
  PetscMalloc1(Interval1,&Qf0);
  for (i = 0; i < Interval1; i++)
  {
    Qf0[i] = user->Qf[i];
  }

  user->setQf(Qf0);

  /* write the current variables */
  user->setCurrentVariables(u);

  /*set parameter*/
  user->setParameterAll();

  /* set boundary coefficient */
  user->setBoundaryCoefficientAll();


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
      rightside = user->getRPrs(i,j,group,phigroup1,phigroup2,phigroup3,phigroup4,keff) ;
      if ( i > 1 )
      {
        Ileft = II - 1 ;
        indexl = Ileft + (group-1)*45*104;
        ac1 = user->getRPa1(i,j,group);
        cc1 = ac1 * u[indexl];
      }
      if ( i < 45 )
      {
        Iright = II + 1 ;
        indexr = Iright + (group-1)*45*104;
        ac2 = user->getRPa2(i,j,group);
        cc2 = ac2 * u[indexr];
      }
      if ( j > 1 )
      {
        Idown = II - 45 ;
        indexd = Idown + (group-1)*45*104 ;
        ac3 = user->getRPa3(i,j,group);
        cc3 = ac3 * u[indexd];
      }
      if ( j < 104 )
      {
        Iup = II + 45 ;
        indexu = Iup + (group-1)*45*104 ;
        ac4 = user->getRPa4(i,j,group);
        cc4 = ac4 * u[indexu];
      }
      indexc = (group-1)*104*45 + II;
      ac1 = user->getRPa1(i,j,group);
      ac2 = user->getRPa2(i,j,group);
      ac3 = user->getRPa3(i,j,group);
      ac4 = user->getRPa4(i,j,group);
      ac0 = user->getRPa0(i,j,group,ac1,ac2,ac3,ac4);
      fu[indexc] = ac0 * u[indexc] +cc1+cc2+cc3+cc4-rightside;
      }  
    }
  }


  /* inner points */
  for (j=12; j<87+1; j++){
    for (i=1; i<22+1; i++){

      /* Tf field */ 
      II = (j-12)*nth + i - 1 + ith ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTfa1(i,j);
        cc1 = ac1 * u[Ileft+THInterval1];
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getTfa2(i,j);
        cc2 = ac2 * u[Iright+THInterval1];
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getTfa3(i,j);
        cc3 = ac3 * u[Idown+THInterval1];
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getTfa4(i,j);
        cc4 = ac4 * u[Iup+THInterval1];
      }
      ac1 = user->getTfa1(i,j);
      ac2 = user->getTfa2(i,j);
      ac3 = user->getTfa3(i,j);
      ac4 = user->getTfa4(i,j);
      ac0 = user->getTfa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getTfrs(i,j,ac3);
      fu[II+THInterval1] = ac0 * u[II+THInterval1] +cc1+cc2+cc3+cc4-rightside;

      /* U field */
      II = (j-12)*nth + i - 1 + ith ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      ac0 = 1.0;
      rightside = user->getUrs(i,j);
      fu[II+THInterval2] = ac0 * u[II+THInterval2] +cc1+cc2+cc3+cc4-rightside;

      /* V field */
      II = (j-12)*nth + i - 1 + ith ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      ac0 = 1.0;
      rightside = user->getVrs(i,j);
      fu[II+THInterval3] = ac0 * u[II+THInterval3] +cc1+cc2+cc3+cc4-rightside;

      /* P field */
      II = (j-12)*nth + i - 1 + ith ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getPa1(i,j);
        cc1 = ac1 * u[Ileft];
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getPa2(i,j);
        cc2 = ac2 * u[Iright];
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getPa3(i,j);
        cc3 = ac3 * u[Idown];
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getPa4(i,j);
        cc4 = ac4 * u[Iup];
      }
      ac1 = user->getPa1(i,j);
      ac2 = user->getPa2(i,j);
      ac3 = user->getPa3(i,j);
      ac4 = user->getPa4(i,j);
      ac0 = user->getPa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getPrs(i,j,ac4);
      fu[II] = ac0 * u[II] +cc1+cc2+cc3+cc4-rightside;
    }
  }

  for (j = 1; j < 104+1; j++)
  {
    for (i = 1; i < 45+1; i++)
    {
      /*Ts field*/
      II = (j-1)*45 + i - 1 + ith ;
      cc1 = 0.0;
      cc2 = 0.0;
      cc3 = 0.0;
      cc4 = 0.0;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTsa1(i,j);
        cc1 = ac1 * u[Ileft+THInterval4];
      }
      if ( i < 45 )
      {
        Iright = II + 1 ;
        ac2 = user->getTsa2(i,j);
        cc2 = ac2 * u[Iright+THInterval4];
      }
      if ( j > 1 )
      {
        Idown = II - 45 ;
        ac3 = user->getTsa3(i,j);
        cc3 = ac3 * u[Idown+THInterval4];
      }
      if ( j < 104 )
      {
        Iup = II + 45 ;
        ac4 = user->getTsa4(i,j);
        cc4 = ac4 * u[Iup+THInterval4];
      }
      ac1 = user->getTsa1(i,j);
      ac2 = user->getTsa2(i,j);
      ac3 = user->getTsa3(i,j);
      ac4 = user->getTsa4(i,j);
      ac0 = user->getTsa0(i,j,ac1,ac2,ac3,ac4);
      rightside = user->getTsrs(i,j,ac2,ac3,ac4);
      fu[II+THInterval4] = ac0 * u[II+THInterval4] +cc1+cc2+cc3+cc4-rightside;
    }
  }

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(FU,&fu);CHKERRQ(ierr);
 
  PetscLogFlops(xm*ym*4.0*27.0);

  PetscFunctionReturn(0);
}

PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat B,void *dummy)
{
  THRP2          *user = (THRP2*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,II,n,nth,ith;
  PetscInt       ys,ym,j;/*y direction index*/
  PetscInt       group;
  PetscInt       start,end,Interval1,Interval2,Interval3;
  PetscInt       THInterval1,THInterval2,THInterval3,THInterval4;
  PetscInt       Ileft,Iright,Idown,Iup,Ilast;
  PetscInt       idcenter,idright,idleft,idup,iddown;
  const PetscScalar    *u;
  PetscScalar    *Ts;
  PetscScalar    rightside,e,b,c,dd,a;
  PetscScalar    ac0,ac1,ac2,ac3,ac4;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  n         = 45;
  xs = 0; ys = 0; xm = 45; ym = 104;
  Interval1 = 4680;
  Interval2 = 9360;
  Interval3 = 14040;
  Ilast     = 18720;
  nth       = 22;
  ith       = 18721;
  THInterval1 = 1672;
  THInterval2 = 3344;
  THInterval3 = 5016;
  THInterval4 = 6688;

  /* abstract RP and TH variables */
  //user->vecGetTHRPArray(u);
  PetscMalloc1(Interval1,&Ts);
  for (i = 0; i < Interval1; i++)
  {
    Ts[i] = 0.0;
  }
  user->updateVoidCrossSection(Ts);

  /* write the current variables */
  user->setCurrentVariables(u);

  /*set parameter*/
  user->setParameterAll();

  /* set boundary coefficient */
  user->setBoundaryCoefficientAll();

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
        b  = user->getRPa1(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idleft,&b,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 45 )
      {
        idright  = II +  1;
        c  = user->getRPa2(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idright,&c,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 1 )
      {
        iddown   = II - 45;
        dd = user->getRPa3(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&iddown,&dd,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 104 )
      {
        idup     = II + 45;
        e  = user->getRPa4(i,j,group) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idup,&e,INSERT_VALUES);CHKERRQ(ierr);
      }        
        a  = user->getRPa0(i,j,group,b,c,dd,e) ;
        ierr  = MatSetValues(B,1,&idcenter,1,&idcenter,&a,INSERT_VALUES);CHKERRQ(ierr);
      }
    }
  }


  /* inner points */
  for (j=12; j<87+1; j++){
    for (i=1; i<22+1; i++){
      
      /* Tf field */ 
      II = (j-12)*nth + i - 1 + THInterval1 +ith ;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTfa1(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Ileft,&ac1,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getTfa2(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iright,&ac2,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getTfa3(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Idown,&ac3,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getTfa4(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iup,&ac4,INSERT_VALUES);CHKERRQ(ierr);
      }
      ac1 = user->getTfa1(i,j);
      ac2 = user->getTfa2(i,j);
      ac3 = user->getTfa3(i,j);
      ac4 = user->getTfa4(i,j);
      ac0 = user->getTfa0(i,j,ac1,ac2,ac3,ac4);
      ierr  = MatSetValues(B,1,&II,1,&II,&ac0,INSERT_VALUES);CHKERRQ(ierr);

      /* U field */
      II = (j-12)*nth + i - 1 + THInterval2 + ith ;
      ac0 = 1.0;
      ierr  = MatSetValues(B,1,&II,1,&II,&ac0,INSERT_VALUES);CHKERRQ(ierr);

      /* V field */
      II = (j-12)*nth + i - 1 + THInterval3 + ith ;
      ac0 = 1.0;
      ierr  = MatSetValues(B,1,&II,1,&II,&ac0,INSERT_VALUES);CHKERRQ(ierr);

      /* P field */
      II = (j-12)*nth + i - 1 + ith;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getPa1(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Ileft,&ac1,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 22 )
      {
        Iright = II + 1 ;
        ac2 = user->getPa2(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iright,&ac2,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 12 )
      {
        Idown = II - 22 ;
        ac3 = user->getPa3(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Idown,&ac3,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 87 )
      {
        Iup = II + 22 ;
        ac4 = user->getPa4(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iup,&ac4,INSERT_VALUES);CHKERRQ(ierr);
      }
      ac1 = user->getPa1(i,j);
      ac2 = user->getPa2(i,j);
      ac3 = user->getPa3(i,j);
      ac4 = user->getPa4(i,j);
      ac0 = user->getPa0(i,j,ac1,ac2,ac3,ac4);
      ierr  = MatSetValues(B,1,&II,1,&II,&ac0,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  for (j = 1; j < 104+1; j++)
  {
    for (i = 1; i < 45+1; i++)
    {
      /*Ts field*/
      II = (j-1)*45 + i - 1 + THInterval4 + ith ;
      if ( i > 1 )
      {
        Ileft = II - 1;
        ac1 = user->getTsa1(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Ileft,&ac1,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( i < 45 )
      {
        Iright = II + 1 ;
        ac2 = user->getTsa2(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iright,&ac2,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j > 1 )
      {
        Idown = II - 45 ;
        ac3 = user->getTsa3(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Idown,&ac3,INSERT_VALUES);CHKERRQ(ierr);
      }
      if ( j < 104 )
      {
        Iup = II + 45 ;
        ac4 = user->getTsa4(i,j);
        ierr  = MatSetValues(B,1,&II,1,&Iup,&ac4,INSERT_VALUES);CHKERRQ(ierr);
      }
      ac1 = user->getTsa1(i,j);
      ac2 = user->getTsa2(i,j);
      ac3 = user->getTsa3(i,j);
      ac4 = user->getTsa4(i,j);
      ac0 = user->getTsa0(i,j,ac1,ac2,ac3,ac4);
      ierr  = MatSetValues(B,1,&II,1,&II,&ac0,INSERT_VALUES);CHKERRQ(ierr);
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
  THRP2           *user = (THRP2*)dummy;
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
