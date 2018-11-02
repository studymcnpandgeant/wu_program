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
#include "Grid.h"

/*

       w - design variables (what we change to get an optimal solution)
       u - state variables (i.e. the PDE solution)
       lambda - the Lagrange multipliers

            U = (w u lambda)

       fu, fw, flambda contain the gradient of L(w,u,lambda)

            FU = (fw fu flambda)

       In this example the PDE is
                             Uxx = 2
                            u(0) = w(0), thus this is the free parameter
                            u(1) = 0
       the function we wish to minimize is
                            \integral u^{2}

       The exact solution for u is given by u(x) = x*x - 1.25*x + .25

       Use the usual centered finite differences.

       Note we treat the problem as non-linear though it happens to be linear

       See ex22.c for the same code, but that interlaces the u and the lambda

*/

// typedef struct {
//   DM          red1,da1,da2,da3,da4;
//   DM          packer;
//   //PetscViewer u_viewer,lambda_viewer;
//   //PetscViewer fu_viewer,flambda_viewer;
// /*  cross section coefficients  */
//   PetscReal   sigmaf_core1_1,sigmaf_core1_2,sigmaf_core2_1,sigmaf_core2_2;
//   PetscReal   sigmaa_core1_1,sigmaa_core1_2,sigmaa_core2_1,sigmaa_core2_2,sigmaa_re_1,sigmaa_re_2;
//   PetscReal   sigmas_core1_1,sigmas_core2_1,sigmas_re_1;
//   PetscReal   d_core1_1,d_core1_2,d_core2_1,d_core2_2,d_re_1,d_re_2;
// } UserCtx;

typedef struct
{
  PetscScalar *corearray,*materialarray;
  PetscScalar *rcor,*zcor;
  PetscScalar *rrcor,*zzcor;
} AppCtx;

extern PetscErrorCode FormInitialGuess(Grid*,Vec,Mat,Mat,Mat,Mat);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  AppCtx         data;
  PetscErrorCode ierr;
  PetscInt       its;
  Vec            U,FU,vlambda,vphi1,vphi2,vphi3,vphi4;
  Vec            u1,u2,u3,u4,u5,u6;
  Mat            A1,A2,A3,A4;/*Matrix for storage initial guess, four groups*/
  SNES           snes;
  Grid           user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Reading vectors in binary format */
  /* Read four vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-sn2dtest5",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u4);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u5);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u6);CHKERRQ(ierr);
  ierr = VecLoad(u1,viewer);CHKERRQ(ierr);/* materials data cross-section infomation */
  ierr = VecLoad(u2,viewer);CHKERRQ(ierr);/* grid's materials infomation (INT format) */
  ierr = VecLoad(u3,viewer);CHKERRQ(ierr);/* dx-infomation */
  ierr = VecLoad(u4,viewer);CHKERRQ(ierr);/* dy-infomation */
  ierr = VecLoad(u5,viewer);CHKERRQ(ierr);/* x-infomation */
  ierr = VecLoad(u6,viewer);CHKERRQ(ierr);/* y-infomation */
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the arrays in user's appctx */
  ierr =VecGetArray(u1,&data.materialarray);
  ierr =VecGetArray(u2,&data.corearray);
  ierr =VecGetArray(u3,&data.rcor);
  ierr =VecGetArray(u4,&data.zcor);
  ierr =VecGetArray(u5,&data.rrcor);
  ierr =VecGetArray(u6,&data.zzcor);

  user.setinitial(data.materialarray,data.corearray,data.rcor,data.zcor,data.rrcor,data.zzcor);

  /*Read the file and store them in the matrix*/
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg)
  {
    SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f option");
  }
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A1);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A1);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A2);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A2);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A3);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A3);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&A4);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A4);CHKERRQ(ierr);
  ierr = MatLoad(A1,viewer);CHKERRQ(ierr);
  ierr = MatLoad(A2,viewer);CHKERRQ(ierr);
  ierr = MatLoad(A3,viewer);CHKERRQ(ierr);
  ierr = MatLoad(A4,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);




  /* Create a global vector that includes a single redundant array and two da arrays */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&user.packer);CHKERRQ(ierr);
  ierr = DMRedundantCreate(PETSC_COMM_WORLD,0,1,&user.red1);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.packer,user.red1);CHKERRQ(ierr);
  /*  Set neutron flux DM   */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,16,35,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da1);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da1);CHKERRQ(ierr);
  ierr = DMSetUp(user.da1);CHKERRQ(ierr);
  /* Set composite DM object */
  ierr = DMCompositeAddDM(user.packer,user.da1);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,16,35,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da2);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da2);CHKERRQ(ierr);
  ierr = DMSetUp(user.da2);CHKERRQ(ierr);
  /* Set composite DM object */
  ierr = DMCompositeAddDM(user.packer,user.da2);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,16,35,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da3);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da3);CHKERRQ(ierr);
  ierr = DMSetUp(user.da3);CHKERRQ(ierr);
  /* Set composite DM object */
  ierr = DMCompositeAddDM(user.packer,user.da3);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,16,35,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da4);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da4);CHKERRQ(ierr);
  ierr = DMSetUp(user.da4);CHKERRQ(ierr);
  /* Set composite DM object */
  ierr = DMCompositeAddDM(user.packer,user.da4);CHKERRQ(ierr);
  
  
  ierr = DMDASetFieldName(user.da1,0,"phi1");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da2,0,"phi2");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da2,0,"phi3");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da2,0,"phi4");CHKERRQ(ierr);
  //ierr = DMDASetFieldName(user.red1,0,"lambda");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user.packer,&U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(&user,U,A1,A2,A3,A4);CHKERRQ(ierr);

  /* create nonlinear solver */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes,Monitor,&user,0);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  /* Output the data */
  ierr = VecCreate(PETSC_COMM_WORLD,&vlambda);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&vphi1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&vphi2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&vphi3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&vphi4);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_testfourgroups",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user.packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user.packer,U,vlambda,vphi1,vphi2,vphi3,vphi4);CHKERRQ(ierr);
  ierr = VecView(vlambda,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi1,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi2,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi3,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi4,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = DMDestroy(&user.red1);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da1);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da2);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da3);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da4);CHKERRQ(ierr);
  ierr = DMDestroy(&user.packer);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&FU);CHKERRQ(ierr);
  ierr = VecDestroy(&vlambda);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi1);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi2);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi3);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi4);CHKERRQ(ierr);
  ierr = VecDestroy(&u1);CHKERRQ(ierr);
  ierr = VecDestroy(&u2);CHKERRQ(ierr);
  ierr = VecDestroy(&u3);CHKERRQ(ierr);
  ierr = VecDestroy(&u4);CHKERRQ(ierr);
  ierr = VecDestroy(&u5);CHKERRQ(ierr);
  ierr = VecDestroy(&u6);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


PetscErrorCode FormInitialGuess(Grid *user,Vec U,Mat A1,Mat A2,Mat A3,Mat A4)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       id;
  PetscScalar    *lambda,**phi1,**phi2,**phi3,**phi4;
  PetscScalar    *initial_phi1,*initial_phi2,*initial_phi3,*initial_phi4;
  Vec            vlambda,vphi1,vphi2,vphi3,vphi4;

  PetscFunctionBeginUser;
  ierr = MatSeqAIJGetArray(A1,&initial_phi1);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(A2,&initial_phi2);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(A3,&initial_phi3);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(A4,&initial_phi4);CHKERRQ(ierr);

  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2,vphi3,vphi4);CHKERRQ(ierr);

  ierr = DMDAGetCorners(user->da1,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);//get two direction index
  ierr = DMDAGetInfo(user->da1,0,&N,&M,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da3,vphi3,&phi3);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da4,vphi4,&phi4);CHKERRQ(ierr);

  if (xs == 0 && ys == 0) lambda[0] = 1.0;
  for ( j=ys ; j < ys+ym ; j++ ){
       for ( i=xs ; i < xs+xm ; i++ ){
          id = j*xm + i ;
          phi1[j][i] = initial_phi1[id] ;
          phi2[j][i] = initial_phi2[id] ;
          phi3[j][i] = initial_phi3[id] ;
          phi4[j][i] = initial_phi4[id] ;
       }
    }

  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da3,vphi3,&phi3);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da4,vphi4,&phi4);CHKERRQ(ierr);

  ierr = DMCompositeGather(user->packer,INSERT_VALUES,U,vlambda,vphi1,vphi2,vphi3,vphi4);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  Grid           *user = (Grid*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N,II,n;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscInt       group1=1,group2=2,group3=3,group4=4;
  PetscInt       chigroup1=1,chigroup2=2,chigroup3=3;
  PetscScalar    *lambda,**phi1,**phi2,**phi3,**phi4,*flambda,**fphi1,**fphi2,**fphi3,**fphi4;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,e,b,c,dd,a;
  PetscReal      sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;
  PetscReal      sigmas12,sigmas13,sigmas14,sigmas23,sigmas24,sigmas34;
  PetscReal      chisigmaf1,chisigmaf2,chisigmaf3,chisigmaf4;
  PetscReal      phi1boundary,phi2boundary,phi3boundary,phi4boundary;
  Vec            vlambda,vphi1,vphi2,vphi3,vphi4,vflambda,vfphi1,vfphi2,vfphi3,vfphi4;

  PetscFunctionBeginUser;
  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->packer,&vflambda,&vfphi1,&vfphi2,&vfphi3,&vfphi4);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2,vphi3,vphi4);CHKERRQ(ierr);

  ierr = DMDAGetCorners(user->da1,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);//get two direction index
  ierr = DMDAGetInfo(user->da1,0,&N,&M,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = VecGetArray(vflambda,&flambda);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vfphi1,&fphi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vfphi2,&fphi2);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da3,vphi3,&phi3);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da3,vfphi3,&fphi3);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da4,vphi4,&phi4);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da4,vfphi4,&fphi4);CHKERRQ(ierr);
  n     = 16;
  phi1boundary = 0.0;
  phi2boundary = 0.0;
  phi3boundary = 0.0;
  phi4boundary = 0.0;


  /* residual f_lambda */
  if (xs == 0 && ys == 0 ) { /* only first processor computes this */
    dotphi = 0.0;
    for ( j=ys ; j < ys+ym ; j++ ){
       for ( i=xs ; i < xs+xm ; i++ ){
          dotphi += phi1[j][i]*phi1[j][i] + phi2[j][i]*phi2[j][i] + phi3[j][i]*phi3[j][i] + phi4[j][i]*phi4[j][i] ;
       }
    }
    flambda[0] =  dotphi - 1.0 ;
  }

  /* bottom boundary */

  if (ys == 0) {
    j     = 0;
    yints = ys + 1;
    /* bottom edge */
    for (i=xs; i<xs+xm; i++) {
      /*phi1 field*/
      rightside = phi1boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;

      /*phi3 field*/
      rightside = phi3boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi3[j][i] = a * phi3[j][i] - rightside;

      /*phi4 field*/
      rightside = phi4boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi4[j][i] = a * phi4[j][i] - rightside;
    }
  }

  /* top boundary */
  if (yints == 1) {
    j     = ys + ym - 1;
    yinte = ys + ym - 1;
    for (i=xs; i<xs+xm; i++) {
      /*phi1 field*/
      rightside = phi1boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;

      /*phi3 field*/
      rightside = phi3boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi3[j][i] = a * phi3[j][i] - rightside;

      /*phi4 field*/
      rightside = phi4boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi4[j][i] = a * phi4[j][i] - rightside;
    }
  }

  /* left boundary */

  if (xs == 0) {
    i     = 0;
    xints = xs + 1;
    /* bottom edge */
    for (j=yints; j<yinte; j++) {
      II = j*n + i ;
      /*phi1 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup1,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup1,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup1,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup1,group4);
      rightside = lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = 0.0 ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group1) ;
      dd = user->getrzas(j,i,II,group1) ;
      e  = user->getrzan(j,i,II,group1) ;
      sigmat1  = user->getrzsigmat(i,II,group1);
      sigmas11 = user->getrzsigmas(i,II,group1,group1);
      a  = (sigmat1 - sigmas11) - b - c - dd - e ;
      fphi1[j][i] = a* phi1[j][i] + c* phi1[j][i+1] + dd* phi1[j-1][i] + e* phi1[j+1][i] - rightside ;

      /*phi2 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup2,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup2,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup2,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup2,group4);
      sigmas12   = user->getrzsigmas(i,II,group1,group2);
      rightside = sigmas12 *phi1[j][i] + lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = 0.0 ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group2) ;
      dd = user->getrzas(j,i,II,group2) ;
      e  = user->getrzan(j,i,II,group2) ;
      sigmat2  = user->getrzsigmat(i,II,group2);
      sigmas22 = user->getrzsigmas(i,II,group2,group2);
      a  = (sigmat2 - sigmas22) - b - c - dd - e ;
      fphi2[j][i] = a* phi2[j][i] + c* phi2[j][i+1] + dd* phi2[j-1][i] + e* phi2[j+1][i] - rightside ;

      /*phi3 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup3,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup3,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup3,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup3,group4);
      sigmas13   = user->getrzsigmas(i,II,group1,group3);
      sigmas23   = user->getrzsigmas(i,II,group2,group3);
      rightside = sigmas13 *phi1[j][i] + sigmas23 *phi2[j][i] + lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = 0.0 ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group3) ;
      dd = user->getrzas(j,i,II,group3) ;
      e  = user->getrzan(j,i,II,group3) ;
      sigmat3  = user->getrzsigmat(i,II,group3);
      sigmas33 = user->getrzsigmas(i,II,group3,group3);
      a  = (sigmat3 - sigmas33) - b - c - dd - e ;
      fphi3[j][i] = a* phi3[j][i] + c* phi3[j][i+1] + dd* phi3[j-1][i] + e* phi3[j+1][i] - rightside ;

      /*phi4 field*/
      sigmas14 = user->getrzsigmas(i,II,group1,group4);
      sigmas24 = user->getrzsigmas(i,II,group2,group4);
      sigmas34 = user->getrzsigmas(i,II,group3,group4);
      rightside = sigmas14 *phi1[j][i] + sigmas24 *phi2[j][i] + sigmas34 *phi3[j][i] ;
      b  = 0.0 ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group4) ;
      dd = user->getrzas(j,i,II,group4) ;
      e  = user->getrzan(j,i,II,group4) ;
      sigmat4  = user->getrzsigmat(i,II,group4);
      sigmas44 = user->getrzsigmas(i,II,group4,group4);
      a  = (sigmat4 - sigmas44) - b - c - dd - e ;
      fphi4[j][i] = a* phi4[j][i] + c* phi4[j][i+1] + dd* phi4[j-1][i] + e* phi4[j+1][i] - rightside ;

    }
  }

  /* right boundary */
  if (xints == 1) {
    i     = xs + xm - 1;
    xinte = xs + xm - 1;
    for (j=yints; j<yinte; j++) {
      /*phi1 field*/
      rightside = phi1boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;

      /*phi3 field*/
      rightside = phi3boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi3[j][i] = a * phi3[j][i] - rightside;

      /*phi4 field*/
      rightside = phi4boundary;
      b         = 0.0;
      c         = 0.0;
      dd        = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi4[j][i] = a * phi4[j][i] - rightside;
    }
  }

  /* inner points */
  for (j=yints; j<yinte; j++){
    for (i=xints; i<xinte; i++){
      II = j*n + i ;
      /*phi1 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup1,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup1,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup1,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup1,group4);
      rightside = lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = user->getrzaw(j,i,II,group1) ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group1) ;
      dd = user->getrzas(j,i,II,group1) ;
      e  = user->getrzan(j,i,II,group1) ;
      sigmat1  = user->getrzsigmat(i,II,group1);
      sigmas11 = user->getrzsigmas(i,II,group1,group1);
      a  = (sigmat1 - sigmas11) - b - c - dd - e ;
      fphi1[j][i] = a* phi1[j][i] + b* phi1[j][i-1] + c* phi1[j][i+1] + dd* phi1[j-1][i] + e* phi1[j+1][i] - rightside ;

      /*phi2 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup2,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup2,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup2,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup2,group4);
      sigmas12   = user->getrzsigmas(i,II,group1,group2);
      rightside = sigmas12 *phi1[j][i] + lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = user->getrzaw(j,i,II,group2) ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group2) ;
      dd = user->getrzas(j,i,II,group2) ;
      e  = user->getrzan(j,i,II,group2) ;
      sigmat2  = user->getrzsigmat(i,II,group2);
      sigmas22 = user->getrzsigmas(i,II,group2,group2);
      a  = (sigmat2 - sigmas22) - b - c - dd - e ;
      fphi2[j][i] = a* phi2[j][i] + b* phi2[j][i-1] + c* phi2[j][i+1] + dd* phi2[j-1][i] + e* phi2[j+1][i] - rightside ;

      /*phi3 field*/
      chisigmaf1 = user->getrzchisigmaf(i,II,chigroup3,group1);
      chisigmaf2 = user->getrzchisigmaf(i,II,chigroup3,group2);
      chisigmaf3 = user->getrzchisigmaf(i,II,chigroup3,group3);
      chisigmaf4 = user->getrzchisigmaf(i,II,chigroup3,group4);
      sigmas13   = user->getrzsigmas(i,II,group1,group3);
      sigmas23   = user->getrzsigmas(i,II,group2,group3);
      rightside = sigmas13 *phi1[j][i] + sigmas23 *phi2[j][i] + lambda[0] * ( chisigmaf1 * phi1[j][i] + chisigmaf2 * phi2[j][i] + chisigmaf3 * phi3[j][i] + chisigmaf4 * phi4[j][i] ) ;
      b  = user->getrzaw(j,i,II,group3) ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group3) ;
      dd = user->getrzas(j,i,II,group3) ;
      e  = user->getrzan(j,i,II,group3) ;
      sigmat3  = user->getrzsigmat(i,II,group3);
      sigmas33 = user->getrzsigmas(i,II,group3,group3);
      a  = (sigmat3 - sigmas33) - b - c - dd - e ;
      fphi3[j][i] = a* phi3[j][i] + b* phi3[j][i-1] + c* phi3[j][i+1] + dd* phi3[j-1][i] + e* phi3[j+1][i] - rightside ;

      /*phi4 field*/
      sigmas14 = user->getrzsigmas(i,II,group1,group4);
      sigmas24 = user->getrzsigmas(i,II,group2,group4);
      sigmas34 = user->getrzsigmas(i,II,group3,group4);
      rightside = sigmas14 *phi1[j][i] + sigmas24 *phi2[j][i] + sigmas34 *phi3[j][i] ;
      b  = user->getrzaw(j,i,II,group4) ;// index opposite of examples in SLEPc
      c  = user->getrzae(j,i,II,group4) ;
      dd = user->getrzas(j,i,II,group4) ;
      e  = user->getrzan(j,i,II,group4) ;
      sigmat4  = user->getrzsigmat(i,II,group4);
      sigmas44 = user->getrzsigmas(i,II,group4,group4);
      a  = (sigmat4 - sigmas44) - b - c - dd - e ;
      fphi4[j][i] = a* phi4[j][i] + b* phi4[j][i-1] + c* phi4[j][i+1] + dd* phi4[j-1][i] + e* phi4[j+1][i] - rightside ;

    }
  }


  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = VecRestoreArray(vflambda,&flambda);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vfphi1,&fphi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vfphi2,&fphi2);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da3,vphi3,&phi3);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da3,vfphi3,&fphi3);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da4,vphi4,&phi4);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da4,vfphi4,&fphi4);CHKERRQ(ierr);

  ierr = DMCompositeGather(user->packer,INSERT_VALUES,FU,vflambda,vfphi1,vfphi2,vfphi3,vfphi4);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vflambda,&vfphi1,&vfphi2,&vfphi3,&vfphi4);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *dummy)
{
  Grid           *user = (Grid*)dummy;
  PetscErrorCode ierr;
  Vec            vlambda,vphi1,vphi2,vphi3,vphi4,U;
  PetscReal      *lambda;
  PetscReal      keff;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2,&vphi3,&vphi4);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2,vphi3,vphi4);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);

  keff = 1.0 / lambda[0] ;

  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, Keff ======= %g\n",its,(double)keff);CHKERRQ(ierr);
  //ierr = VecView(vphi1,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
