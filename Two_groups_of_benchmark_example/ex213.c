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

typedef struct {
  DM          red1,da1,da2;
  DM          packer;
  //PetscViewer u_viewer,lambda_viewer;
  //PetscViewer fu_viewer,flambda_viewer;
/*  cross section coefficients  */
  PetscReal   sigmaf_core1_1,sigmaf_core1_2,sigmaf_core2_1,sigmaf_core2_2;
  PetscReal   sigmaa_core1_1,sigmaa_core1_2,sigmaa_core2_1,sigmaa_core2_2,sigmaa_re_1,sigmaa_re_2;
  PetscReal   sigmas_core1_1,sigmas_core2_1,sigmas_re_1;
  PetscReal   d_core1_1,d_core1_2,d_core2_1,d_core2_2,d_re_1,d_re_2;
} UserCtx;

extern PetscErrorCode FormInitialGuess(UserCtx*,Vec,Mat,Mat);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       its;
  Vec            U,FU,vlambda,vphi1,vphi2;
  Mat            A1,A2;/*Matrix for storage initial guess*/
  SNES           snes;
  UserCtx        user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* set cross-section coefficient */
  user.sigmaf_core1_1     = 0.0085;
  user.sigmaf_core1_2     = 0.1851;
  user.sigmaf_core2_1     = 0.006;
  user.sigmaf_core2_2     = 0.150;
  user.sigmaa_core1_1     = 0.0121;
  user.sigmaa_core1_2     = 0.121;
  user.sigmaa_core2_1     = 0.01;
  user.sigmaa_core2_2     = 0.1;
  user.sigmaa_re_1        = 0.0004;
  user.sigmaa_re_2        = 0.020;
  user.sigmas_core1_1     = 0.0241;
  user.sigmas_core2_1     = 0.016;
  user.sigmas_re_1        = 0.0493;
  user.d_core1_1          = 1.267;
  user.d_core1_2          = 0.354;
  user.d_core2_1          = 1.280;
  user.d_core2_2          = 0.400;
  user.d_re_1             = 1.130;
  user.d_re_2             = 0.166;
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
  ierr = MatLoad(A1,viewer);CHKERRQ(ierr);
  ierr = MatLoad(A2,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);




  /* Create a global vector that includes a single redundant array and two da arrays */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&user.packer);CHKERRQ(ierr);
  ierr = DMRedundantCreate(PETSC_COMM_WORLD,0,1,&user.red1);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.packer,user.red1);CHKERRQ(ierr);
/*  Set neutron flux DM   */
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,20,24,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da1);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da1);CHKERRQ(ierr);
  ierr = DMSetUp(user.da1);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.packer,user.da1);CHKERRQ(ierr);
  ierr = DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,20,24,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,&user.da2);CHKERRQ(ierr);
  ierr = DMSetFromOptions(user.da2);CHKERRQ(ierr);
  ierr = DMSetUp(user.da2);CHKERRQ(ierr);
  ierr = DMCompositeAddDM(user.packer,user.da2);CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da1,0,"phi1");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user.da2,0,"phi2");CHKERRQ(ierr);
  //ierr = DMDASetFieldName(user.red1,0,"lambda");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user.packer,&U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(&user,U,A1,A2);CHKERRQ(ierr);

  /* create nonlinear solver */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);
  ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes,Monitor,&user,0);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  /* Output the data */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_testtwogroups",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user.packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user.packer,U,vlambda,vphi1,vphi2);CHKERRQ(ierr);
  ierr = VecView(vlambda,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi1,viewer);CHKERRQ(ierr);
  ierr = VecView(vphi2,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = DMDestroy(&user.red1);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da1);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da2);CHKERRQ(ierr);
  ierr = DMDestroy(&user.packer);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&FU);CHKERRQ(ierr);
  ierr = VecDestroy(&vlambda);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi1);CHKERRQ(ierr);
  ierr = VecDestroy(&vphi2);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

PETSC_STATIC_INLINE PetscReal d(const UserCtx *ctx,PetscInt j,PetscInt i,PetscInt g)
{
  //if ( i>2 && i<17 && j>2 && j<11 ) return g % 2 ? ctx->d_core1_1 : ctx->d_core1_2;
  //else if ( i>2 && i<17 && j>10 && j<21 ) return g % 2 ? ctx->d_core2_1 : ctx->d_core2_2;
  return g % 2 ? ctx->d_core1_1 : ctx->d_core1_2 ;
}

PETSC_STATIC_INLINE PetscReal sigmaa(const UserCtx *ctx,PetscInt j,PetscInt i,PetscInt g)
{
  //if ( i>2 && i<17 && j>2 && j<11 ) return g % 2 ? ctx->sigmaa_core1_1 : ctx->sigmaa_core1_2;
  //else if ( i>2 && i<17 && j>10 && j<21 ) return g % 2 ? ctx->sigmaa_core2_1 : ctx->sigmaa_core2_2;
  return g % 2 ? ctx->sigmaa_core1_1 : ctx->sigmaa_core1_2 ;
}

PETSC_STATIC_INLINE PetscReal sigmaf(const UserCtx *ctx,PetscInt j,PetscInt i,PetscInt g)
{
  //if ( i>2 && i<17 && j>2 && j<11 ) return g % 2 ? ctx->sigmaf_core1_1 : ctx->sigmaf_core1_2;
  //else if ( i>2 && i<17 && j>10 && j<21 ) return g % 2 ? ctx->sigmaf_core2_1 : ctx->sigmaf_core2_2;
  return g % 2 ? ctx->sigmaf_core1_1 : ctx->sigmaf_core1_2 ;
}

PETSC_STATIC_INLINE PetscReal sigmas(const UserCtx *ctx,PetscInt j,PetscInt i,PetscInt g)
{
  //if ( i>2 && i<17 && j>2 && j<11 ) return g % 2 ? ctx->sigmas_core1_1 : 0.0 ;
  //else if ( i>2 && i<17 && j>10 && j<21 ) return g % 2 ? ctx->sigmas_core2_1 : 0.0 ;
  return g % 2 ? ctx->sigmas_core1_1 : 0.0 ;
}


PetscErrorCode FormInitialGuess(UserCtx *user,Vec U,Mat A1,Mat A2)
{
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       id;
  PetscScalar    *lambda,**phi1,**phi2;
  PetscScalar    *initial_phi1,*initial_phi2;
  Vec            vlambda,vphi1,vphi2;

  PetscFunctionBeginUser;
  ierr = MatSeqAIJGetArray(A1,&initial_phi1);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(A2,&initial_phi2);CHKERRQ(ierr);

  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2);CHKERRQ(ierr);

  ierr = DMDAGetCorners(user->da1,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);//get two direction index
  ierr = DMDAGetInfo(user->da1,0,&N,&M,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);

  if (xs == 0 && ys == 0) lambda[0] = 1.0;
  for ( j=ys ; j < ys+ym ; j++ ){
       for ( i=xs ; i < xs+xm ; i++ ){
          id = j*xm + i ;
          phi1[j][i] = initial_phi1[id] ;
          phi2[j][i] = initial_phi2[id] ;
       }
    }

  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);

  ierr = DMCompositeGather(user->packer,INSERT_VALUES,U,vlambda,vphi1,vphi2);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  UserCtx        *user = (UserCtx*)dummy;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N;
  PetscInt       ys,ym,j,M;/*y direction index*/
  PetscInt       xints,xinte,yints,yinte;
  PetscScalar    *lambda,**phi1,**phi2,*flambda,**fphi1,**fphi2;
  PetscScalar    dx,dy;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,e,b,c,dd,a;
  PetscReal      phi1boundary,phi2boundary;
  Vec            vlambda,vphi1,vphi2,vflambda,vfphi1,vfphi2;

  PetscFunctionBeginUser;
  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->packer,&vflambda,&vfphi1,&vfphi2);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2);CHKERRQ(ierr);

  ierr = DMDAGetCorners(user->da1,&xs,&ys,NULL,&xm,&ym,NULL);CHKERRQ(ierr);//get two direction index
  ierr = DMDAGetInfo(user->da1,0,&N,&M,0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = VecGetArray(vflambda,&flambda);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da1,vfphi1,&fphi1);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da2,vfphi2,&fphi2);CHKERRQ(ierr);
  dx    = 5.0;
  dy    = 5.0;
  phi1boundary = 0.0;
  phi2boundary = 0.0;


  /* residual f_lambda */
  if (xs == 0 && ys == 0 ) { /* only first processor computes this */
    dotphi = 0.0;
    for ( j=ys ; j < ys+ym ; j++ ){
       for ( i=xs ; i < xs+xm ; i++ ){
          dotphi += phi1[j][i]*phi1[j][i] + phi2[j][i]*phi2[j][i] ;
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
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;
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
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;
    }
  }

  /* left boundary */

  if (xs == 0) {
    i     = 0;
    xints = xs + 1;
    /* bottom edge */
    for (j=yints; j<yinte; j++) {
      /*phi1 field*/
      rightside = phi1boundary;
      b         = 0.0;
      c         = 0.0;
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;
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
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi1[j][i] = a * phi1[j][i] - rightside;

      /*phi2 field*/
      rightside = phi2boundary;
      b         = 0.0;
      c         = 0.0;
      dd         = 0.0;
      e         = 0.0;
      a         = 1.0;
      fphi2[j][i] = a * phi2[j][i] - rightside;
    }
  }

  /* inner points */
  for (j=yints; j<yinte; j++){
    for (i=xints; i<xinte; i++){
      // if (phi1[j][i]<0.0)
      // {
      //   phi1[j][i] = 0.0;
      // }
      // if (phi2[j][i]<0.0)
      // {
      //   phi2[j][i] = 0.0;
      // }
      /*phi1 field*/
      rightside = - (0.0121 + 0.0241) * phi1[j][i] + lambda[0] * 0.0085 * phi1[j][i] + lambda[0] * 0.1851 * phi2[j][i] ;
      b = 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      c = 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      dd = 1.0/ ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      e = 1.0 / ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      a = b + c + dd + e ;
      fphi1[j][i] = a* phi1[j][i] - b* phi1[j][i-1] - c* phi1[j][i+1] - dd* phi1[j-1][i] - e* phi1[j+1][i] - rightside ;

      /*phi2 field*/
      rightside = - 0.121 * phi2[j][i] + 0.0241 * phi1[j][i] ;
      b = 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      c = 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      dd = 1.0/ ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      e = 1.0 / ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      a = b + c + dd + e ;
      fphi2[j][i] = a* phi2[j][i] - b* phi2[j][i-1] - c* phi2[j][i+1] - dd* phi2[j-1][i] - e* phi2[j+1][i] - rightside ;

    }
  }


  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = VecRestoreArray(vflambda,&flambda);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vphi1,&phi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da1,vfphi1,&fphi1);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vphi2,&phi2);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da2,vfphi2,&fphi2);CHKERRQ(ierr);

  ierr = DMCompositeGather(user->packer,INSERT_VALUES,FU,vflambda,vfphi1,vfphi2);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  ierr = DMCompositeRestoreLocalVectors(user->packer,&vflambda,&vfphi1,&vfphi2);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *dummy)
{
  UserCtx        *user = (UserCtx*)dummy;
  PetscErrorCode ierr;
  Vec            vlambda,vphi1,vphi2,U;
  PetscReal      *lambda;
  PetscReal      keff;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);
  ierr = DMCompositeGetLocalVectors(user->packer,&vlambda,&vphi1,&vphi2);CHKERRQ(ierr);
  ierr = DMCompositeScatter(user->packer,U,vlambda,vphi1,vphi2);CHKERRQ(ierr);
  ierr = VecGetArray(vlambda,&lambda);CHKERRQ(ierr);

  keff = 1.0 / lambda[0] ;

  ierr = VecRestoreArray(vlambda,&lambda);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, Keff ======= %g\n",its,(double)keff);CHKERRQ(ierr);
  //ierr = VecView(vphi1,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
