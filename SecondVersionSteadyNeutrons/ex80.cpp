/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2018, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] = "Solves a generalized eigensystem Ax=kBx with matrices loaded from a file.\n"
  "The command line options are:\n"
  "  -f1 <filename> -f2 <filename>, PETSc binary files containing A and B.\n"
  "  -evecs <filename>, output file to save computed eigenvectors.\n"
  "  -ninitial <nini>, number of user-provided initial guesses.\n"
  "  -finitial <filename>, binary file containing <nini> vectors.\n"
  "  -nconstr <ncon>, number of user-provided constraints.\n"
  "  -fconstr <filename>, binary file containing <ncon> vectors.\n\n";

#include <slepceps.h>
#include "Grid.h"

/*
   User-defined routines
*/
typedef struct
{
  PetscScalar *corearray,*materialarray;
  PetscScalar *rcor,*zcor;
  PetscScalar *rrcor,*zzcor;
} AppCtx;

PetscErrorCode MatXFModel(PetscInt n,PetscInt m,Mat XF,Grid *grid);
PetscErrorCode MatLModel(PetscInt n,PetscInt m,Mat L,Grid *grid);


int main(int argc,char **argv)
{
  AppCtx         user;
  Mat            XF,L;             /* matrices */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;
  KSP            ksp;
  EPSType        type;
  PetscReal      tol;
  Vec            xr,xi;
  Vec            u1,u2,u3,u4,u5,u6;    /* These vectors used for reading binary file */
  PetscViewer    viewer;
  PetscInt       nev,maxit,its,lits;
  PetscInt       n=16,m=35,N;   /* x-16grids   y-35grids */
  PetscBool      terse;
  PetscErrorCode ierr;
  Grid           grid;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* Reading vectors in binary format */
  /* Read four vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
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
  ierr =VecGetArray(u1,&user.materialarray);
  ierr =VecGetArray(u2,&user.corearray);
  ierr =VecGetArray(u3,&user.rcor);
  ierr =VecGetArray(u4,&user.zcor);
  ierr =VecGetArray(u5,&user.rrcor);
  ierr =VecGetArray(u6,&user.zzcor);

  grid.setinitial(user.materialarray,user.corearray,user.rcor,user.zcor,user.rrcor,user.zzcor);



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, XF x = \lambda L x 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);

  N = n*m*4;/* There are four groups */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThis is a neutron eigenvalue problem.\n\n");CHKERRQ(ierr);

  /* Create the fission matrix */
  ierr = PetscPrintf(PETSC_COMM_WORLD," Creating the XF Matrix...\n");CHKERRQ(ierr);
  
  ierr = MatCreate(PETSC_COMM_WORLD,&XF);CHKERRQ(ierr);
  ierr = MatSetSizes(XF,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(XF);CHKERRQ(ierr);
  ierr = MatSetUp(XF);CHKERRQ(ierr);
  ierr = MatXFModel(n,m,XF,&grid);CHKERRQ(ierr);


  /* Create the scattering matrix */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nCreating the L matrix.\n");CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&L);CHKERRQ(ierr);
  ierr = MatSetSizes(L,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L);CHKERRQ(ierr);
  ierr = MatSetUp(L);CHKERRQ(ierr);
  ierr = MatLModel(n,m,L,&grid);CHKERRQ(ierr);


  ierr = MatCreateVecs(XF,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(XF,NULL,&xi);CHKERRQ(ierr);

  /*
     Read initial guesses if available
  */
  // ierr = PetscOptionsGetInt(NULL,NULL,"-ninitial",&nini,&flg);CHKERRQ(ierr);
  // if (flg) {
  //   if (nini<=0) SETERRQ(PETSC_COMM_WORLD,1,"The number of initial vectors must be >0");
  //   ierr = PetscOptionsGetString(NULL,NULL,"-finitial",filename,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  //   if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must specify the name of the file containing the initial vectors");
  //   ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  //   ierr = VecDuplicateVecs(xr,nini,&Iv);CHKERRQ(ierr);
  //   for (i=0;i<nini;i++) {
  //     ierr = VecLoad(Iv[i],viewer);CHKERRQ(ierr);
  //   }
  //   ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  // }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a generalized eigenvalue problem
  */
  ierr = EPSSetOperators(eps,XF,L);CHKERRQ(ierr);

  /*
     If the user provided initial guesses or constraints, pass them here
  */
  // ierr = EPSSetInitialSpace(eps,nini,Iv);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STGetKSP(st,&ksp);CHKERRQ(ierr);
  ierr = KSPGetTotalIterations(ksp,&lits);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %D\n",lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Show detailed info unless -terse option is given by user
   */
  ierr = PetscOptionsHasName(NULL,NULL,"-terse",&terse);CHKERRQ(ierr);
  if (terse) {
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,NULL);CHKERRQ(ierr);
  } else {
    ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL);CHKERRQ(ierr);
    ierr = EPSReasonView(eps,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  /*
     Save eigenvectors, if requested
  */
//   ierr = PetscOptionsGetString(NULL,NULL,"-evecs",filename,PETSC_MAX_PATH_LEN,&evecs);CHKERRQ(ierr);
//   ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
//   if (nconv>0 && evecs) {
//     ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
//     ierr = EPSIsHermitian(eps,&ishermitian);CHKERRQ(ierr);
//     for (i=0;i<nconv;i++) {
//       ierr = EPSGetEigenvector(eps,i,xr,xi);CHKERRQ(ierr);
//       ierr = VecView(xr,viewer);CHKERRQ(ierr);
// #if !defined(PETSC_USE_COMPLEX)
//       if (!ishermitian) { ierr = VecView(xi,viewer);CHKERRQ(ierr); }
// #endif
//     }
//     ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
//   }

  /*
     Free work space
  */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&XF);CHKERRQ(ierr);
  ierr = MatDestroy(&L);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  /* Free reading structures */
  ierr = VecRestoreArray(u1,&user.materialarray);
  ierr = VecRestoreArray(u2,&user.corearray);
  ierr = VecRestoreArray(u3,&user.rcor);
  ierr = VecRestoreArray(u4,&user.zcor);
  ierr = VecDestroy(&u1);CHKERRQ(ierr);
  ierr = VecDestroy(&u2);CHKERRQ(ierr);
  ierr = VecDestroy(&u3);CHKERRQ(ierr);
  ierr = VecDestroy(&u4);CHKERRQ(ierr);
  // if (nini > 0) {
  //   ierr = VecDestroyVecs(nini,&Iv);CHKERRQ(ierr);
  // }
  // if (ncon > 0) {
  //   ierr = VecDestroyVecs(ncon,&Cv);CHKERRQ(ierr);
  // }
  ierr = SlepcFinalize();
  return ierr;
}




PetscErrorCode MatXFModel(PetscInt n,PetscInt m,Mat XF,Grid *grid)
{
  PetscInt        II,Iend,Interval1,Interval2,Interval3;
  PetscErrorCode  ierr;
  PetscReal       chisigmaf1,chisigmaf2,chisigmaf3,chisigmaf4;
  PetscInt        chigroup1=1,chigroup2=2,chigroup3=3;
  PetscInt        group1=1,group2=2,group3=3,group4=4;

  PetscFunctionBeginUser;
  Iend = n*m ;
  Interval1 = n*m ;
  Interval2 = n*m*2 ;
  Interval3 = n*m*3 ;
  /* define the chi */
  for (II = 0; II < Iend ; II++ ) {
    chisigmaf1 = grid->getchisigmaf(II,chigroup1,group1);
    chisigmaf2 = grid->getchisigmaf(II,chigroup1,group2);
    chisigmaf3 = grid->getchisigmaf(II,chigroup1,group3);
    chisigmaf4 = grid->getchisigmaf(II,chigroup1,group4);
    ierr = MatSetValue(XF,II,II,chisigmaf1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II,II+Interval1,chisigmaf2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II,II+Interval2,chisigmaf3,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II,II+Interval3,chisigmaf4,INSERT_VALUES);CHKERRQ(ierr);
  }

  for (II = 0; II < Iend ; II++ ) {
    chisigmaf1 = grid->getchisigmaf(II,chigroup2,group1);
    chisigmaf2 = grid->getchisigmaf(II,chigroup2,group2);
    chisigmaf3 = grid->getchisigmaf(II,chigroup2,group3);
    chisigmaf4 = grid->getchisigmaf(II,chigroup2,group4);
    ierr = MatSetValue(XF,II+Interval1,II,chisigmaf1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval1,II+Interval1,chisigmaf2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval1,II+Interval2,chisigmaf3,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval1,II+Interval3,chisigmaf4,INSERT_VALUES);CHKERRQ(ierr);
  }

  for (II = 0; II < Iend ; II++ ) {
    chisigmaf1 = grid->getchisigmaf(II,chigroup3,group1);
    chisigmaf2 = grid->getchisigmaf(II,chigroup3,group2);
    chisigmaf3 = grid->getchisigmaf(II,chigroup3,group3);
    chisigmaf4 = grid->getchisigmaf(II,chigroup3,group4);
    ierr = MatSetValue(XF,II+Interval2,II,chisigmaf1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval2,II+Interval1,chisigmaf2,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval2,II+Interval2,chisigmaf3,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II+Interval2,II+Interval3,chisigmaf4,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(XF,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(XF,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatLModel(PetscInt n, PetscInt m, Mat L,Grid *grid)
{
  PetscReal       a,d,b,c,e,a1,b1,c1,d1,e1;
  PetscReal       a2,d2,b2,c2,e2,a3,b3,c3,d3,e3;
  PetscInt        i,j,II,Iend,Interval1,Interval2,Interval3;
  PetscReal       sigmat1,sigmas11,sigmat2,sigmas22,sigmat3,sigmas33,sigmat4,sigmas44;
  PetscReal       sigmas12,sigmas13,sigmas14,sigmas23,sigmas24,sigmas34;
  PetscInt        group1=1,group2=2,group3=3,group4=4;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  Iend = n*m ;
  Interval1 = n*m ;
  Interval2 = 2*n*m ;
  Interval3 = 3*n*m ;
  for (II = 0; II<Iend ; II++) {
    i = II/n; j = II-i*n;
    if (i>0) {
      a = grid->getas(i,j,II,group1);
      ierr = MatSetValue(L,II,II-n,a,INSERT_VALUES);CHKERRQ(ierr); 
    }
    else
    {
      a = grid->getas(i,j,II,group1);
    }

    if (i<m-1) {
      d = grid->getan(i,j,II,group1);
      ierr = MatSetValue(L,II,II+n,d,INSERT_VALUES);CHKERRQ(ierr); 
    }
    else
    {
      d = grid->getan(i,j,II,group1);
    }

    if (j>0) { 
      b = grid->getaw(i,j,II,group1);
      ierr = MatSetValue(L,II,II-1,b,INSERT_VALUES);CHKERRQ(ierr); 
    }
    else
    {
      b = grid->getaw(i,j,II,group1);
    }

    if (j<n-1) { 
      c = grid->getae(i,j,II,group1);
      ierr = MatSetValue(L,II,II+1,c,INSERT_VALUES);CHKERRQ(ierr); 
    }
    else
    {
      c = grid->getae(i,j,II,group1);
    }

    sigmat1 = grid->getsigmat(II,group1);
    sigmas11 = grid->getsigmas(II,group1,group1);
    e = (sigmat1 - sigmas11) - a - d - b - c;
    ierr = MatSetValue(L,II,II,e,INSERT_VALUES);CHKERRQ(ierr);
  }
   for (II = 0 ; II<Iend; II++) {
     i = II/n; j = II-i*n;
     if (i>0) {
       a1 = grid->getas(i,j,II,group2);
       ierr = MatSetValue(L,II+Interval1,II-n+Interval1,a1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
       a1 = grid->getas(i,j,II,group2);
     }

     if (i<m-1) { 
       d1 = grid->getan(i,j,II,group2);
       ierr = MatSetValue(L,II+Interval1,II+n+Interval1,d1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
       d1 = grid->getan(i,j,II,group2);
     }

     if (j>0) { 
       b1 = grid->getaw(i,j,II,group2);
       ierr = MatSetValue(L,II+Interval1,II-1+Interval1,b1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
       b1 = grid->getaw(i,j,II,group2);
     }

     if (j<n-1) { 
       c1 = grid->getae(i,j,II,group2);
       ierr = MatSetValue(L,II+Interval1,II+1+Interval1,c1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
       c1 = grid->getae(i,j,II,group2);
     }

     sigmat2 = grid->getsigmat(II,group2);
     sigmas22 = grid->getsigmas(II,group2,group2);
     e1 = (sigmat2 - sigmas22) - a1 - d1 - b1 - c1;
     ierr = MatSetValue(L,II+Interval1,II+Interval1,e1,INSERT_VALUES);CHKERRQ(ierr);
   }
  for (II = 0 ; II<Iend; II++) {
     i = II/n; j = II-i*n;
     if (i>0) {
       a2 = grid->getas(i,j,II,group3);
       ierr = MatSetValue(L,II+Interval2,II-n+Interval2,a2,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	a2 = grid->getas(i,j,II,group3);
     }

     if (i<m-1) { 
       d2 = grid->getan(i,j,II,group3);
       ierr = MatSetValue(L,II+Interval2,II+n+Interval2,d2,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	d2 = grid->getan(i,j,II,group3);
     }

     if (j>0) { 
       b2 = grid->getaw(i,j,II,group3);
       ierr = MatSetValue(L,II+Interval2,II-1+Interval2,b2,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	b2 = grid->getaw(i,j,II,group3);
     }

     if (j<n-1) { 
       c2 = grid->getae(i,j,II,group3);
       ierr = MatSetValue(L,II+Interval2,II+1+Interval2,c2,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	c2 = grid->getae(i,j,II,group3);
     }

     sigmat3 = grid->getsigmat(II,group3);
     sigmas33 = grid->getsigmas(II,group3,group3);
     e2 = (sigmat3 - sigmas33) - a2 -d2 - b2 - c2 ;
     ierr = MatSetValue(L,II+Interval2,II+Interval2,e2,INSERT_VALUES);CHKERRQ(ierr);
   }
  for (II = 0 ; II<Iend; II++) {
     i = II/n; j = II-i*n;
     if (i>0) {
       a3 = grid->getas(i,j,II,group4);
       ierr = MatSetValue(L,II+Interval3,II-n+Interval3,a3,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	a3 = grid->getas(i,j,II,group4);
     }

     if (i<m-1) {
       d3 = grid->getan(i,j,II,group4);
       ierr = MatSetValue(L,II+Interval3,II+n+Interval3,d3,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	d3 = grid->getan(i,j,II,group4);
     }

     if (j>0) {
       b3 = grid->getaw(i,j,II,group4);
       ierr = MatSetValue(L,II+Interval3,II-1+Interval3,b3,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	b3 = grid->getaw(i,j,II,group4);
     }

     if (j<n-1) {  
       c3 = grid->getae(i,j,II,group4);
       ierr = MatSetValue(L,II+Interval3,II+1+Interval3,c3,INSERT_VALUES);CHKERRQ(ierr); 
     }
     else
     {
     	c3 = grid->getae(i,j,II,group4);
     }

     sigmat4 = grid->getsigmat(II,group4);
     sigmas44 = grid->getsigmas(II,group4,group4);
     e3 = (sigmat4 - sigmas44) - a3 - d3 - b3 - c3 ;
     ierr = MatSetValue(L,II+Interval3,II+Interval3,e3,INSERT_VALUES);CHKERRQ(ierr);
   }
   for (II = 0; II < Iend ; II++ ) {
     i = II/n; j = II-i*n;
     sigmas12 = - grid->getsigmas(II,group1,group2); 
     ierr = MatSetValue(L,II+Interval1,II,sigmas12,INSERT_VALUES);CHKERRQ(ierr);
   }
  for (II = 0; II < Iend ; II++ ) {
     i = II/n; j = II-i*n;
     sigmas13 = - grid->getsigmas(II,group1,group3);
     sigmas23 = - grid->getsigmas(II,group2,group3); 
     ierr = MatSetValue(L,II+Interval2,II,sigmas13,INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(L,II+Interval2,II+Interval1,sigmas23,INSERT_VALUES);CHKERRQ(ierr);
   }
  for (II = 0; II < Iend ; II++ ) {
     i = II/n; j = II-i*n;
     sigmas14 = - grid->getsigmas(II,group1,group4);
     sigmas24 = - grid->getsigmas(II,group2,group4);
     sigmas34 = - grid->getsigmas(II,group3,group4); 
     ierr = MatSetValue(L,II+Interval3,II,sigmas14,INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(L,II+Interval3,II+Interval1,sigmas24,INSERT_VALUES);CHKERRQ(ierr);
     ierr = MatSetValue(L,II+Interval3,II+Interval2,sigmas34,INSERT_VALUES);CHKERRQ(ierr);
   }
  ierr = MatAssemblyBegin(L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*TEST

   test:
      suffix: 1
      args: -f1 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62a.petsc -f2 ${SLEPC_DIR}/share/slepc/datafiles/matrices/bfw62b.petsc -eps_nev 4 -terse
      requires: double !complex !define(PETSC_USE_64BIT_INDICES)

   test:
      suffix: ciss_1
      args: -f1 ${DATAFILESPATH}/matrices/complex/mhd1280a.petsc -f2 ${DATAFILESPATH}/matrices/complex/mhd1280b.petsc -eps_type ciss -eps_ciss_usest 0 -eps_ciss_quadrule chebyshev -rg_type ring -rg_ring_center 0 -rg_ring_radius .5 -rg_ring_width 0.2 -rg_ring_startangle .25 -rg_ring_endangle .5 -terse
      requires: complex datafilespath
      timeoutfactor: 2

TEST*/

