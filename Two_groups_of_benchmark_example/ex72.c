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

/*
   User-defined routines
*/
PetscErrorCode MatXFModel(PetscInt n,PetscInt m,Mat XF);
PetscErrorCode MatLModel(PetscInt n,PetscInt m,Mat L);


int main(int argc,char **argv)
{
  Mat            XF,L;             /* matrices */
  EPS            eps;             /* eigenproblem solver context */
  ST             st;
  KSP            ksp;
  EPSType        type;
  PetscReal      tol;
  Vec            xr,xi;
  PetscInt       nev,maxit,its,lits;
  PetscInt       n=20,m=24,N;
  PetscBool      terse;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Load the matrices that define the eigensystem, XF x = \lambda L x 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);

  N = n*m*2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThis is a neutron eigenvalue problem.\n\n");CHKERRQ(ierr);

  /* Create the fission matrix */
  ierr = PetscPrintf(PETSC_COMM_WORLD," Creating the XF Matrix...\n");CHKERRQ(ierr);
  
  ierr = MatCreate(PETSC_COMM_WORLD,&XF);CHKERRQ(ierr);
  ierr = MatSetSizes(XF,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(XF);CHKERRQ(ierr);
  ierr = MatSetUp(XF);CHKERRQ(ierr);
  ierr = MatXFModel(n,m,XF);CHKERRQ(ierr);


  /* Create the scattering matrix */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nCreating the L matrix.\n");CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&L);CHKERRQ(ierr);
  ierr = MatSetSizes(L,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L);CHKERRQ(ierr);
  ierr = MatSetUp(L);CHKERRQ(ierr);
  ierr = MatLModel(n,m,L);CHKERRQ(ierr);


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
  // if (nini > 0) {
  //   ierr = VecDestroyVecs(nini,&Iv);CHKERRQ(ierr);
  // }
  // if (ncon > 0) {
  //   ierr = VecDestroyVecs(ncon,&Cv);CHKERRQ(ierr);
  // }
  ierr = SlepcFinalize();
  return ierr;
}




PetscErrorCode MatXFModel(PetscInt n,PetscInt m,Mat XF)
{
  PetscReal       sigmaf1,sigmaf2;
  PetscInt        i,j,II,Iend,Interval;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  Iend = n*m ;
  Interval = n*m ;
  for (II = 0; II < Iend ; II++ ) {
    i = II/n; j = II-i*n;
    sigmaf1 = 0.0085;
    sigmaf2 = 0.1851; 
    ierr = MatSetValue(XF,II,II,sigmaf1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValue(XF,II,II+Interval,sigmaf2,INSERT_VALUES);CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(XF,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(XF,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MatLModel(PetscInt n, PetscInt m, Mat L)
{
  PetscReal       a,d,b,c,e,a1,b1,c1,d1,e1;
  PetscReal       dx,dy;
  PetscReal       sigmas;
  PetscInt        i,j,II,Iend,Interval;
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;
  Iend = n*m ;
  Interval = n*m ;
  dx = 5.0;
  dy = 5.0;
  for (II = 0; II<Iend ; II++) {
    i = II/n; j = II-i*n;
    if (i>0) {
      a =  - 1.0 / ( dy*dy / (2.0 * 1.267) + dy*dy / ( 2.0 * 1.267 )) ;
      ierr = MatSetValue(L,II,II-n,a,INSERT_VALUES);CHKERRQ(ierr); 
    }
    if (i<m-1) { 
      d = - 1.0 / ( dy*dy / (2.0 * 1.267) + dy*dy / ( 2.0 * 1.267 )) ;
      ierr = MatSetValue(L,II,II+n,d,INSERT_VALUES);CHKERRQ(ierr); 
    }
    if (j>0) { 
      b = - 1.0 / ( dx*dx / (2.0 * 1.267) + dx*dx / ( 2.0 * 1.267 )) ;
      ierr = MatSetValue(L,II,II-1,b,INSERT_VALUES);CHKERRQ(ierr); 
    }
    if (j<n-1) { 
      c = - 1.0 / ( dx*dx / (2.0 * 1.267) + dx*dx / ( 2.0 * 1.267 )) ;
      ierr = MatSetValue(L,II,II+1,c,INSERT_VALUES);CHKERRQ(ierr); 
    }

    e = 0.0121 + 0.0241 + 4.0 / ( dx*dx / (2.0 * 1.267) + dx*dx / ( 2.0 * 1.267 ));
    ierr = MatSetValue(L,II,II,e,INSERT_VALUES);CHKERRQ(ierr);
  }
   for (II = 0 ; II<Iend; II++) {
     i = II/n; j = II-i*n;
     if (i>0) {
       a1 =  - 1.0 / ( dy*dy / (2.0 * 0.354) + dy*dy / ( 2.0 * 0.354)) ;
       ierr = MatSetValue(L,II+Interval,II-n+Interval,a1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     if (i<m-1) { 
       d1 = - 1.0 / ( dy*dy / (2.0 * 0.354) + dy*dy / ( 2.0 * 0.354)) ;
       ierr = MatSetValue(L,II+Interval,II+n+Interval,d1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     if (j>0) { 
       b1 = - 1.0 / ( dx*dx / (2.0 * 0.354) + dx*dx / ( 2.0 * 0.354)) ;
       ierr = MatSetValue(L,II+Interval,II-1+Interval,b1,INSERT_VALUES);CHKERRQ(ierr); 
     }
     if (j<n-1) { 
       c1 = - 1.0 / ( dx*dx / (2.0 * 0.354) + dx*dx / ( 2.0 * 0.354)) ;
       ierr = MatSetValue(L,II+Interval,II+1+Interval,c1,INSERT_VALUES);CHKERRQ(ierr); 
     }

     e1 = 0.121 + 4.0 / ( dx*dx / (2.0 * 0.354) + dx*dx / ( 2.0 * 0.354)) ;
     ierr = MatSetValue(L,II+Interval,II+Interval,e1,INSERT_VALUES);CHKERRQ(ierr);
   }
   for (II = 0; II < Iend ; II++ ) {
     i = II/n; j = II-i*n;
     sigmas = -0.0241; 
     ierr = MatSetValue(L,II+Interval,II,sigmas,INSERT_VALUES);CHKERRQ(ierr);
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
