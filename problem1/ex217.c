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
} UserCtx;

extern PetscErrorCode FormInitialGuess(UserCtx*,Vec,Mat,Mat);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);
extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode Monitor(SNES,PetscInt,PetscReal,void*);


int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  PetscInt       its;
  Vec            vlambda,vphi1,vphi2;
  Mat            A1,A2;/*Matrix for storage initial guess*/
  SNES           snes;
  UserCtx        user;
  PetscViewer    viewer;
  char           file[PETSC_MAX_PATH_LEN];/* input file name */
  PetscBool      flg;
  PetscInt       nx=20,ny=24,nglobal;
  Vec            U,FU;
  Mat            B;


  ierr = PetscInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

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

  nglobal = 2 * nx * ny + 1 ;
  ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
  ierr = VecSetSizes(U, PETSC_DECIDE, nglobal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(U);CHKERRQ(ierr);
  ierr = VecDuplicate(U,&FU);CHKERRQ(ierr);

  /* Form initial Guess lambda, phi1, phi2 */
  ierr = FormInitialGuess(&user,U,A1,A2);CHKERRQ(ierr);

  /* create nonlinear solver */
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,nglobal,nglobal);CHKERRQ(ierr);
  ierr = MatSetType(B,MATSEQAIJ);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(B,5,NULL);CHKERRQ(ierr);

  ierr = SNESSetFunction(snes,FU,FormFunction,&user);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,NULL,B,FormJacobian,&user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESMonitorSet(snes,Monitor,&user,0);CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,U);CHKERRQ(ierr);
  ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

  /* Output the data */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu_testtwogroups",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
  ierr = VecView(U,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = SNESDestroy(&snes);CHKERRQ(ierr);

  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&FU);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


PetscErrorCode FormInitialGuess(UserCtx *user,Vec U,Mat A1,Mat A2)
{
  PetscErrorCode ierr;
  PetscInt       start,end,row;
  PetscInt       n,ii,jj;/*y direction index*/
  PetscScalar    *u;
  PetscScalar    *initial_phi1,*initial_phi2;

  PetscFunctionBeginUser;
  ierr = MatSeqAIJGetArray(A1,&initial_phi1);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(A2,&initial_phi2);CHKERRQ(ierr);

  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);
  for (row = 0; row < 480; row++) {
    n  = row%(20*24);
    ii = n%20;
    jj = (n-ii)/20;
    u[row]     = initial_phi1[n];
    u[row+480] = initial_phi2[n];
  }
  u[960] = 1.0;
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}



/*
      Evaluates FU = Gradiant(L(w,u,lambda))

*/
PetscErrorCode FormFunction(SNES snes,Vec U,Vec FU,void *dummy)
{
  UserCtx        *user = (UserCtx*)dummy;
  PetscErrorCode ierr;
  PetscInt       row,n,ii,jj;
  PetscInt       start,end,idleft,idright,idup,iddown;
  const PetscScalar *u;
  PetscScalar       *fu;
  PetscScalar    dx,dy;
  PetscScalar    dotphi;//wu-define
  PetscScalar    rightside,e,b,c,dd,a;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecGetArray(FU,&fu);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(U, &start, &end);CHKERRQ(ierr);

  dx    = 5.0;
  dy    = 5.0;

  for (row = 0; row < 480; row++) {
    n  = row%(20*24);
    ii = n%20;
    jj = (n-ii)/20;
    if (ii==0 || jj== 0 || ii==19 || jj==23)
    {
      fu[row]     = u[row] - 0;
      fu[row+480] = u[row+480] - 0;
    }
    else
    {
      idleft  = row - 1 ;
      idright = row + 1 ;
      idup    = row + 20;
      iddown  = row - 20;

      /*phi1 field*/
      rightside = - (0.0121 + 0.0241) * u[row] + u[960] * 0.0085 * u[row] + u[960] * 0.1851 * u[row+480] ;
      b = 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      c = 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      dd = 1.0/ ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      e = 1.0 / ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      a = b + c + dd + e ;
      fu[row] = a* u[row] - b* u[idleft] - c* u[idright] - dd* u[iddown] - e* u[idup] - rightside ;

      /*phi2 field*/
      rightside = - 0.121 * u[row+480] + 0.0241 * u[row] ;
      b = 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      c = 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      dd = 1.0/ ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      e = 1.0 / ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      a = b + c + dd + e ;
      fu[row+480] = a* u[row+480] - b* u[idleft+480] - c* u[idright+480] - dd* u[iddown+480] - e* u[idup+480] - rightside ;
    }
  }

  dotphi = 0.0;
  for (row = 0; row < 960; row++) {
    dotphi += u[row]*u[row];
  }
  fu[960] = dotphi - 1.0 ;

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = VecRestoreArray(FU,&fu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FormJacobian(SNES snes,Vec U,Mat J,Mat B,void *ctx)
{
  UserCtx        *user = (UserCtx*)ctx;
  PetscErrorCode ierr;
  PetscInt       xs,xm,i,N;
  PetscInt       ys,ym,j,M;
  PetscInt       row,row2,n,ii,jj,end;
  Mat            Bll,B11,B22,B12,B21;
  IS             *is;
  const PetscScalar    *u;
  PetscScalar    e,b,c,dd,a,dx,dy;
  PetscScalar    unit;
  PetscScalar    v[5],col[5];
  PetscInt       idleft,idright,idup,iddown;
  PetscInt       idleft2,idright2,idup2,iddown2;
  Vec            vlambda,vphi1,vphi2;

  PetscFunctionBeginUser;
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);
  unit = 1.0;
  end  = 960;
  dx    = 5.0;
  dy    = 5.0;
  for (row = 0; row < 480; row++) {
    n  = row%(20*24);
    ii = n%20;
    jj = (n-ii)/20;
    row2 = row + 480 ;
    if (ii==0 || jj== 0 || ii==19 || jj==23)
    {
      ierr  = MatSetValues(B,1,&row,1,&row,&unit,INSERT_VALUES);CHKERRQ(ierr);
      ierr  = MatSetValues(B,1,&row2,1,&row2,&unit,INSERT_VALUES);CHKERRQ(ierr);
    }
    else
    {
      idleft  = row - 1 ;
      idright = row + 1 ;
      idup    = row + 20;
      iddown  = row - 20;

      /*phi1 field*/
      v[1] = - 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ; col[1] = idleft;
      v[2] = - 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ; col[2] = idright;
      v[3] = - 1.0/  ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ; col[3] = iddown;
      v[4] = - 1.0 / ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ; col[4] = idup;
      v[0] = v[1] + v[2] + v[3] + v[4] + v[0] ;   col[0] = row;
      ierr  = MatSetValues(B,1,&row,5,col,v,INSERT_VALUES);CHKERRQ(ierr);
      // /*phi1 field*/
      // b = - 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      // ierr  = MatSetValues(B,1,&row,1,&idleft,&b,INSERT_VALUES);CHKERRQ(ierr);
      // c = - 1.0 / ( dx * dx / (2* 1.267) + dx * dx / (2* 1.267 ) ) ;
      // ierr  = MatSetValues(B,1,&row,1,&idright,&c,INSERT_VALUES);CHKERRQ(ierr);
      // dd = - 1.0/  ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      // ierr  = MatSetValues(B,1,&row,1,&iddown,&dd,INSERT_VALUES);CHKERRQ(ierr);
      // e = - 1.0 / ( dy * dy / (2* 1.267) + dy * dy / (2* 1.267 ) ) ;
      // ierr  = MatSetValues(B,1,&row,1,&idup,&e,INSERT_VALUES);CHKERRQ(ierr);
      // a = 0.0 - b - c - dd - e ;
      // ierr  = MatSetValues(B,1,&row,1,&row,&a,INSERT_VALUES);CHKERRQ(ierr);

      idleft2  = row - 1  + 480 ;
      idright2 = row + 1  + 480 ;
      idup2    = row + 20 + 480 ;
      iddown2  = row - 20 + 480 ;

      /*phi2 field*/
      b = - 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      ierr  = MatSetValues(B,1,&row2,1,&idleft2,&b,INSERT_VALUES);CHKERRQ(ierr);
      c = - 1.0 / ( dx * dx / (2* 0.354) + dx * dx / (2* 0.354 ) ) ;
      ierr  = MatSetValues(B,1,&row2,1,&idright2,&c,INSERT_VALUES);CHKERRQ(ierr);
      dd = - 1.0/  ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      ierr  = MatSetValues(B,1,&row2,1,&iddown2,&dd,INSERT_VALUES);CHKERRQ(ierr);
      e = - 1.0 / ( dy * dy / (2* 0.354) + dy * dy / (2* 0.354 ) ) ;
      ierr  = MatSetValues(B,1,&row2,1,&idup2,&e,INSERT_VALUES);CHKERRQ(ierr);
      a = 0.0 - b - c - dd - e ;
      ierr  = MatSetValues(B,1,&row2,1,&row2,&a,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr  = MatSetValues(B,1,&end,1,&end,&unit,INSERT_VALUES);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd  (B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  if (J != B) {
    ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);

}

PetscErrorCode Monitor(SNES snes,PetscInt its,PetscReal fnorm,void *dummy)
{
  UserCtx        *user = (UserCtx*)dummy;
  PetscErrorCode ierr;
  Vec            U;
  const PetscScalar *u;
  PetscReal      keff;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes,&U);CHKERRQ(ierr);
  ierr = VecGetArrayRead(U,&u);CHKERRQ(ierr);

  keff = 1.0 / u[960] ;

  ierr = VecRestoreArrayRead(U,&u);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"iter = %D, Keff ======= %g\n",its,(double)keff);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

