
static char help[] = "Basic vector routines.\n\n";

/*T
   Concepts: vectors^basic routines;
   Processors: n
T*/



/*
  Include "petscvec.h" so that we can use vectors.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscis.h     - index sets
     petscviewer.h - viewers
*/

#include <petscvec.h>
#include "Point.h"

typedef struct
{
  PetscScalar *corearray,*materialarray;
  PetscScalar *rcor,*zcor;
} AppCtx;

PetscErrorCode foo(void*);

int main(int argc,char **args)
{
  AppCtx         user;
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  PetscInt       m = 10,a3;
  // PetscScalar    *v1;
  // PetscScalar    *v2;
  Vec            u1,u2,u3,u4;
  PetscViewer    viewer;
  PetscReal      a1;
  PetscReal      a2,a4,a5,a6;
  Point          p1;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);

  /* Read in vector in binary format */

  /* Read two vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-vec4",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u2);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u3);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u4);CHKERRQ(ierr);
  ierr = VecLoad(u1,viewer);CHKERRQ(ierr);
  ierr = VecLoad(u2,viewer);CHKERRQ(ierr);
  ierr = VecLoad(u3,viewer);CHKERRQ(ierr);
  ierr = VecLoad(u4,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the two array */
  ierr =VecGetArray(u1,&user.materialarray);
  ierr =VecGetArray(u2,&user.corearray);
  ierr =VecGetArray(u3,&user.rcor);
  ierr =VecGetArray(u4,&user.zcor);
  a1 = user.materialarray[1];
  a2 = user.corearray[17];
  a3 = a2;
  a4 = user.rcor[5];
  a5 = user.zcor[5];
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A1 element is %g \n", a1 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A2 element is %d \n", a3 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The r 6 element is %g \n", a4 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The z 6 element is %g \n", a5 );CHKERRQ(ierr);
  //ierr = foo(&user);CHKERRQ(ierr);
  p1.setdx(a1);
  a6 = p1.getdx();
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The Point p1 is %g \n", a6 );CHKERRQ(ierr);

  /* Free data structures */
  ierr = VecRestoreArray(u1,&user.materialarray);
  ierr = VecRestoreArray(u2,&user.corearray);
  ierr = VecRestoreArray(u3,&user.rcor);
  ierr = VecRestoreArray(u4,&user.zcor);
  ierr = VecDestroy(&u1);CHKERRQ(ierr);
  ierr = VecDestroy(&u2);CHKERRQ(ierr);
  ierr = VecDestroy(&u3);CHKERRQ(ierr);
  ierr = VecDestroy(&u4);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}


/*TEST

   test:

   test:
      suffix: 2
      nsize: 2
      output_file: output/ex1_1.out

   test:
      suffix: 2_cuda
      nsize: 2
      args: -vec_type cuda
      output_file: output/ex1_1.out
      requires: veccuda

   test:
      suffix: cuda
      args: -vec_type cuda
      output_file: output/ex1_1.out
      requires: veccuda

TEST*/
