
static char help[] = "Tests binary I/O of vectors and illustrates the use of user-defined event logging.\n\n";

#include <petscvec.h>

/* Note:  Most applications would not read and write a vector within
  the same program.  This example is intended only to demonstrate
  both input and output. */

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size;
  PetscInt       m = 10,a3;
  PetscScalar    *v1;
  PetscScalar    *v2;
  Vec            u1,u2;
  PetscViewer    viewer;
  PetscReal      a1;
  PetscReal      a2;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);

  /* Read in vector in binary format */

  /* Read two vectors in binary format */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from vector.dat ...\n \n");CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"wu-vec",FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u1);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,&u2);CHKERRQ(ierr);
  ierr = VecLoad(u1,viewer);CHKERRQ(ierr);
  ierr = VecLoad(u2,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* Initialize the two array */
  ierr =VecGetArray(u1,&v1);
  ierr =VecGetArray(u2,&v2);
  a1 = v1[1];
  a2 = v2[17];
  a3 = a2;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A1 element is %g \n", a1 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A2 element is %d \n", a3 );CHKERRQ(ierr);


  /* Free data structures */
  ierr = VecRestoreArray(u1,&v1);
  ierr = VecRestoreArray(u2,&v2);
  ierr = VecDestroy(&u1);CHKERRQ(ierr);
  ierr = VecDestroy(&u2);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

     test:
       nsize: 1
       requires: mpiio

     test:
       suffix: 2
       nsize: 2
       requires: mpiio

TEST*/
