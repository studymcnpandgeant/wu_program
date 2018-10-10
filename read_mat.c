
static char help[] = "Tests MatCreateComposite()\n\n";

/*T
   Concepts: Mat^composite matrices
   Processors: n
T*/

/*
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h    - vectors
     petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers
*/
#include <petscmat.h>

int main(int argc,char **args)
{
  Mat            A,B;                       /* matrix */
  PetscViewer    fd;                      /* viewer */
  char           file[PETSC_MAX_PATH_LEN];            /* input file name */
  PetscErrorCode ierr;
  PetscBool      flg;
  PetscScalar    *a,*b;
  PetscReal      a1,a2;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  /*
     Determine files from which we read the two linear systems
     (matrix and right-hand-side vector).
  */
  ierr = PetscOptionsGetString(NULL,NULL,"-f",file,PETSC_MAX_PATH_LEN,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary file with the -f option");

  /*
     Open binary file.  Note that we use FILE_MODE_READ to indicate
     reading from this file.
  */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file,FILE_MODE_READ,&fd);CHKERRQ(ierr);

  /*
     Load the matrix; then destroy the viewer.
  */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetType(A,MATSEQAIJ);
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetType(B,MATSEQAIJ);
  ierr = MatLoad(A,fd);CHKERRQ(ierr);
  ierr = MatLoad(B,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  /* Transform the data to a array */
  ierr = MatSeqAIJGetArray(A,&a);CHKERRQ(ierr);
  ierr = MatSeqAIJGetArray(B,&b);CHKERRQ(ierr);
  a1   = a[22];
  a2   = b[22];



  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A1 element is %g \n", a1 );CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"The A2 element is %g \n", a2 );CHKERRQ(ierr);


  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */

  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}





/*TEST

   test:
      args: -f ${PETSCPATH}/src/mat/examples/tutorials/wu-read2

TEST*/
