#ifndef GRID_H
#define GRID_H
#include <petsc.h>

class Grid
{
private:
	int n1,n2,n3,n4,n5,n6;
	PetscScalar *materialarray,*corearray;
    PetscScalar *dr,*dz;
    PetscScalar *rcor,*zcor;
public:
	Grid();
	~Grid();
	void setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6);
	PetscScalar getmaterialarray(PetscInt a); 
	PetscScalar getae(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getaw(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getan(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getas(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	
};






#endif