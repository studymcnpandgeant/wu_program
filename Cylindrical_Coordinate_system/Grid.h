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
	PetscScalar getchisigmaf(PetscInt II, PetscInt chi,PetscInt group);
	PetscScalar getsigmat(PetscInt II,PetscInt group);
	PetscScalar getsigmas(PetscInt II,PetscInt group11,PetscInt group22);
	PetscScalar getrzae(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getrzaw(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getrzan(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getrzas(PetscInt i,PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getrzchisigmaf(PetscInt j,PetscInt II, PetscInt chi,PetscInt group);
	PetscScalar getrzsigmat(PetscInt j,PetscInt II,PetscInt group);
	PetscScalar getrzsigmas(PetscInt j,PetscInt II,PetscInt group11,PetscInt group22);
	
};






#endif