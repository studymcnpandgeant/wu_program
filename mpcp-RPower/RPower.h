#ifndef RPOWER_H
#define RPOWER_H
#include <petsc.h>
class RPower
{
private:
	int n1,n2,n3,n4,n5,n6,n7,n8,n9;
	PetscScalar *CVal,*CVbl,*CVdxl;
    PetscScalar *CVdyl,*CVsigdl,*CVsigs0l;
    PetscScalar *CVsigtl,*CVv,*CVvsigfl;
public:
	PetscScalar *Qf;
	PetscScalar *QQQ;
	RPower();
	~RPower();
	void setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9);
	double getnum();
	PetscScalar geta1(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar geta2(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar geta3(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar geta4(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar geta0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4);
	PetscScalar getrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff);
	PetscScalar getsigf(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getcvv(PetscInt i, PetscInt j);
	void iniQf();
	void addQQQ();
};


#endif