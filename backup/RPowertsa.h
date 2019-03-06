#ifndef RPOWERTSA_H
#define RPOWERTSA_H
#include <petsc.h>
class RPowertsa
{
private:
	int n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11;
	PetscScalar *CVal,*CVbl,*CVdxl;
    PetscScalar *CVdyl,*CVsigdl,*CVsigs0l;
    PetscScalar *CVsigtl,*CVvsiga,*CVv,*CVvsigfl;
    
public:
	PetscScalar *Qf;
	PetscScalar *QQQ;
	PetscScalar *Cdflux,*Cdfluxold,*fluxold;
	PetscScalar beta[6],lamda[6],dt,keff;
	PetscScalar v[4];
	PetscReal   oshift;
	RPowertsa();
	~RPowertsa();
	void setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9);
	double getnum();
	PetscScalar gettsa1(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar gettsa2(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar gettsa3(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar gettsa4(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar gettsa0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4);
	PetscScalar gettsrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff);
	PetscScalar getsigf(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getcvv(PetscInt i, PetscInt j);
	PetscScalar getdflux(PetscInt i, PetscInt j, PetscInt gdl, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4);

	void iniQf();
	void addQQQ();
	void perturbation(PetscReal time);
};


#endif