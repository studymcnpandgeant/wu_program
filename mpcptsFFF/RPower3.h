#ifndef RPOWER3_H
#define RPOWER3_H
#include <petsc.h>
#include <cmath>

class RPower3
{
private:
	int n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12;
	PetscScalar *CVal,*CVbl,*CVdxl;
    PetscScalar *CVdyl,*CVsigdl,*CVsigs0l;
    PetscScalar *CVsigtl,*CVvsiga,*CVv,*CVvsigfl;
    PetscScalar *CVsigtlfixed,*CVsigdlfixed,*CVvsigflfixed,*CVsigs0lfixed;
    PetscScalar *Bco;
public:
	PetscScalar *Qf;
	PetscScalar *QQQ;
	PetscScalar *Cdflux,*Cdfluxold,*fluxold;
	PetscScalar beta[6],lamda[6],dt,keff;
	PetscScalar v[4];
	PetscReal   oshift;
	RPower3();
	~RPower3();
	void setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10);
	double getnum();
	void computeQfth(const PetscScalar *u);
	void updateCrossSection(PetscScalar *Ts);
	void updateVoidCrossSection(PetscScalar *Ts);
	//void updateArtificialCrossSection(PetscScalar *Ts);
	PetscScalar getRPa1(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getRPa2(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getRPa3(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getRPa4(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getRPa0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4);
	PetscScalar getRPrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff, PetscReal time);
	PetscScalar getsigf(PetscInt i, PetscInt j, PetscInt group);
	PetscScalar getcvv(PetscInt i, PetscInt j);
	PetscScalar getdflux(PetscInt i, PetscInt j, PetscInt gdl, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4);
	void iniQf();
	void addQQQ();
	void perturbation(PetscReal time);
};


#endif