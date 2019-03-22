#ifndef TH3_H
#define TH3_H
#include <petsc.h>
#include <cmath>
class TH3
{
private:
	int n1,n2,n3,n4,n5,n6,n7,n8,n9,n21;
	PetscScalar *THCVal,*THCVbl,*THCVdxl;
    PetscScalar *THCVdyl,*THCVvl,*THd;
    PetscScalar *THee,*THg,*THl;
    PetscInt  *mtth;
    PetscScalar *Kseff,*Kf,*ROUf,*Alpha,*KK;
    PetscScalar *Kseffx,*Kseffy,*Kx,*Kfx,*ROUfx,*Ky,*ROUfy,*Kfy;
    PetscScalar *ROUs,*Cps;
public:
	TH3();
	~TH3();
	PetscScalar *THP,*THTf,*THTs,*THU,*THV,*THQfth;
	void setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9,PetscScalar *array10);
	void setTHinitial(PetscScalar *array21,PetscScalar *array22,PetscScalar *array23,PetscScalar *array24,PetscScalar *array25,PetscScalar *array26);
	void setparameter(PetscInt i, PetscInt j,PetscScalar Ts, PetscScalar Tf, PetscScalar P, PetscScalar U, PetscScalar V);
	inline PetscScalar max(PetscScalar a1, PetscScalar a2, PetscScalar a3)
	{
		PetscScalar maxvalue1,maxvalue2;

		maxvalue1 = (a1>a2)?a1:a2;
		maxvalue2 = (maxvalue1>a3)?maxvalue1:a3;
        
        return maxvalue2;
	}
	void setboundarycoefficient(PetscInt i, PetscInt j);
	void setCurrentVariables(const PetscScalar *u);
	void setParameterAll();
	void setBoundaryCoefficientAll();
	void setQf(PetscScalar *Qf);
	PetscScalar getTsa1(PetscInt i, PetscInt j);
	PetscScalar getTsa2(PetscInt i, PetscInt j);
	PetscScalar getTsa3(PetscInt i, PetscInt j);
	PetscScalar getTsa4(PetscInt i, PetscInt j);
	PetscScalar getTsa0(PetscInt i, PetscInt j,PetscScalar a1,PetscScalar a2,PetscScalar a3,PetscScalar a4);
	PetscScalar getTsrs(PetscInt i, PetscInt j,PetscScalar a2,PetscScalar a3,PetscScalar a4);
	PetscScalar getTfa1(PetscInt i, PetscInt j);
	PetscScalar getTfa2(PetscInt i, PetscInt j);
	PetscScalar getTfa3(PetscInt i, PetscInt j);
	PetscScalar getTfa4(PetscInt i,PetscInt j);
	PetscScalar getTfa0(PetscInt i, PetscInt j,PetscScalar a1,PetscScalar a2,PetscScalar a3,PetscScalar a4);
	PetscScalar getTfrs(PetscInt i, PetscInt j,PetscScalar a3);
	PetscScalar getUrs(PetscInt i, PetscInt j);
	PetscScalar getVrs(PetscInt i, PetscInt j);
	PetscScalar getPa1(PetscInt i, PetscInt j);
	PetscScalar getPa2(PetscInt i,PetscInt j);
	PetscScalar getPa3(PetscInt i,PetscInt j);
	PetscScalar getPa4(PetscInt i,PetscInt j);
	PetscScalar getPa0(PetscInt i,PetscInt j, PetscScalar a1,PetscScalar a2,PetscScalar a3, PetscScalar a4);
	PetscScalar getPrs(PetscInt i,PetscInt j,PetscScalar a4);
	PetscScalar getTsTimeCoefficient(PetscInt i, PetscInt j);
	PetscScalar getTfTimeCoefficient(PetscInt i, PetscInt j);
	//PetscScalar getnum(PetscInt i, PetscInt j);
};


#endif