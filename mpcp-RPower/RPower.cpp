
#include "RPower.h"
#include <iostream>
#include <petsc.h>

RPower::RPower()
{
	n1 = 4784;
	CVal = new PetscScalar [n1];
	n2 = 4725;
	CVbl = new PetscScalar [n2];
	n3 = 4680;
	CVdxl = new PetscScalar [n3];
	n4 = 4680;
	CVdyl = new PetscScalar [n4];
	n5 = 18720;
	CVsigdl = new PetscScalar [n5];
	n6 = 74880;
	CVsigs0l = new PetscScalar [n6];
	n7 = 18720;
	CVsigtl = new PetscScalar [n7];
	n8 = 4680;
	CVv = new PetscScalar [n8];
	n9 = 18720;
	CVvsigfl = new PetscScalar [n9];

	Qf = new PetscScalar[n3];
	QQQ = new PetscScalar[n3];
}

RPower::~RPower()
{
    delete []CVal;
	delete []CVbl;
	delete []CVdxl;
	delete []CVdyl;
	delete []CVsigdl;
	delete []CVsigs0l;
	delete []CVsigtl;
	delete []CVv;
	delete []CVvsigfl;
	delete []Qf;
	delete []QQQ;
}

void RPower::setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9)
{

	for (int i = 0; i < 4784; i++)
	{
		CVal[i] = array1[i];
	}

	for (int i = 0; i < 4725; i++)
	{
		CVbl[i] = array2[i];
	}

	for (int i = 0; i < 4680; i++)
	{
		CVdxl[i] = array3[i];
	}

	for (int i = 0; i < 4680; i++)
	{
		CVdyl[i] = array4[i];
	}

	for (int i = 0; i < 18720; i++)
	{
		CVsigdl[i] = array5[i];
	}

	for (int i = 0; i < 74880; i++)
	{
		CVsigs0l[i] = array6[i];
	}
	for (int i = 0; i < 18720; i++)
	{
		CVsigtl[i] = array7[i];
	}

	for (int i = 0; i < 4680; i++)
	{
		CVv[i] = array8[i];
	}

	for (int i = 0; i < 18720; i++)
	{
		CVvsigfl[i] = array9[i];
	}

	for (int i = 0; i < 4680; i++)
	{
		Qf[i] = 0.0;
	}
	
}

PetscScalar RPower::geta1(PetscInt i, PetscInt j, PetscInt group)
{
	PetscScalar v11,v21,c11,ac1;
	PetscInt index1,index2,index3;

    if ( i > 1 )
    {
    	index1 = 45*(j-1)+i-1;
	    v11 = CVv[index1-1];
	    v21 = CVv[index1];
	    index2 = 45*104*(group-1)+45*(j-1)+i-1;
	    c11 = (CVsigdl[index2-1]*v11+CVsigdl[index2]*v21)/(v11+v21);
	    index3 = 46*(j-1)+i-1;
	    ac1 = -2*c11*CVal[index3]/(CVdxl[index1]+CVdxl[index1-1]);

	    return ac1;
    }
    else
    {
    	ac1 = 0.0;

    	return ac1;
    }

}

PetscScalar RPower::geta2(PetscInt i, PetscInt j, PetscInt group)
{
	PetscScalar v11,v21,c11,ac2;
	PetscInt index1,index2,index3;

	if ( i < 45 )
	{
		index1 = 45*(j-1)+i-1;
		v11 = CVv[index1+1];
		v21 = CVv[index1];
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		c11 = (CVsigdl[index2+1]*v11+CVsigdl[index2]*v21)/(v11+v21);
		index3 = 46*(j-1)+i-1;
		ac2 = -2*c11*CVal[index3+1]/(CVdxl[index1]+CVdxl[index1+1]);

		return ac2;
	}
	else
	{
		index1 = 45*(j-1)+i-1;
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		index3 = 46*(j-1)+i-1;
		ac2 = -2*CVsigdl[index2]*CVal[index3+1]/CVdxl[index1];

		return ac2;
	}
}

PetscScalar RPower::geta3(PetscInt i, PetscInt j, PetscInt group)
{
	PetscScalar v11,v21,c11,ac3;
	PetscInt index1,index2,index3;

	if ( j > 1 )
	{
		index1 = 45*(j-1)+i-1;
		v11 = CVv[index1-45];
		v21 = CVv[index1];
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		c11 = (CVsigdl[index2-45]*v11+CVsigdl[index2]*v21)/(v11+v21);
        index3 = 45*(j-1)+i-1;
		ac3 = -2*c11*CVbl[index3]/(CVdyl[index1]+CVdyl[index1-45]);

		return ac3;
	}
	else
	{
		index1 = 45*(j-1)+i-1;
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		ac3 = -2*CVsigdl[index2]*CVbl[index1]/CVdyl[index1];

		return ac3;
	}
}

PetscScalar RPower::geta4(PetscInt i, PetscInt j, PetscInt group)
{
	PetscScalar v11,v21,c11,ac4;
	PetscInt index1,index2,index3;

	if ( j < 104 )
	{
		index1 = 45*(j-1)+i-1;
		v11 = CVv[index1+45];
		v21 = CVv[index1];
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		c11 = (CVsigdl[index2+45]*v11+CVsigdl[index2]*v21)/(v11+v21);
		index3 = 45*(j-1)+i-1;
		ac4 = -2*c11*CVbl[index3+45]/(CVdyl[index1]+CVdyl[index1+45]);

		return ac4;
	}
	else
	{
		index1 = 45*(j-1)+i-1;
		index2 = 45*104*(group-1)+45*(j-1)+i-1;
		ac4 = -2*CVsigdl[index2]*CVbl[index1+45]/CVdyl[index1];

		return ac4;
	}
}

PetscScalar RPower::geta0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4)
{
	PetscInt index1,index2,index3;
	PetscScalar ac0;

    index1 = 45*(j-1)+i-1;
	index2 = 45*104*(group-1)+45*(j-1)+i-1;
	index3 = 45*104*4*(group-1)+45*104*(group-1)+45*(j-1)+i-1;
    ac0 = (CVsigtl[index2]-CVsigs0l[index3])*CVv[index1] - a1 - a2 - a3 - a4;

    return ac0;
}

PetscScalar RPower::getrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff)
{
	PetscScalar Qfd,Qsd,CVsd;
    PetscInt index1,index21,index22,index23,index24,index31,index32,index33,index34;
    PetscScalar chi[4],chil;
	
	chi[0] = 0.984390; chi[1] = 0.015616; chi[2] = 0.67690e-7; chi[3] = 0.0;
	index1 = 45*(j-1)+i-1;
	index21 = 45*104*(1-1)+45*(j-1)+i-1;
	index22 = 45*104*(2-1)+45*(j-1)+i-1;
	index23 = 45*104*(3-1)+45*(j-1)+i-1;
	index24 = 45*104*(4-1)+45*(j-1)+i-1;
	switch(group)
	{
		case 1: chil = chi[0];
		        break;
		case 2: chil = chi[1];
		        break;
		case 3: chil = chi[2];
		        break;
		case 4: chil = chi[3];
		        break;
	}
	Qfd = flux1*CVvsigfl[index21]*chil+flux2*CVvsigfl[index22]*chil+flux3*CVvsigfl[index23]*chil+flux4*CVvsigfl[index24]*chil;
	index31 = 4*45*104*(group-1)+45*104*(1-1)+45*(j-1)+i-1;
	index32 = 4*45*104*(group-1)+45*104*(2-1)+45*(j-1)+i-1;
	index33 = 4*45*104*(group-1)+45*104*(3-1)+45*(j-1)+i-1;
	index34 = 4*45*104*(group-1)+45*104*(4-1)+45*(j-1)+i-1;
	switch(group)
	{
		case 1: Qsd = flux2*CVsigs0l[index32]+flux3*CVsigs0l[index33]+flux4*CVsigs0l[index34];
		        break;
		case 2: Qsd = flux1*CVsigs0l[index31]+flux3*CVsigs0l[index33]+flux4*CVsigs0l[index34];
		        break;
		case 3: Qsd = flux2*CVsigs0l[index32]+flux1*CVsigs0l[index31]+flux4*CVsigs0l[index34];
		        break;
		case 4: Qsd = flux2*CVsigs0l[index32]+flux3*CVsigs0l[index33]+flux1*CVsigs0l[index31];
		        break;
	}
	CVsd = (Qfd+Qsd+Qfd*(1.0-keff)/(keff+0.0000000001))*CVv[index1];

	return CVsd;
}

double RPower::getnum()
{
    PetscInt ii=20,jj=40,gg=3;
    PetscInt index23;	
     
	index23 = 45*104*(gg-1)+45*(jj-1)+ii-1;
	return CVvsigfl[index23] ;
}

PetscScalar RPower::getsigf(PetscInt i, PetscInt j, PetscInt group)
{
	PetscInt index;

	index = 45*104*(group-1)+45*(j-1)+i-1;
	return CVvsigfl[index];
}

PetscScalar RPower::getcvv(PetscInt i, PetscInt j)
{
	PetscInt index;

	index = 45*(j-1)+i-1;
	return CVv[index];
}

void RPower::iniQf()
{
	for (int i = 0; i < 4680; i++)
	{
		Qf[i] = 0.0;
	}
}

void RPower::addQQQ()
{
	for (int i = 0; i < 4680; i++)
	{
		Qf[i] = Qf[i] + QQQ[i];
	}
}