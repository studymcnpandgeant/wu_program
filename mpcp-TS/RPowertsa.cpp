
#include "RPowertsa.h"
#include <iostream>
#include <petsc.h>
#include <cmath>

RPowertsa::RPowertsa()
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
	n7 = 18720;
	CVvsiga = new PetscScalar [n7];
	n8 = 4680;
	CVv = new PetscScalar [n8];
	n9 = 18720;
	CVvsigfl = new PetscScalar [n9];

	Qf     = new PetscScalar[n3];
	QQQ    = new PetscScalar[n3];
	n10 = 4680*6;
	Cdflux    = new PetscScalar[n10];
	Cdfluxold = new PetscScalar[n10];
	n11 = 4680*4;
	fluxold   = new PetscScalar[n11];
}

RPowertsa::~RPowertsa()
{
    delete []CVal;
	delete []CVbl;
	delete []CVdxl;
	delete []CVdyl;
	delete []CVsigdl;
	delete []CVsigs0l;
	delete []CVsigtl;
	delete []CVvsiga;
	delete []CVv;
	delete []CVvsigfl;
	delete []Qf;
	delete []QQQ;
	delete []Cdflux;
	delete []Cdfluxold;
	delete []fluxold;
}

void RPowertsa::setinitial(PetscScalar *array1,PetscScalar *array2,PetscScalar *array3,PetscScalar *array4,PetscScalar *array5,PetscScalar *array6,PetscScalar *array7,PetscScalar *array8,PetscScalar *array9)
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

	for (int i = 0; i < 18720; i++)
	{
		CVvsiga[i] = 0.0;
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
	for (int i = 0; i < 28080; i++)
	{
		Cdflux[i] = 0.0;
	}
	for (int i = 0; i < 28080; i++)
	{
		Cdfluxold[i] = 0.0;
	}
	for (int i = 0; i < 18720; i++)
	{
		fluxold[i] = 0.0;
	}
	//initial some parameters
	v[0] = 1.0e9;
    v[1] = 5.7e7;
    v[2] = 4.8e6;
    v[3] = 4.3e5;
    keff = 0.999883;
	dt = 0.1;//very important variable 
    beta[0]=4.9588439e-4;
	beta[1]=3.0678140e-3;
	beta[2]=2.6332880e-3;
	beta[3]=5.4992171e-3;
	beta[4]=1.8051409e-3;
	beta[5]=3.5614735e-4;
	lamda[0]=1.2729429e-2;
	lamda[1]=3.1470820e-2;
	lamda[2]=0.11841680;
	lamda[3]=0.31652680;
	lamda[4]=1.40925502;
	lamda[5]=3.76775407;
	oshift = 1.0e-5;
}

PetscScalar RPowertsa::gettsa1(PetscInt i, PetscInt j, PetscInt group)
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

PetscScalar RPowertsa::gettsa2(PetscInt i, PetscInt j, PetscInt group)
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

PetscScalar RPowertsa::gettsa3(PetscInt i, PetscInt j, PetscInt group)
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

PetscScalar RPowertsa::gettsa4(PetscInt i, PetscInt j, PetscInt group)
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

PetscScalar RPowertsa::gettsa0(PetscInt i, PetscInt j, PetscInt group, PetscScalar a1, PetscScalar a2, PetscScalar a3, PetscScalar a4)
{
	PetscInt index1,index2,index3;
	PetscScalar ac0;

    index1 = 45*(j-1)+i-1;
	index2 = 45*104*(group-1)+45*(j-1)+i-1;
	index3 = 45*104*4*(group-1)+45*104*(group-1)+45*(j-1)+i-1;
    ac0 = (CVsigtl[index2]-CVsigs0l[index3])*CVv[index1] - a1 - a2 - a3 - a4;

    return ac0;
}

PetscScalar RPowertsa::gettsrs(PetscInt i, PetscInt j, PetscInt group, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4, PetscScalar keff, PetscReal time)
{
	PetscScalar Qfd,Qsd,CVsd;
    PetscInt index1,index21,index22,index23,index24,index31,index32,index33,index34;
    PetscInt indexold1,indexold2,indexold3,indexold4;
    PetscScalar chi[4],chil;
    PetscScalar c11,beta0,c31,c41,c21;
    PetscInt    gdl,indexcdfluxold,indexfluxold;
	
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

	//dt    = 0.1;//very important variable
	if (time > 0.0)
	{
	
	beta0 = 1.3857492e-2;
	//beta[0]=4.9588439e-4;
	//beta[1]=3.0678140e-3;
	//beta[2]=2.6332880e-3;
	//beta[3]=5.4992171e-3;
	//beta[4]=1.8051409e-3;
	//beta[5]=3.5614735e-4;
	//lamda[0]=1.2729429e-2;
	//lamda[1]=3.1470820e-2;
	//lamda[2]=0.11841680;
	//lamda[3]=0.31652680;
	//lamda[4]=1.40925502;
	//lamda[5]=3.76775407;
	c11 = (flux1*CVvsigfl[index21]+flux2*CVvsigfl[index22]+flux3*CVvsigfl[index23]+flux4*CVvsigfl[index24])*beta0*chil;
	c31 = flux1*CVvsigfl[index21]+flux2*CVvsigfl[index22]+flux3*CVvsigfl[index23]+flux4*CVvsigfl[index24];
    indexold1 = 45*104*(1-1)+45*(j-1)+i-1;
	indexold2 = 45*104*(2-1)+45*(j-1)+i-1;
	indexold3 = 45*104*(3-1)+45*(j-1)+i-1;
	indexold4 = 45*104*(4-1)+45*(j-1)+i-1;
	c41 = fluxold[indexold1]*CVvsigfl[index21]+fluxold[indexold2]*CVvsigfl[index22]+fluxold[indexold3]*CVvsigfl[index23]+fluxold[indexold4]*CVvsigfl[index24];
	c21 = 0.0;
	for (gdl = 0; gdl < 6; gdl++)
	{
		indexcdfluxold = 45*104*gdl + 45*(j-1) + i-1;
		c21 = c21+ beta[gdl]*(((1-(1-exp(-lamda[gdl]*dt))/(lamda[gdl]*dt))*c31)+((1-exp(-lamda[gdl]*dt))/(lamda[gdl]*dt)-exp(-lamda[gdl]*dt))*c41)
		      +lamda[gdl]*Cdfluxold[indexcdfluxold]*exp(-lamda[gdl]*dt);
	}
	c21 = c21*chil;
	indexfluxold = 45*104*(group-1)+45*(j-1)+i-1;
	CVsd = CVsd + CVv[index1] * (0.0 - c11 + c21);
    }

	return CVsd ;
}

double RPowertsa::getnum()
{
    PetscInt ii=20,jj=40,gg=3;
    PetscInt index23;	
     
	index23 = 45*104*(gg-1)+45*(jj-1)+ii-1;
	return CVvsigfl[index23] ;
}

PetscScalar RPowertsa::getsigf(PetscInt i, PetscInt j, PetscInt group)
{
	PetscInt index;

	index = 45*104*(group-1)+45*(j-1)+i-1;
	return CVvsigfl[index];
}

PetscScalar RPowertsa::getcvv(PetscInt i, PetscInt j)
{
	PetscInt index;

	index = 45*(j-1)+i-1;
	return CVv[index];
}

void RPowertsa::iniQf()
{
	for (int i = 0; i < 4680; i++)
	{
		Qf[i] = 0.0;
	}
	for (int i = 0; i < 28080; i++)
	{
		Cdflux[i] = 0.0;
	}
}

void RPowertsa::addQQQ()
{
	for (int i = 0; i < 4680; i++)
	{
		Qf[i] = Qf[i] + QQQ[i];
	}
}

void RPowertsa::perturbation(PetscReal time)
{
	PetscInt     j,i,gl,gll,indexs0,index3,n;
	PetscScalar  c11;


	n = 45;
	for (j = 1; j < 104+1; j++)
    {
      for (i = 1; i < 45+1; i++)
      {
        for (gl = 1; gl < 5 ; gl++)
        {
          index3 = 45*104*(gl-1) + (j-1)*n + i - 1 ;
          c11 = 0.0;
          for (gll = 1; gll < 5 ; gll++)
          {
          	indexs0 = 45*104*4*(gll-1)+45*104*(gl-1)+45*(j-1)+i-1;
          	c11 = c11 + CVsigs0l[indexs0];
          }
          CVvsiga[index3] = CVsigtl[index3] - c11;
          if (i>=25 && i<=28 && j>=1 && j<=87 && gl==4)
           {
           	CVvsiga[index3] = CVvsiga[index3] * 0.9058 + 0.0* time;  	
           } 
          CVsigtl[index3] = CVvsiga[index3] + c11;
        }
      }
    }
}

PetscScalar RPowertsa::getdflux(PetscInt i, PetscInt j, PetscInt gdl, PetscScalar flux1, PetscScalar flux2, PetscScalar flux3, PetscScalar flux4)
{
	PetscScalar Qfd,dflux;
	PetscInt    index21,index22,index23,index24;
	
	index21 = 45*104*(1-1)+45*(j-1)+i-1;
	index22 = 45*104*(2-1)+45*(j-1)+i-1;
	index23 = 45*104*(3-1)+45*(j-1)+i-1;
	index24 = 45*104*(4-1)+45*(j-1)+i-1;
	Qfd = flux1*CVvsigfl[index21] + flux2*CVvsigfl[index22] + flux3*CVvsigfl[index23] + flux4*CVvsigfl[index24];
	dflux = beta[gdl-1]*Qfd/lamda[gdl-1];

	return dflux;
}