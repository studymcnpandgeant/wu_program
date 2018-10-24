
#include "Point.h"
#include <iostream>
#include <petsc.h>

Point::Point()
{
	dx = 0.0;
}

Point::~Point()
{
   std::cout << "\n exittttttttttttttttt" ;
}

void Point::setdx(double dxx)
{
	dx = dxx ; 
}

double Point::getdx()
{
	return dx ;
}