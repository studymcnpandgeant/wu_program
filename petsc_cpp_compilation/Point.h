#ifndef POINT_H
#define POINT_H
#include <petsc.h>
class Point
{
private:
	PetscReal dx;
public:
	Point();
	~Point();
	void setdx( double dxx);
	double getdx();
};


#endif