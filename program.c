#include "Point.h"
#include <iostream>

int main()
{
	Point p1;
	double a = 2.34,b;

	p1.setdx(a);
    b = p1.getdx();
    cout << "The double b is " << b ;


	return 0;
}