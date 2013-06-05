////////// MoDAP //////////
// atom.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <iostream>
#include "lattice.h"
#include "atom.h"
using namespace std;

//============================================================
// members of atom class (constructors, methods, etc.)
//------------------------------------------------------------

//== constructors ==//

atom::atom() : typ("X"), a(0), b(0), c(0) {}

atom::atom(double l, double m, double n) : typ("X"), a(l), b(m), c(n) {}

atom::atom(std::string s, double l, double m, double n) :
                                           typ(s), a(l), b(m), c(n) {}

/*
//== destructor ==//

atom::~atom()
{
	delete &typ;
	delete &ind;
	delete &ai, &bi, &ci;
}
*/

//== methods ==//

// cartesian x-coord
double atom::x(const lattice lat) const
{
	return this->a*lat.ax + this->b*lat.bx + this->c*lat.cx;
}

// cartesian y-coord
double atom::y(const lattice lat) const
{
	return this->a*lat.ay + this->b*lat.by + this->c*lat.cy;
}

// cartesian z-coord
double atom::z(const lattice lat) const
{
	return this->a*lat.az + this->b*lat.bz + this->c*lat.cz;
}

// reduced coord vector
vector<double> atom::rRed() const
{
	double array[3] = {this->a, this->b, this->c};
	vector<double> rvec(array, array+3);
	return rvec;
}

// cartesian coord vector
vector<double> atom::rCart(const lattice lat) const
{
	double array[3] = {this->x(lat), this->y(lat), this->z(lat)};
	vector<double> rvec(array, array+3);
	return rvec;
}

// print atom info
void atom::print() const
{
	cout << "[" << type << ", " << a << ", " << b << ", " << c << "]";
}

//== operators ==//

bool atom::operator==(const atom& at) const
{
	return this->type==at.type && this->a==at.a && this->b==at.b && this->c==at.c;
}

//============================================================
/*
int main()
{
    // initialize atom and print it
	atom at("N", 0.123, 0.456, 0.789);
	at.print(); cout << endl;
	return 0;
}
*/
