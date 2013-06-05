#ifndef GUARD_lattice_h
#define GUARD_lattice_h
////////// MoDAP //////////
// lattice.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>

//================================
// lattice class
//--------------------------------
class lattice
{
public:
	// variables
	const double ax, ay, az;
	const double bx, by, bz;
	const double cx, cy, cz;

	// constructors
	lattice();
	lattice(double);
	lattice(double, double, double);
	lattice(double, double, double, double, double, double, double, double, double);

	// methods
	std::vector<double> a() const;
	std::vector<double> b() const;
	std::vector<double> c() const;
	double amag() const;
	double bmag() const;
	double cmag() const;
	std::vector<double> angles() const;
	double volume() const;
	void print() const;
};

//================================

#endif