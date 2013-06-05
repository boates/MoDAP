#ifndef GUARD_atom_h
#define GUARD_atom_h
////////// MoDAP //////////
// atom.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include "lattice.h"

//================================
// atom class
//--------------------------------
class atom
{
public:
	// variables
	const double a, b, c;   // reduced coords
	const std::string type; // atom type

	// constructors
    atom();
    atom(double, double, double);
    atom(std::string, double, double, double);
	
	// destructor
//	~atom();
	
	// methods
	std::string type() const;
	double x(const lattice) const;
	double y(const lattice) const;
	double z(const lattice) const;
	std::vector<double> rRed() const;
	std::vector<double> rCart(const lattice) const;
	void print() const;
	
	// operators
	bool operator==(const atom&) const;
};
//================================

#endif