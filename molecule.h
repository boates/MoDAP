#ifndef GUARD_molecule_h
#define GUARD_molecule_h
////////// MoDAP //////////
// molecule.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>

//================================
// molecule class
//--------------------------------
class molecule {
public:
    // constructors
	molecule();
	
	// fields
	std::string label;
	std::vector<atom> atoms;
	std::vector<bond> bonds;
	
	// methods
	double natom() const;
	double nbond() const;
	std::string species();
	std::string ID() const;
	void print();
	void print(std::string);

	// operators
	bool operator==(const molecule& m) const;
	bool operator> (const molecule& m) const;
	bool operator>=(const molecule& m) const;
	bool operator< (const molecule& m) const;
	bool operator<=(const molecule& m) const;
	molecule  operator+(const molecule& m);
	molecule& operator+=(const molecule& m);
	molecule  operator+(const atom& a);
	molecule& operator+=(const atom& a);
	molecule  operator+(const bond& b);
	molecule& operator+=(const bond& b);
};
//================================

#endif