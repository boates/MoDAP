#ifndef GUARD_bond_h
#define GUARD_bond_h
////////// MoDAP //////////
// bond.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>

//================================
// bond class
//--------------------------------
class bond
{
public:
	// variables
//	const atom at1, at2;
	const int at1, at2;
	const double length;
	
	// constructor
	bond(int, int, double);
//	bond(atom, atom, double);
	
	// methods
	void print() const;
	void print(std::string) const;
	
	// operators
	bool operator==(const bond&) const;
	bool operator> (const bond&) const;
	bool operator>=(const bond&) const;
	bool operator< (const bond&) const;
	bool operator<=(const bond&) const;
};
//================================

#endif