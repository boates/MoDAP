#ifndef GUARD_Student_info
#define GUARD_Student_info
////////// MoDAP //////////
// physics.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>

double distance(atom, atom, const lattice&);
std::vector<double> displacement(atom, atom, const lattice&);
double distance_slow(const atom&, const atom&, const lattice&);
int pbc_round(const double&);
double distance_ortho(const atom&, const atom&, const lattice&);
std::vector<double> displacement_ortho(const atom&, const atom&, const lattice&);
double angle_atoms(const atom&, const atom&, const atom&, const lattice&);

#endif





