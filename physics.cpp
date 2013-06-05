////////// MoDAP //////////
// physics.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include "config.h"
#include "lattice.h"
#include "utils.h"
#include "physics.h"
using namespace std;

// compute the distance between two atoms (arbitrary lattice)
double distance(atom at1, atom at2, const lattice& lat){
	
	// return distance between two atoms using the minimum image convention //

	// make sure at1 and at2 reduced coordinates are wrapped
	at1.a += -floor(at1.a); at2.a += -floor(at2.a);
	at1.b += -floor(at1.b); at2.b += -floor(at2.b);
	at1.c += -floor(at1.c); at2.c += -floor(at2.c);
	
	// get difference in reduced a coordinate within 0.5
	double da = at1.a - at2.a;
	if (da > 0.5) {at2.a += 1; da += -1;}
	else if (da < -0.5) {at1.a += 1; da += 1;}

	// get difference in reduced b coordinate within 0.5
	double db = at1.b - at2.b;
	if (db > 0.5) {at2.b += 1; db += -1;}
	else if (db < -0.5) {at1.b += 1; db += 1;}

	// get difference in reduced c coordinate within 0.5
	double dc = at1.c - at2.c;
	if (dc > 0.5) {at2.c += 1; dc += -1;}
	else if (dc < -0.5) {at1.c += 1; dc += 1;}
	
	// cartesian vector of the reduced differences
	vector<double> dr_cart;
	dr_cart.push_back( da*lat.ax + db*lat.bx + dc*lat.cx );
	dr_cart.push_back( da*lat.ay + db*lat.by + dc*lat.cy );
	dr_cart.push_back( da*lat.az + db*lat.bz + dc*lat.cz );
	
	// return magnitude of "difference" cartesian vector (= distance)
	return magnitude( dr_cart );
}

// compute the displacement vector between two atoms (arbitrary lattice)
vector<double> displacement(atom at1, atom at2, const lattice& lat){
	
	// return dispalcement vector between two atoms w/ minimum image convention //

	// make sure at1 and at2 reduced coordinates are wrapped
	at1.a += -floor(at1.a); at2.a += -floor(at2.a);
	at1.b += -floor(at1.b); at2.b += -floor(at2.b);
	at1.c += -floor(at1.c); at2.c += -floor(at2.c);
	
	// get difference in reduced a coordinate within 0.5
	double da = at1.a - at2.a;
	if (da > 0.5) {at2.a += 1; da += -1;}
	else if (da < -0.5) {at1.a += 1; da += 1;}

	// get difference in reduced b coordinate within 0.5
	double db = at1.b - at2.b;
	if (db > 0.5) {at2.b += 1; db += -1;}
	else if (db < -0.5) {at1.b += 1; db += 1;}

	// get difference in reduced c coordinate within 0.5
	double dc = at1.c - at2.c;
	if (dc > 0.5) {at2.c += 1; dc += -1;}
	else if (dc < -0.5) {at1.c += 1; dc += 1;}
	
	// cartesian vector of the reduced differences
	vector<double> dr_cart;
	dr_cart.push_back( da*lat.ax + db*lat.bx + dc*lat.cx );
	dr_cart.push_back( da*lat.ay + db*lat.by + dc*lat.cy );
	dr_cart.push_back( da*lat.az + db*lat.bz + dc*lat.cz );
	
	// "difference" cartesian vector (= displacement)
	return dr_cart;
}

// compute the distance between two atoms via loop over supercell (arbitrary lattice)
double distance_slow(const atom& at1, const atom& at2, const lattice& lat){
	
	// return distance between two atoms using the minimum image convention //

    // make temporary config file with only atom 1
	config c; c += at1; c.wrap();
	c.lat = lat; c.natom = c.get_natom(); // natom = 1
	
	// make 3x3x3 supercell
	config s = c.supercell(3,3,3);

	// define iterator for vector<atom> and distance variable
	typedef vector<atom>::const_iterator iter;
	double d; double dmin = 100.0; // something large
	
	// shift at2 to middle box of 3x3x3 supercell
	vector<double> at2rcart;
	at2rcart.push_back(at2.x + lat.ax + lat.bx + lat.cx);
	at2rcart.push_back(at2.y + lat.ay + lat.by + lat.cy);
	at2rcart.push_back(at2.z + lat.az + lat.bz + lat.cz);
	
	// loop over all at1's in 3x3x3 supercell
	for (iter a1=s.atoms.begin();a1!=s.atoms.end();++a1){		

		// compute distance between a2 and all at1 replicas
		// as magnitude of vector difference in cartesian coordinates
		d = magnitude( vec_subtract( at2rcart, a1->rcart(s.lat) ) );
		
		// check if it is the smallest distance so far, if so set dmin
		if (d < dmin){
			dmin = d;
		}
	}
	return dmin;
}

// pbc_round function
int pbc_round(const double& v){
	int i = v;
	if (abs(v)-i >= 0.5){
		if (v > 0) {i+=1;}
		if (v < 0) {i-=1;}
	}
	return i;
}

// compute the distance between two atoms (orthrhombic cell)
double distance_ortho(const atom& at1, const atom& at2, const lattice& lat){
	// compute cartesian distance
	double dx = at1.x - at2.x;
	double dy = at1.y - at2.y;
	double dz = at1.z - at2.z;
	dx -= lat.amag * pbc_round( dx/lat.amag );
	dy -= lat.bmag * pbc_round( dy/lat.bmag );
	dz -= lat.cmag * pbc_round( dz/lat.cmag );
	return sqrt( pow(dx,2)+pow(dy,2)+pow(dz,2) );
}

// compute the displacement between two atoms (orthrhombic cell)
vector<double> displacement_ortho(const atom& at1, const atom& at2, const lattice& lat){
	// compute cartesian displacement vector
	double dx = at1.x - at2.x;
	double dy = at1.y - at2.y;
	double dz = at1.z - at2.z;
	dx -= lat.amag * pbc_round( dx/lat.amag );
	dy -= lat.bmag * pbc_round( dy/lat.bmag );
	dz -= lat.cmag * pbc_round( dz/lat.cmag );
	vector<double> dr;
	dr.push_back(dx); dr.push_back(dy); dr.push_back(dz);
	return dr;
}

// compute the a-b-c angle (in degrees!)
double angle_atoms(const atom& a, const atom& b, const atom& c, const lattice& lat){
	return angle( displacement_ortho(b,a,lat), displacement_ortho(b,c,lat) );
}






