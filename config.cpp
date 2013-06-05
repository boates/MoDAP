////////// MoDAP //////////
// config.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include "lattice.h"
#include "physics.h"
#include "utils.h"
#include "config.h"
using namespace std;

//============================================================
// members of config class (constructors, methods, etc.)
//------------------------------------------------------------

//== constructors ==//

config::config()
{
    nAtom = 0;  // total number of atoms
    nType = 0;  // number of types of atoms
    nBond = 0;  // total number of bonds
    nMole = 0;  // number of molecules
	tStep = 0;  // timestep
	lat   = null;
}

//== methods ==//

iterator config::begin() const
{
	return this->atoms.begin();
}

iterator config::end() const
{
	return this->atoms.end();
}

// number of atoms
int config::natoms() const
{
//	int N = 0; // total number of atoms
	// loop over all atom types
//	for (map<string, int>::iterator it = nAtomType.begin(); it != nAtomType.end(); ++it) {
//		N += it->second; // sum the number of each atom type
//	}
//	return N;
	return nAtom;
}

// number of atoms of a given type
int config::natoms(const string& typat) const
{
	return nAtomType[typat];
/*
    // number of atoms in cfg of type typat //
	int count = 0;
	
	// define iterator type for atoms vector
	typedef vector<atom>::const_iterator iter;
	
    // loop over atoms in cfg and count the ones of desired type
    for (iter i=this->atoms.begin();i!=this->atoms.end();++i){
	    if (i->type == typat){
			count += 1;
		}
	}
	return count;
*/
}

// number of bonds
int config::nbonds() const
{
	return nBond;
}

// number of different atom types
int config::ntypes() const
{
	return nType;
}

// number of molecules
int config::nmoles() const
{
	return nMole;
}

//== insertion methods ==//

void config::insert(atom a)
{
    // append to nodes vector
	this->atoms.push_back(a);

	// update nAtomType map (initialize key for new species)
	if (nAtomType.find(t) == nAtomType.end()) {
	    nAtomType[t] = 1;
    } else {
        nAtomType[t] += 1;
    }
    // update total number of atoms
	nAtom += 1;
}

void config::insert(molecule m)
{
    // append to molecules vector
	this->moles.push_back(m);
}


//== I/O methods ==//

int config::trjRead(ifstream& trj)
{
    //// read in a single trj configuration            ////
    //// returns 0 if successful config read, 1 if eof ////

    // read in lattice info
	double ax, ay, az, bx, by, bz, cx, cy, cz;
	trj >> ax >> ay >> az;
	trj >> bx >> by >> bz;
	trj >> cx >> cy >> cz;
//	lattice tmp_lat (ax, ay, az, bx, by, bz, cx, cy, cz);
//	this->lat = tmp_lat;
	this->lat = new lattice tmp_lat(ax, ay, az, bx, by, bz, cx, cy, cz);	
	
	// exit if eof and return 1
	if (trj.eof()) { return 1; }

	// clear remaining white-space from previous line
	string line;
	getline(trj, line);

    // read natoms line
	getline(trj, line);
	
	// split the line into a vector<string> of natoms 
	vector<string> natoms_vstr = split(line);
	
	// convert to vector of ints
	vstr_to_vnum(natoms_vstr, new vector<int> natoms_vec);

	// get total natom as sum
	int natom_sum = 0;
	for (vector<int>::const_iterator i=natoms_vec.begin(); i!=natoms_vec.end(); ++i) {
		natom_sum += *i;
	}

	// read in the atom types and reduced coordinates
	for (int i=0; i!=natom_sum; ++i){
		
		// initialize atom object and read in its data
		string t; double aRed, bRed, cRed;
		trj >> t >> aRed >> bRed >> cRed;
		atom at(t, aRed, bRed, cRed);

		// add atom to the config object
		this->insert(at);
	}
	return 0;
}

void config::trjWrite(ofstream& trj)
{
    //// write a single trj configuration ////

    // write lattice vectors in natom and tstep
	trj << this->lat.ax << " " << this->lat.ay << " " << this->lat.az << endl;
	trj << this->lat.bx << " " << this->lat.by << " " << this->lat.bz << endl;
	trj << this->lat.cx << " " << this->lat.cy << " " << this->lat.cz << endl;

	// write each natom
	typedef map<string, int>::iterator iter;
	for (iter it=nAtomType.begin(); it!=nAtomType.end(); ++it) {
        trj << it->second << " ";
    }
	trj << endl;

	// write current atomic configuration
	typedef config::iterator iter;
	for (iter n=this->begin(); n!=this->end(); ++n){
		// write each atom's type and reduced coordinates
		trj << n->at->type << " " << n->at->a << " " << n->at->b << " " << n->at->c << endl;
	}
}
/*
void config::xyzRead(ifstream& xyz)
{
    //// read in a single xyz configuration ////

    // read in natom and tstep
//	int natom_total;
	xyz >> new int natom_total;
	xyz >> this->tstep;
	
	// read in current atomic configuration
	for (int i=0; i!=natom_total; ++i){
		
		// initialize atom object and read in its data
		// initialize atom object and read in its data
		string t; double aRed, bRed, cRed;
		trj >> t >> aRed >> bRed >> cRed;
		atom at(t, aRed, bRed, cRed);
		
		atom at;
		xyz >> at.type >> at.x >> at.y >> at.z;
		at.index = i+1; // atom indexing starts at 1!
		
		// add atom to the config object
		*this += at;
	}
}
*/
void config::xyzWrite(ofstream& xyz)
{
    //// write a single xyz configuration ////

    // read in natom and tstep
	xyz << this->nAtom << endl;
	xyz << this->tStep << endl;

	// write current atomic configuration
	typedef config::iterator iter;
	for (iter n=this->begin(); n!=this->end(); ++n){
		// write each atom's type and cartesian coordinates
		xyz << n->at->type << " " << n->at->x << " " << n->at->y << " " << n->at->z << endl;
	}	
}
/*
int config::poscarRead(ifstream& poscar)
{
    //// read in a poscar configuration ////
    //// either vasp4 or vasp5          ////
    //// return version number (4 or 5) ////

    // version number to return
	int version = 0;

    // read header line
	string header; getline(poscar,header);

    // initialize parameters
	double alat, ax, ay, az, bx, by, bz, cx, cy, cz; string line;

    // read alat
	getline(poscar,line); str_to_num(line, alat);
	
    // read in lattice vectors
	getline(poscar,line);
	vector<double> avec; vstr_to_vnum( split(line), avec);
	ax = avec[0]*alat; ay = avec[1]*alat; az = avec[2]*alat;
	
	getline(poscar,line);
	vector<double> bvec; vstr_to_vnum( split(line), bvec);
	bx = bvec[0]*alat; by = bvec[1]*alat; bz = bvec[2]*alat;
	
	getline(poscar,line);
	vector<double> cvec; vstr_to_vnum( split(line), cvec);
	cx = cvec[0]*alat; cy = cvec[1]*alat; cz = cvec[2]*alat;
	
    // build lattice
//	lattice tmp_lat (ax, ay, az, bx, by, bz, cx, cy, cz);
//	this->lat = tmp_lat;
	this->lat = new lattice tmp_lat(ax, ay, az, bx, by, bz, cx, cy, cz);
	
    // read typat(version 5) or natoms(version 4) line
	string line1;
	getline(poscar, line1);

    // read natoms(version 5) or Direct(version 4) line
	string line2;
	getline(poscar, line2);

	// read Direct(version 5) or first atom's reduced coordinates(version 4)
	string line3;
	getline(poscar,line3);
	
	// line3 will be the Direct line (length=1) if vasp version 5
	if (split(line3).size() == 1) {

		// version 5
		version = 5;

		// assign the atom types from line1
		this->types = split(line1);

	    // assigna natoms from line2
        vstr_to_vnum(split(line2), this->natoms);

	} else {
	    // version 4
		version = 4;
		
	    // assign natoms from line1
        vstr_to_vnum(split(line1), this->natoms);

	    // give the atoms generic types: A1, A2, A3, ...
		for (vector<int>::size_type j=0; j!=this->natoms.size(); ++j) {
			string count; num_to_str(j+1, count);
			this->types.push_back("A"+count);
	    }
	}
	
	// get total natom as sum
	for (vector<int>::const_iterator i=this->natoms.begin();i!=this->natoms.end();++i) {
		this->natom += *i;
	}

    // create natsum vector for future atom type naming convenience
	vector<int> natsum = this->natoms;
	for (vector<int>::size_type i=1; i!=natsum.size(); ++i) { 
		natsum[i] += natsum[i-1];
	}

	// read in the atom types and reduced coordinates
	vector<double> red, empty;
	for (int i=0;i!=this->natom;++i){
		
		/// if version 4, the first atom was already read into line3
        if (i == 0 && version == 4) {
    		vstr_to_vnum( split(line3), red);
        } else {
    		// initialize atom object and read in its data
	    	getline(poscar,line); vstr_to_vnum( split(line), red);
    	}

        // initialize atom and assign its reduced coordinates
		atom at; at.a = red[0]; at.b = red[1]; at.c = red[2];		

		// atom indexing starts at 1!
		at.index = i+1;
		
		// assign the correct atom type
		int j = 0;
		while ( (i-natsum[j]) >= 0) { ++j; }
		at.type = this->types[j];
		
		// reset red to empty vector<double>
		red = empty;
		
		// set cartesian coordinates
		at.set_cart(this->lat);
		
		// add atom to the config object
		*this += at;
	}
	
	// set ntypes
	this->ntypes = this->get_ntypes();
	
	return version;
}

int config::xdatcarRead(ifstream& xdatcar, const config& cfg)
{
    //// read in a xdatcar configuration  ////
    //// assuming header already stripped ////

    //// reads blank line BEFORE reading coords ////

    //// returns 0 if successful config read, 1 if eof ////

    // read blank line
	string line; getline(xdatcar,line);
	
    //// cfg should contain lattice/natoms/types info    ////
    //// i.e. information from a POSCAR file for example ////
	this->lat    = cfg.lat;
	this->natom  = cfg.natom;
	this->natoms = cfg.natoms;
	this->ntypes = cfg.ntypes;
	this->types  = cfg.types;
	
	// if end of XDATCAR file, exit function
	if (xdatcar.eof()) { return 1; }
	
    // create natsum vector for future atom type naming convenience
	vector<int> natsum = this->natoms;
	for (vector<int>::size_type i=1; i!=natsum.size(); ++i) { 
		natsum[i] += natsum[i-1];
	}

	// read in the atom types and reduced coordinates
	vector<double> red, empty;
	for (int i=0; i!=this->natom; ++i){
		
		// initialize atom object and read in its data
		getline(xdatcar,line);
		vstr_to_vnum(split(line), red);

        // initialize atom and assign its reduced coordinates
		atom at; at.a = red[0]; at.b = red[1]; at.c = red[2];		

		// atom indexing starts at 1!
		at.index = i+1;
		
		// assign the correct atom type
		int j = 0;
		while ( (i-natsum[j]) >= 0) { ++j; }
		at.type = this->types[j];
		
		// reset red to empty vector<double>
		red = empty;
		
		// set cartesian coordinates
		at.set_cart(this->lat);
		
		// add atom to the config object
		*this += at;
	}
	return 0;
}

*/

void config::wrap()
{
    //// wrap the reduced coordinates of config ////
    //// and update the cartesian ones          ////

    // loop over all atoms in config
	typedef config::iterator iter;
	for (iter at=this->begin(); at!=this->end(); ++at) {
		
        // keep in mind, floor(-4.021) = -5
		double anew = at->a - floor(at->a); // + (1-sgn(at->a))/2;
		double bnew = at->b - floor(at->b); // + (1-sgn(at->b))/2;
		double cnew = at->c - floor(at->c); // + (1-sgn(at->c))/2;

        // replace old atom with new wrapped one
		*at = new atom(at->type, anew, bnew, cnew);
	}

}

config config::supercell(const int& K, const int& L, const int& M) const
{
    //// return a K x L x M supercell of cfg ////
    //// does NOT modify original cfg object ////
	config super;
	vector<double> a_sup = doubleXvec(this->lat.a(), K);
	vector<double> b_sup = doubleXvec(this->lat.b(), L);
	vector<double> c_sup = doubleXvec(this->lat.c(), M);
//	lattice slat(a_sup, b_sup, c_sup);
//	super.lat = slat;
	super.lat = new lattice slat(a_sup, b_sup, c_sup);

    // loop over all atoms in config
	typedef config::iterator iter;
	atom q; int i = 1;
	for (iter at=this->begin(); at!=this->end(); ++at) {
		
		// assign short-hand variable names for current atom
		string t = at->type;
		double a = at->a;
		double b = at->b;
		double c = at->c;
		
		// adjust and append current atom to supercell
		super += new atom(t, a/K, b/L, c/M);
		
		// loop and create supercell duplicates
		for (int k=1; k!=K; ++k){
			super += new atom(t, (a+k)/K, b/L, c/M);
		}
		for (int l=1; l!=L; ++l){
			super += new atom(t, a/K, (b+l)/L, c/M);
		}
		for (int m=1; m!=M; ++m){
			super += new atom(t, a/K, b/L, (c+m)/M);
		}
		for (int k=1; k!=K; ++k){
			for (int l=1; l!=L; ++l){
    			super += new atom(t, (a+k)/K, (b+l)/L, c/M);
			}
		}
		for (int k=1; k!=K; ++k){
			for (int m=1; m!=M; ++m){
    			super += new atom(t,(a+k)/K, b/L, (c+m)/M);
			}
		}
		for (int l=1; l!=L; ++l){
			for (int m=1; m!=M; ++m){
    			super += new atom(t, a/K, (b+l)/L, (c+m)/M);
			}
		}
		for (int k=1; k!=K; ++k){
    		for (int l=1; l!=L; ++l){
	    		for (int m=1; m!=M; ++m){
    	    		super += new atom(t, (a+k)/K, (b+l)/L, (c+m)/M);
				}
			}
		}
	}
    return super;	
}

list<double> config::distances(const string& typ1, const string& typ2)
{	
	// create list for atomic distances
	list<double> d;
	
	// define iterator for vector<atom>
	typedef config::iterator iter;

    // loop over all atoms for double loop
	for (iter a1=this->begin(); a1!=this->end(); ++a1) {
		
		// loop over atoms up to a1 (to avoid double counting)
		for (iter a2=this->begin(); a2!=a1; ++a2) {
			
			// if a1 & a2 are of the desired atom pair
			if ( (a1->type == typ1 and a2->type == typ2) or
			     (a1->type == typ2 and a2->type == typ1)) {
				
    			// compute a1-a2 distance and append to d list
	    		d.push_back( distance(*a1, *a2, this->lat) );
	        }
		}
	}
	return d;	
}

void config::molecules(map<string,double>& cutoffs)
{
    //// find all molecules in a configuration ////

	// double loop over atoms in configuration
	typedef config::const_iterator iter;
	for (iter a1=this->begin(); a1!=this->end(); ++a1) {
			
		// create new molecule beginning with a1
		molecule m;
		m += *a1;
		
		// second loop over atoms in configuration
   		for (iter a2=this->begin(); a2!=this->end(); ++a2) {
	
            // compute distance between a1 and a2
			double d = distance(*a1, *a2, this->lat);

            // create atom pair names for cutoffs
			string p1 = a1->type+"-"+a2->type;
			string p2 = a2->type+"-"+a1->type;
			
            // determine relevant cutoff (zero if not in cutoffs)
			double cutoff = max(cutoffs[p1],cutoffs[p2]);

			// if a1 within cutoff distance of a2 (and not the same atom)
			if (d <= cutoff && a1!=a2){
				
				// add a2 to a1's molecule
				m += *a2;
				
				// add an a1-a2 bond to the molecule
				// a1 is the first element in m.atoms
				// a2 is the last element in m.atoms
				// d is the bond length
				m += bond(m.atoms[0],m.atoms[m.atoms.size()-1],d);
			}
      	}
    	// add to molecules list
		(*this) += m;
	}

	// merge all molecules with common atoms
	merge_molecules();
	
	// update object nmol and nbond values
	this->nmol = this->get_nmol();
	this->nbond = this->get_nbond();
}

void config::merge_molecules() {

    //// merge molecules that share common atom(s) recursively ////

	// define type for looping
	typedef list<molecule>::iterator iter;

	// loop over molecules in mols
	for (iter i=this->mols.begin();i!=this->mols.end();++i){

		// loop over molecules in mols up to i (i!=j imposed here)
		iter j = this->mols.begin();
		while (j != i) {

	        // check to see if i and j molecules share common atoms
			if ( common_value(i->atoms,j->atoms) ){
				
				// merge molecule j into i
				*i += *j;

				// erase j from mols list and adjust iterator
				// subtlty here: once erasing j, i iterator not affected?
				j = this->mols.erase(j);
				
			} else
				++j;
	    }
	}

    // make sure sum of atoms in molecules = total natom
	// loop over molecules list to sum all atoms in all molecules
	int asum = 0;
	for (iter a=this->mols.begin();a!=this->mols.end();++a){
		asum += (a->atoms).size();
	}
	// if asum is inconsistent with natom, print info and throw error
	if (asum != this->natom){
		cout << "total atoms in mols = " << asum << ", natom = " << this->natom << endl;
		cout << "nmols = " << this->get_nmol() << endl;
		throw invalid_argument("total number of atoms in molecules != system natom");
	}
};

//== operators ==//

config& config::operator+=(const atom& a)
{
//	if (!this.contains(a)) this->insert(a);
	this->insert(a);
	return *this;
/*
	this->atoms.push_back(a);
	return *this;
*/
}
/*
config& config::operator+=(const bond& b)
{
	this->bonds.push_back(b);
	return *this;
}
*/
config& config::operator+=(const molecule& m)
{
	this->insert(m);
//	this->moles.push_back(m);
	return *this;
}

//============================================================











