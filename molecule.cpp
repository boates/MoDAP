////////// MoDAP //////////
// molecule.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "atom.h"
#include "bond.h"
#include "utils.h"
#include "alphabetical.h"
#include "molecule.h"
using namespace std;

//============================================================
// members of molecule class (constructors, methods, etc.)
//------------------------------------------------------------
// constructors
molecule::molecule(){
	label = "";
};
// methods
double molecule::natom() const{
	return atoms.size();	
};
double molecule::nbond() const{
	return bonds.size();	
};
string molecule::species() {
	
	//// return species name as string (i.e. "CO2" or "N2H4") ////
	
	// define map to hold atom types and counts
	// use alphabetical::operator() for key comparions
	//     ---> alphabetical by first character of atom type
	map<string,int,alphabetical> typat;
	
	// define iterator type for atoms vector
	typedef vector<atom>::const_iterator iter;
	
	// loop over atoms in molecule
	for (iter i=this->atoms.begin();i!=this->atoms.end();++i){
		// initialize all detected type's counts to zero
		typat[(*i).type] = 0;
	}
	// loop over atoms in molecule
	for (iter i=this->atoms.begin();i!=this->atoms.end();++i){
		// increase count by one for each type found
		typat[(*i).type] += 1;
	}
	
	// create species string
	string spec = "";
	
	// loop over <key,value> pairs in semi-alphabetized typat map
	for (map<string,int>::iterator i=typat.begin();i!=typat.end();++i){
		// write key (atom type) to species string
		spec += (*i).first;

        // only write amount if > 1
        if ((*i).second > 1) {
    		// write value (amount) to species string
		    ostringstream convert; // stream for int to string conversion
	    	convert << (*i).second;
    		spec += convert.str();
        }
	}
	return spec;
};
string molecule::ID() const{
	
	//// create a unique ID for the molecule based ////
	//// on atom types and their indicies          ////
	
	// define map to hold atom indices and types
	// soft sorted by atom index (map default)
	map<int,string> indices;
	
	// define iterator type for atoms vector
	typedef vector<atom>::const_iterator iter;
	
	// loop over atoms in molecule
	for (iter i=this->atoms.begin();i!=this->atoms.end();++i){
		// key = atom index, value = atom type
		indices[(*i).index] = (*i).type;
	}
	
	// create ID string
	string id = "";
	
	// loop over <key,value> pairs in indices map
	for (map<int,string>::iterator i=indices.begin();i!=indices.end();++i){
		// write value (atom type) to ID string
		id += (*i).second;

   		// write key (index) to ID string
	    ostringstream convert; // stream for int to string conversion
    	convert << (*i).first;
   		id += convert.str();
	}
	return id;
};
void molecule::print() {
	// print molecule
	cout << "#" << endl << "# ====== printing molecule ====== #" << endl;
	
	// print molecule species
    cout << "# species: " << this->species() << endl;

	// print molecule ID
	cout << "# ID: " << this->ID() << endl;
	
	// print molecule label
	cout << "# label: " << label << endl;
	
	// print list of atoms in molecule
	cout << "# atoms (" << (*this).natom() << "):" << endl;
	typedef vector<atom>::const_iterator a_it;
	for (a_it a=this->atoms.begin();a!=this->atoms.end();++a){
		cout << "#  "; (*a).print(); cout << endl;
	}
	
	// print list of bonds in molecule
	cout << "# bonds (" << (*this).nbond() << "):" << endl;
	typedef vector<bond>::const_iterator b_it;
	for (b_it b=this->bonds.begin();b!=this->bonds.end();++b){
		cout << "#  "; (*b).print();
	}
	cout << "# ====== end of molecule ======== #" << endl << "#" << endl;
};
void molecule::print(string str) {
	// print molecule
	cout << "#" << endl << "# === printing molecule; "+str+" only === #" << endl;
	
	// print molecule species
	if (str == "species"){
    	cout << "# species: " << this->species() << endl;
    }
	// print molecule ID
	else if (str == "id" || str == "ID"){
    	cout << "# ID: " << this->ID() << endl;
    }
	// print molecule label
	else if (str == "label"){
    	cout << "# label: " << label << endl;
    }
	// print list of atoms in molecule
	else if (str == "atoms") {
    	cout << "# atoms (" << (*this).natom() << "):" << endl;
	    typedef vector<atom>::const_iterator a_it;
    	for (a_it a=this->atoms.begin();a!=this->atoms.end();++a){
	    	cout << "#  "; (*a).print(); cout << endl;
    	}
    }
	// print list of bonds in molecule
	else if (str == "bonds") {
    	cout << "# bonds (" << (*this).nbond() << "):" << endl;
	    typedef vector<bond>::const_iterator b_it;
    	for (b_it b=this->bonds.begin();b!=this->bonds.end();++b){
	    	cout << "#  "; (*b).print();
    	}
    }
    // print options
    else
		cout << "# options: \"species\", \"ID\", \"label\" \"atoms\", \"bonds\""
    	     << "(all: give no arg)" << endl;
	cout << "# ======== end of molecule ======== #" << endl << "#" << endl;
};
// operators
bool molecule::operator==(const molecule& m) const{
	return this->atoms==m.atoms && this->bonds==m.bonds;
};
bool molecule::operator> (const molecule& m) const{
	return this->natom() > m.natom();
};
bool molecule::operator>=(const molecule& m) const{
	return this->natom() >= m.natom();
};
bool molecule::operator< (const molecule& m) const{
	return this->natom() < m.natom();
};
bool molecule::operator<=(const molecule& m) const{
	return this->natom() <= m.natom();
};
molecule molecule::operator+(const molecule& m) {
	
    // initialize new molecule for the sum
	molecule mol;

    // insert all atoms from 'this' and 'm'	and remove duplicates
	mol.atoms.insert(mol.atoms.end(), this->atoms.begin(), this->atoms.end());
	mol.atoms.insert(mol.atoms.end(), m.atoms.begin(), m.atoms.end());
	remove_duplicates(mol.atoms);

    // insert all atoms from 'this' and 'm'	and remove duplicates
	mol.bonds.insert(mol.bonds.end(), this->bonds.begin(), this->bonds.end());
	mol.bonds.insert(mol.bonds.end(), m.bonds.begin(), m.bonds.end());
	remove_duplicates(mol.bonds);
	
	return mol;
};
molecule& molecule::operator+=(const molecule& m) {

    // insert all of atoms from 'm' into 'this' and remove duplicates
	this->atoms.insert(this->atoms.end(), m.atoms.begin(), m.atoms.end());
	remove_duplicates(this->atoms);

    // insert all of bonds from 'm' into 'this' and remove duplicates
	this->bonds.insert(this->bonds.end(), m.bonds.begin(), m.bonds.end());
	remove_duplicates(this->bonds);

	return *this;
};
molecule molecule::operator+(const atom& a) {
	molecule mol;
	mol.atoms.push_back(a);
	return mol;
};
molecule& molecule::operator+=(const atom& a) {
	this->atoms.push_back(a);
	return *this;
};
molecule molecule::operator+(const bond& b) {
	molecule mol;
	mol.bonds.push_back(b);
	return mol;
};
molecule& molecule::operator+=(const bond& b) {
	this->bonds.push_back(b);
	return *this;
};
//============================================================


