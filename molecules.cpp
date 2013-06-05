////////// MoDAP //////////
// molecules.cpp
// Author: Brian Boates
///////////////////////////
#include <fstream>
#include <cmath>
// #include <iomanip>

#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include "lattice.h"
#include "config.h"
#include "alphabetical.h"
#include "physics.h"
#include "utils.h"
using namespace std;

// function to read and parse the input file
void read_input(string& ftrj, map<string,double>& cutoffs) {
    
	// open MOLIN file (naming fixed)
//  cout << "Name of input file:" << endl;
// 	string fin; cin >> fin;
 	string fin = "MOLIN";
	ifstream input(fin.c_str());
	
	// check for successful file opening
	if (input.is_open()){

	    // read in trj filename
		input >> ftrj;

        // read rest of input
        while ( !input.eof() ) {

	        // read atom pair and cutoff
			string p = ""; double rcut = 0;
			input >> p >> rcut;
			
            // add to cutoff map
			cutoffs[p] = rcut;			
        }
		input.close();
	}
	else
	    cout << "Unable to open file: " << fin << endl;
}

int main() {

	// read input file
	string ftrj; map<string,double> cutoffs;
	read_input(ftrj,cutoffs);
	
	// open trj file
	ifstream trj(ftrj.c_str());
	
	// create map to hold species counts
	map<string,double> mol_counter;
	
	// create timestep counter for later averaging
	int tsteps = 0;
	
	// loop through xyz file
	while(!trj.eof()){

		// initilizle config (a snapshot of atoms)
		config cfg;

        // read in a single trj configuration
		cfg.read_trj(trj);
		
		// check for end of file and break if reached
		if (trj.eof()) {break;}

		// increment number of steps
		tsteps += 1;

        // find list of molecules in cfg (into cfg.mols)
		cfg.molecules(cutoffs);
		
		// loop over molecules in cfg
		typedef list<molecule>::iterator iter;
		for (iter i=cfg.mols.begin();i!=cfg.mols.end();++i){
			
			// retrieve current molecular species name
			string spec = i->species();
			
			// check to see if already in map
			// (for careful treatment of initialized count value)
			if (mol_counter.count(spec) > 0){
    			// add a count for each molecule detected
				mol_counter[spec] += 1.0;
			} else
			    // otherwise it is the first appearance of this species
				mol_counter[spec] = 1.0;
		}		
	}

    // close trajectory file
	trj.close();

	// mol_counter has now counted all species in the system //
    
    // create inverted key/value ---> value/key map for sorting
	map<int,string> sorted;

    // define iterator for the map
	typedef map<string,double>::iterator m_iter;
	
	// for avoiding duplicate counts in map
	double x = 0.01;

    // loop through mol_counter to create sorted map
	for (m_iter i=mol_counter.begin();i!=mol_counter.end();++i){

		// in case two species have identical counts
		if (sorted.count(i->second) > 0) {
			
			// slightly change the count for the map
			i->second += x;
			
			// increment small x value to make unique again
			x += 0.01;
		}
		// should auto-sort by < for keys
		sorted[ i->second ] = i->first;
	}

	// open output file and write header
	ofstream out("MOLOUT");
	out << "# species, count/tsteps, total count" << endl;
	
    // loop through sorted map and write to file
	typedef map<int,string>::iterator s_iter;
	for (s_iter i=sorted.begin();i!=sorted.end();++i){
		// write molecule species, count per tstep, and total count
		out << i->second << " " << (i->first)/tsteps << " " << i->first << endl;
	}
	
	// close output file
	out.close();
	
	return 0;
}


















