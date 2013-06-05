////////// MoDAP //////////
// trjvasp.cpp
// Author: Brian Boates
///////////////////////////
#include <fstream>
#include <cmath>
#include <algorithm>

#include "atom.h"
#include "bond.h"
#include "molecule.h"
#include "lattice.h"
#include "config.h"
#include "utils.h"
using namespace std;

// read in POSCAR and XDATCAR and create a .trj file (XDATCAR.trj)
int main() {

	// open PCFIN file (naming fixed)
 	string fposcar = "POSCAR";   ifstream poscar(fposcar.c_str());
	string fxdatcar = "XDATCAR"; ifstream xdatcar(fxdatcar.c_str());
	
	// check for successful file opening
	if (poscar.is_open() && xdatcar.is_open()){

   		// read in POSCAR file
		config cfg_poscar;
   		int version = cfg_poscar.read_poscar(poscar);
			
		// strip the XDATCAR header info (5 lines)
		string header;
		for (int i = 0; i!=5; ++i) { getline(xdatcar,header); }
		
		// timestep counter
		int tsteps = 0;
		
		// open .trj output file
		string ftrj = "XDATCAR.trj"; ofstream trj(ftrj.c_str());
		
    	// loop through trj file
	    while(!xdatcar.eof()){

    		// initilizle config (a snapshot of atoms)
	    	config cfg;
		
            // read in a single trj configuration
	    	cfg.read_xdatcar(xdatcar, cfg_poscar);
		
    		// check for end of file and break if reached
	    	if (xdatcar.eof()) {break;}

    		// increment number of steps and set for current cfg
	    	tsteps += 1;
			cfg.tstep = tsteps;

	        // write to .trj file
			cfg.write_trj(trj);
       	}

        // close POSCAR, XDATCAR, and .trj files
    	poscar.close();
	    xdatcar.close();
    	trj.close();
		
	} else
		cout << "Could not find POSCAR and/or XDATCAR --- exiting..." << endl;
		
	return 0;
}
