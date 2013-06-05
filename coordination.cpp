////////// MoDAP //////////
// coordination.cpp
// Author: Brian Boates
///////////////////////////
#include <fstream>
#include <cmath>
#include <algorithm>
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
void read_input(string& ftrj, vector< pair<string,string> >& pairs) {
    
	// open COORDIN file (naming fixed)
 	string fin = "COORDIN";
	ifstream input(fin.c_str());
	
	// check for successful file opening
	if (input.is_open()){

	    // read in trj filename
		input >> ftrj;

        // read rest of input
        while ( !input.eof() ) {

	        // read atom types
			string t1 = ""; string t2 = "";
			input >> t1 >> t2;
			
			// if end-of-file, break to avoid reading empty atom pair
			if (input.eof()) {break;}

   			// create pair object for atom types
       		pair<string,string> p;
       		p.first = t1; p.second = t2;
			
            // add pair to pairs vector
    		pairs.push_back(p);
        }
		input.close();
	}
	else
	    cout << "Unable to open file: " << fin << endl;
}

int main() {

	// read input file
	string ftrj; vector< pair<string,string> > pairs;
	read_input(ftrj,pairs);
	
	// make map to hold a list of coordinations for each pair
	map< string, list<double> > dmap;
	
	// open trj file
	ifstream trj(ftrj.c_str());
	
	// variables for later pcf normalization
	int tsteps = 0;
	int natom = 0;
	double vol = 0;
	double pi = acos(-1);

	// loop through xyz file
	while(!trj.eof()){

		// initilizle config (a snapshot of atoms)
		config cfg;

        // read in a single trj configuration
		cfg.read_trj(trj);
		
		// check for end of file and break if reached
		if (trj.eof()) {break;}

		// set natom and vol
		natom = cfg.get_natom();
		vol = cfg.lat.volume();
		
		// increment number of steps
		tsteps += 1;
		
		// define variables
		string pname; list<double> d;
		
		// loop over desired atom pairs (different pcf's)
		typedef vector< pair<string,string> >::const_iterator iter;
		for (iter p=pairs.begin();p!=pairs.end();++p){

			// shorthand for pair
			pname = p->first + p->second;
			
			// get list of distances for given pair from current cfg
			d = cfg.distances(p->first, p->second);

            // append current cfg distances to end of overall list for pair
			dmap[pname].insert(dmap[pname].end(), d.begin(), d.end());
		}
	}
    // close trajectory file
	trj.close();

    //====================================================//
	// at this point dmap has complete distance lists for //
	// all desired pairs over the entire trajectory       //
    //====================================================//
	
	// default (hard-coded) bin_size = 0.02
	// should never really be a reason to change this for a pcf
	double bin_size = 0.05;
	
	// loop over desired atom pairs (different pcf's)
	typedef vector< pair<string,string> >::const_iterator iter;
	for (iter p=pairs.begin();p!=pairs.end();++p){
		
    	// compute a histogram for each pair's distance list
		string pname = p->first+p->second;
		
		// create pcf filename for pair and open output file
		string fname = "PCF_"+pname+".dat";
		ofstream out(fname.c_str());
		
		// number of bins is rmax/bin_size
		int nbins = *max_element(dmap[pname].begin(), dmap[pname].end()) / bin_size;
		
		// retrieve histogram for current atom pair's distances
		pair< vector<double>,vector<double> > xyhist;
		xyhist = histogram(dmap[pname], nbins);
		
		// assign r and hist as first and second from xyhist pair
		vector<double> r = xyhist.first;
		vector<double> hist = xyhist.second;
		
		// loop through the histogram and write pcf to file
		for (vector<double>::size_type i=0;i!=r.size();++i){
			
			// normalize histogram value for pcf
			double norm = (natom-1)*(natom/vol)*(4*pi*pow(r[i],2)*3bin_size) *tsteps;
			double pcf = hist[i] / norm;
			
			// write r and pcf to file
			out << r[i] << " " << pcf << endl;
		}
		out.close();
	}
	return 0;
}



