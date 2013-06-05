////////// MoDAP //////////
// pcf.cpp
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
void read_input(string& ftrj, double& rmax, vector< pair<string,string> >& pairs) {
    
	// open PCFIN file (naming fixed)
 	string fin = "PCFIN";
	ifstream input(fin.c_str());
	
	// check for successful file opening
	if (input.is_open()){

	    // read in trj filename
		input >> ftrj;
		
		// read in max distance to compute pcf
		input >> rmax;

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
	string ftrj; double rmax; vector< pair<string,string> > pairs;
	read_input(ftrj,rmax,pairs);
	
	// make map to hold a list of distances for each pair
	map< string, list<double> > dmap;
	
	// make map to hold natom for each atom type
	map< string, int > natoms;
	
	// open trj file
	ifstream trj(ftrj.c_str());
	
	// variables for later pcf normalization
	int tsteps = 0;
	int natom = 0;
	double vol = 0; double a = 0; double b = 0; double c = 0;
	double pi = acos(-1);

	// loop through trj file
	while(!trj.eof()){

		// initilizle config (a snapshot of atoms)
		config cfg;

        // read in a single trj configuration
		cfg.read_trj(trj);
		
		// check for end of file and break if reached
		if (trj.eof()) {break;}

		// set natom, vol, and lattice vector magnitudes
		natom = cfg.get_natom();
		vol = cfg.lat.volume();
		a = cfg.lat.amag;
		b = cfg.lat.bmag;
		c = cfg.lat.cmag;
		
		// get natom for each type in the pair
		typedef vector<string>::const_iterator typit;
		for (typit i=cfg.types.begin();i!=cfg.types.end();++i) {
			natoms[*i] = cfg.get_natom(*i);
		}
		
		// increment number of steps and set for current cfg
		tsteps += 1;
		cfg.tstep = tsteps;
		
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
	double bin_size = 0.02;

    // for proper normalization
	int natom1 = 0;
	int natom2 = 0;
	
	// loop over desired atom pairs (different pcf's)
	typedef vector< pair<string,string> >::const_iterator iter;
	for (iter p=pairs.begin();p!=pairs.end();++p){
		
    	// compute a histogram for each pair's distance list
		string pname = p->first+p->second;
		
		// create pcf filename for pair and open output file
		string fname = "PCF_"+pname+".dat";
        ofstream out(fname.c_str());
		
		// get natom1 and natom2 for current pair
		natom1 = natoms[p->first];
		natom2 = natoms[p->second];
		
		// choose max and min values for the PCF histogram
		// double mx = *max_element(dmap[pname].begin(), dmap[pname].end());
		double mx = rmax;
		double mn = 0;
		
		// retrieve histogram for current atom pair's distances
		pair< vector<double>, vector<double> > xyhist;
		xyhist = histogram(dmap[pname], bin_size, mx, mn);

		// assign r and hist as first and second from xyhist pair
		vector<double> r = xyhist.first;
		vector<double> hist = xyhist.second;
		
		// define normalization factor based on double counting and like-atoms
		double factor = natom1*(natom2/vol);
		if (p->first == p->second) {
			factor = (natom1-1)*(natom2/vol) / 2.0;
		}
		
		// loop through the histogram and write pcf to file
		for (vector<double>::size_type i=0;i!=r.size();++i){
			
			// normalize histogram value for pcf //
			
			// normalization exceptions for r > a/2, b/2, c/2 (orthrhombic) (from schwegler's RDF)
			double cx1 = 0.0; double cx2 = 0.0;
			double cy1 = 0.0; double cy2 = 0.0;
			double cz1 = 0.0; double cz2 = 0.0;
			double rdr = r[i] + bin_size;
	        if (r[i] > a/2.0) { cx1 = pi/12.0 * (pow(a,3) - 12.0*a*pow(r[i],2) + 16.0*pow(r[i],3)); }
	        if (r[i]+bin_size > a/2.0) { cx2 = pi/12.0 * (pow(a,3) - 12.0*a*pow(rdr,2) + 16.0*pow(rdr,3)); }
	        if (r[i] > b/2.0) {	cy1 = pi/12.0 * (pow(b,3) - 12.0*b*pow(r[i],2) + 16.0*pow(r[i],3)); }
	        if (r[i]+bin_size > b/2.0) { cy2 = pi/12.0 * (pow(b,3) - 12.0*b*pow(rdr,2) + 16.0*pow(rdr,3)); }
	        if (r[i] > c/2.0) {	cz1 = pi/12.0 * (pow(c,3) - 12.0*c*pow(r[i],2) + 16.0*pow(r[i],3)); }
	        if (r[i]+bin_size > c/2.0) { cz2 = pi/12.0 * (pow(c,3) - 12.0*c*pow(rdr,2) + 16.0*pow(rdr,3)); }

            // put it all together for complete normalization		
			// double shell = 4.0*pi/3.0 * ( pow(rdr,3) - pow(r[i],3) ) - dc;
			double dc = (cx2-cx1) + (cy2-cy1) + (cz2-cz1);
			double shell = 4.0*pi*pow(r[i],2)*bin_size - dc;
			double norm = factor * shell * tsteps;
			
			// finally compute normalized pcf value
			double pcf = hist[i] / norm;
			
			// write r and pcf to file
			out << r[i] << " " << pcf << endl;
		}
		// close the pcf output file
		out.close();
	}
	return 0;
}


















