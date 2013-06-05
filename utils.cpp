////////// MoDAP //////////
// utils.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <cctype>
#include "utils.h"
using namespace std;

// factorial function
int factorial(const int& n){
	int f=1;
	for (int i=1;i!=n+1;++i) {
		f *= i;
	}
	return f;
}

// m choose n
int choose(const int& m, const int& n){
	return factorial(m) / (factorial(n)*factorial(m-n));
}

// sign function
int sgn(const double x){
	// return +1 for x >= 0
	if (x >= 0) {return 1;}
	else if (x < 0) {return -1;}
}

// multiply n-dimensional vector by a double x: i.e. x*vector
vector<double> doubleXvec(const vector<double> vec, const double x){

	    // make sure vec is not empty
		if (vec.size()==0){
			throw length_error("attempted multiplication with an empty vector");
		}

		// otherwise, compute x*vec
		else {
			// initialize new vector as copy of vec
		    vector<double> v (vec);
		    // loop through vec and compute x*vec
	    	for (vector<double>::size_type i=0;i!=vec.size();++i){
			    v[i] = x*vec[i];
		    }
	    	return v;
	    }
}

// n-dimensional vector addition: a + b
vector<double> vec_add(const vector<double> a, const vector<double> b){

		// make sure vectors a and b are the same size
		if (a.size()!=b.size()){
			throw length_error("addition of differently sized vectors");
		}

	    // make sure vectors a and b are not empty
		else if (a.size()==0){
			throw length_error("addition of empty vectors");
		}

		// compute a + b
		else {
		    vector<double> v;
	    	for (vector<double>::size_type i=0;i!=a.size();++i){
				v.push_back( a[i] + b[i] );
		    }
	    	return v;
	    }
}

// n-dimensional vector subtraction: a - b
vector<double> vec_subtract(const vector<double> a, const vector<double> b){
	
		// make sure vectors a and b are the same size
		if (a.size()!=b.size()){
			throw length_error("subtraction of differently sized vectors");
		}

	    // make sure vectors a and b are not empty
		else if (a.size()==0){
			throw length_error("subtraction of empty vectors");
		}

		// compute a - b
		else {
		    vector<double> v;
	    	for (vector<double>::size_type i=0;i!=a.size();++i){
			    v.push_back( a[i] - b[i] );
		    }
	    	return v;
	    }
}

// dot product of two n-dimensional vectors
double dot(const vector<double>& a, const vector<double>& b){

	// make sure vectors a and b are the same size
	if (a.size()!=b.size()){
		throw length_error("dot product of differently sized vectors");
	}

    // make sure vectors a and b are not empty
	else if (a.size()==0){
		throw length_error("dot product of empty vectors");
	}
	
	// compute the dot product of vectors a and b
	else {
	    double v = 0;
    	for (vector<double>::size_type i=0;i!=a.size();++i){
		    v += a[i]*b[i];
	    }
    	return v;
    }
}

// cross product of two 3-dimensional vectors
vector<double> cross(const vector<double>& a, const vector<double>& b){

	// make sure vectors a and b are the same size = 3
	if (a.size()!=3 || b.size()!=3){
		throw length_error("cross product vectors not of length 3");
	}

	// compute the cross product of vectors a and b
	else {
		vector<double> vec;
		vec.push_back( a[1]*b[2] - a[2]*b[1] );
		vec.push_back( a[2]*b[0] - a[0]*b[2] );
		vec.push_back( a[0]*b[1] - a[1]*b[0] );
    	return vec;
    }
}

// absolute magnitude of n-dimensional vector
double magnitude(const vector<double>& v){

	// make sure vector v is not empty
	if (v.size()==0){
		throw length_error("magnitude of an empty vector");
	}

    // compute the absolute magnitude of vector v	
	else{
        double mag = 0;
    	for (vector<double>::size_type i=0;i!=v.size();++i){
	    	mag += pow(v[i],2);
        }
	    return pow(mag,0.5);
    }
}

// angle of two 3-dimensional vectors
double angle(const vector<double>& a, const vector<double>& b){

	// make sure vectors a and b are the same size = 3
	if (a.size()!=3 || b.size()!=3){
		throw length_error("angle of vectors not of length 3");
	}
	
    // compute a-b angle and return in degrees
	double theta = acos( dot(a,b) / (magnitude(a)*magnitude(b)) );
	return (180.0 / acos(-1)) * theta;	
}

// true if arg is whitespace, false otherwise
bool space(const char& c) {
	return isspace(c);
}

// false if the arg is whitespace, true otherwise
bool not_space(const char& c) {
	return !isspace(c);
}

// split a string by spaces
vector<string> split(const string& str) {

    // define iterator and return vector
	typedef string::const_iterator iter;
	vector<string> ret;
	
	iter i = str.begin();
	while (i != str.end()) {

		// ignore leading blanks
		i = find_if(i, str.end(), not_space);

		// find the end of the next word
		iter j = find_if(i, str.end(), space);
		
		// copy the characters in [i,j)
		if (i != str.end()) {
			ret.push_back(string(i,j));
		}
		i = j;
	}
	return ret;
}












