#ifndef GUARD_utils_h
#define GUARD_utils_h
////////// MoDAP //////////
// utils.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <sstream>

int factorial(const int&);
int choose(const int&, const int&);
int sgn(const double);

std::vector<double> vec_add(const std::vector<double>, const std::vector<double>);
std::vector<double> vec_subtract(const std::vector<double>, const std::vector<double>);
std::vector<double> doubleXvec(const std::vector<double>, const double);
double dot(const std::vector<double>&, const std::vector<double>&);
std::vector<double> cross(const std::vector<double>&, const std::vector<double>&);
double magnitude(const std::vector<double>&);
double angle(const std::vector<double>&, const std::vector<double>&);
bool space(const char&);
bool not_space(const char&);
std::vector<std::string> split(const std::string&);
std::vector<int> vstr_to_vint(const std::vector<std::string>&);

//============================================================
// template class function must be entirely in header file
//------------------------------------------------------------

// if two vectors<T> have a common value
template <class T>
bool common_value(const std::vector<T>& a,const std::vector<T>& b){
	
	// define iterator type
	typedef typename std::vector<T>::const_iterator iter;

	// loop through vector<T> a
	for (iter it=a.begin(); it!=a.end(); ++it){
		
		// see if current a value is in b
		if ( find(b.begin(), b.end(), *it) != b.end() ){
			
			// if it is, return true (a & b share a common value)
			return true;
		}
	}	
	// otherwise return false
	return false;
}

// remove duplicate values from a vector
template <class T>
void remove_duplicates(std::vector<T>& v){
	
	// create empty vector<T>
	std::vector<T> vec;

	// define iterator type
	typedef typename std::vector<T>::const_iterator iter;
	
	// loop through given vector
	for (iter it=v.begin();it!=v.end();++it){
		
		// if current item not yet in vec
		if ( find(vec.begin(),vec.end(),*it) == vec.end() ){
			
			// append new item
			vec.push_back(*it);
		}
	}
	// replace original vector with new duplicate removed one
	v = vec;
}

// shorthand command to print a vector (like in python)
template <class T>
void printvec(const std::vector<T>& v){

    // define iterator type
	typedef typename std::vector<T>::const_iterator iter;

	// loop over vector elements and print values
	std::cout << "[ ";
	for (iter i=v.begin();i!=v.end();++i){
		std::cout << *i;
		if (i!=v.end()-1) {    // comma after every value
			std::cout << ", "; // except the last
		}
	}
	std::cout << " ]" << std::endl;
}

// shorthand command to print a list (like in python)
template <class T>
void printlist(const std::list<T>& l){

    // define iterator type
	typedef typename std::list<T>::const_iterator iter;

	// loop over vector elements and print values
	std::cout << "[ ";
	for (iter i=l.begin();i!=l.end();++i){
		std::cout << *i;
		std::cout << ", ";
	}
	std::cout << " ]" << std::endl;
}

// make a histogram: T should be either vector<#> or list<#>  // const int& nbins
template <class T>
std::pair< std::vector<double>,std::vector<double> > histogram(const T& c, const double& bin_size, double mx = -999, double mn = -999, const double norm = 1.0){

    // define size_type to use throughout
	typedef std::vector<double>::size_type vint;

	// If mx and mn not specified, set defaults as max and min values found in c
	if (mx == -999) { mx = *max_element(c.begin(), c.end()); } // max
	if (mn == -999) { mn = *min_element(c.begin(), c.end()); } // min

    // total number of bins	
	vint nbins = (mx - mn) / bin_size;

	// define the "minimum bin" (bin for smallest value)
	vint mbin = mn / bin_size;
	
	// define allocation parameter
	vint alloc = mbin + nbins;
	
    // return the histogram as a vector of doubles (allocate for nbins values = 0)
	std::vector<double> hist(alloc);

	// define "x-axis" for histogram as a vector of doubles (allocate for nbins values)
	std::vector<double> x(alloc);
	
	// set "x-axis" values and initialize hist values to 0
	for (std::vector<double>::size_type i=mbin; i!=alloc; ++i){

		// "x-axis" values
		x[i] = i*bin_size + bin_size/2.0;

        // initialize histogram values to zero
		hist[i] = 0;
	}
	
	// define iterator for c
//	typedef typename T::const_iterator iter;
	typedef typename T::const_iterator iter;
//	typedef std::list<double>::const_iterator iter;
	
	// loop through c and put values in their respective bins
	std::vector<double>::size_type bin;

	printvec(x);

	for (iter j=c.begin(); j!=c.end(); ++j){
		
		printvec(x);

		// determine bin for current value
		bin = *j / bin_size;
		
		// add 1 to respective bin (normalized if desired, default: norm=1)
		hist[bin] += 1.0 / norm;
	}
	
	printvec(x);
	
	// define pair to pass back "x" and "y" axis for histogram
	std::pair< std::vector<double>, std::vector<double> > xy;
	xy.first = x;
	xy.second = hist;

	return xy;
}

// convert a string to a number (int/double)
template <typename T>
void str_to_num(const std::string& str, T& num) {
	std::istringstream convert(str); // stream for int to string conversion
	convert >> num;
}

// convert a number (int/double) to a string
template <typename T>
void num_to_str(const T& num, std::string& str) {
	std::ostringstream convert; // stream for int to string conversion
	convert << num;
	str += convert.str();
}

// convert a vector<string> to a vector<T>
template <typename T>
void vstr_to_vnum(const std::vector<std::string>& vstr, std::vector<T>& vnum) {
	
	// loop over values in vstr
	typedef std::vector<std::string>::const_iterator iter;
	for (iter i=vstr.begin(); i!=vstr.end(); ++i) {

		// convert string to int
		T n; str_to_num(*i,n);
		
		// append int to vint
		vnum.push_back( n );
	}	
}

//============================================================

#endif