#ifndef GUARD_alphabetical_h
#define GUARD_alphabetical_h
////////// MoDAP //////////
// alphabetical.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>

// class for alphabetical sorting strings by first character only
// needed for map sorting in molecule class
class alphabetical {
    public:
	    bool operator() (const std::string& a, const std::string& b){
        	return tolower(a[0]) < tolower(b[0]);
        }
};

#endif


