#ifndef GUARD_config_h
#define GUARD_config_h
////////// MoDAP //////////
// config.h
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <list>
#include <map>
#include <iostream>

//================================
// config class (graph-like)
//--------------------------------

class config
{
private:
	// variables
	int nAtom, nType, nBond, nMole;
	std::map<std::string, int> nAtomType;
	vector<atom> atoms;
	vector<molecule> moles;
/*
	std::vector<node> nodes;
    // classes
    class node
    {
    private:
		atom at;                  // node's atom
		std::vector<bond> bonds;  // atom's bonds (weighted graph edges)
		
		node(atom a)
		{
			this->at = a;
		}
	};
*/
	// methods
	void merge_molecules();
	
public:
    // classes
    class iterator
	{
        private:
            pointer ptr_;
        public:
//            typedef iterator self_type;
            typedef atom value_type;
            typedef atom& reference;
            typedef atom* pointer;
            typedef std::forward_iterator_tag iterator_category;
            typedef int difference_type;
            iterator(atom* ptr) : ptr_(ptr) { }
            iterator operator++() { iterator i = *this; ptr_++; return i; }
            iterator operator++(int junk) { ptr_++; return *this; }
            reference operator*() { return *ptr_; }
            pointer operator->() { return ptr_; }
            bool operator==(const iterator& rhs) { return ptr_ == rhs.ptr_; }
            bool operator!=(const iterator& rhs) { return ptr_ != rhs.ptr_; }
    };

    // variables
	int tstep;
    lattice lat;
	
	// constructors
	config();
	
	// methods
	iterator begin();
	iterator end();
	
	int natoms() const;
	int natoms(const std::string&) const;
	int ntypes() const;
	int nbonds() const;
	int nmoles() const;
	
	void insert(atom);
	void insert(molecule);
	
	int trjRead(std::ifstream& trj);
	void trjWrite(std::ofstream& trj);
//	void xyzRead(std::ifstream& xyz);
	void xyzWrite(std::ofstream& xyz);
//	int poscarRead(std::ifstream& poscar);
//	int xdatcarRead(std::ifstream& xdatcar, const config& cfg);

	void wrap();
	config supercell(const int&, const int&, const int&) const;
//	std::list<double> distances(const std::string&);
	std::list<double> distances(const std::string&, const std::string&);
	void molecules(std::map<std::string,double>& cutoffs);
	
	// operators
	config& operator+=(const atom& a);
//	config& operator+=(const bond& b);
	config& operator+=(const molecule& m);
};
//================================

#endif