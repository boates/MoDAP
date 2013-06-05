////////// MoDAP //////////
// bond.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <iostream>
//#include "atom.h"
#include "bond.h"
using namespace std;

//============================================================
// members of bond class (constructors, methods, etc.)
// this is basically like a weighted edge for a graph object
//------------------------------------------------------------

//== constructors ==//

bond::bond(int a1, int a2, double len) : at1(a1), at2(a2), length(len) {}

//== methods ==//

// print atom indices and bond length
void bond::print() const
{
	cout << "{ " at1 << ", " << at2 << ", length=" << length << " }" << endl;
//	cout << "{ ";
//	at1.print();
//	cout << ", ";
//	at2.print();
//	cout << ", length=" << length << " }" << endl;
}

// alternate print function
void bond::print(string str) const
{
    // print atom indices
    if (str == "atoms") {
		cout << "{ " << at1 << ", " << at2 << " }" << endl;
//		cout << "{ ";
//		at1.print();
//		cout << ", ";
//		at2.print();
//		cout << " }" << endl;
    }
    // print length
    else if (str == "length") cout << "length=" << length << endl;
    // print options
    else cout << "options: \"atoms\", \"length\" (all: give no arg)" << endl;
}

//== operators ==//

// compare atom pair and length (also check reversed atom pair)
bool bond::operator==(const bond& b) const
{
	if (at1 == b.at1 && at2 == b.at2 && length == b.length) return true;
	else return (at1 == b.at2 && at2 == b.at1 && length == b.length);
}
bool bond::operator> (const bond& b) const
{
	return length > b.length;
}
bool bond::operator>=(const bond& b) const
{
	return length >= b.length;
}
bool bond::operator< (const bond& b) const
{
	return length < b.length;
}
bool bond::operator<=(const bond& b) const
{
	return length <= b.length;
}

//============================================================
/*
int main()
{
	// create a bond and print it
	atom a1("N", 1, 1.2, 3.4, 5.6);
	atom a2("C", 2, 2.1, 4.3, 6.5);
	bond b(a1, a2, 1.456);
	b.print();
	b.print("atoms");
	b.print("length");
	cout << b.length << endl;
    return 0;
}
*/






