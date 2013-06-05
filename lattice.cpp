////////// MoDAP //////////
// lattice.cpp
// Author: Brian Boates
///////////////////////////
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "lattice.h"
#include "utils.h"
using namespace std;

//============================================================
// members of lattice class (constructors, methods, etc.)
//------------------------------------------------------------

//== constructors ==//

lattice::lattice() : ax(0), ay(0), az(0), bx(0), by(0), bz(0), cx(0), cy(0), cz(0) {}

lattice::lattice(double l) : ax(l), ay(0), az(0), bx(0), by(l), bz(0), cx(0), cy(0), cz(l) {}

lattice::lattice(double l, double m, double n) : ax(l), ay(0), az(0),
                                                 bx(0), by(m), bz(0),
                                                 cx(0), cy(0), cz(n) {}

lattice::lattice(double lx, double ly, double lz,
                 double mx, double my, double mz,
                 double nx, double ny, double nz) :
                 ax(lx), ay(ly), az(lz),
                 bx(mx), by(my), bz(mz),
                 cx(nx), cy(ny), cz(nz) {}

//== methods ==//

// lattice vector a
vector<double> lattice::a() const
{
	double array[3] = {this->ax, this->ay, this->az};
	vector<double> avec(array, array+3);
	return avec;
}

// lattice vector b
vector<double> lattice::b() const
{
	double array[3] = {this->bx, this->by, this->bz};
	vector<double> bvec(array, array+3);
	return bvec;
}

// lattice vector c
vector<double> lattice::c() const
{
	double array[3] = {this->cx, this->cy, this->cz};
	vector<double> cvec(array, array+3);
	return cvec;
}

// magnitude of lattice vector a
double lattice::amag() const
{
	return magnitude(this->a());
}

// magnitude of lattice vector b
double lattice::bmag() const
{
	return magnitude(this->b());
}

// magnitude of lattice vector c
double lattice::cmag() const
{
	return magnitude(this->c());
}

// lattice vector angles alpha = cb, beta = ac, gamma = ab
vector<double> lattice::angles() const
{
    // build angle vector to return
	vector<double> angs;
	angs.push_back( angle(this->b(),this->c()) );  // alpha
	angs.push_back( angle(this->a(),this->c()) );  // beta
	angs.push_back( angle(this->a(),this->b()) );  // gamma
	return angs;
}

// lattice volume
double lattice::volume() const
{
	return abs( dot(this->a(), cross(this->b(),this->c())) );
}

// print lattice info
void lattice::print() const
{
	// print lattice
	cout << "#" << endl << "# ====== printing lattice ====== #" << endl;
    cout << "#" << endl;

	// print magnitudes
	cout << "# amag, bmag, cmag = ";
	cout << this->amag() << ", " << this->bmag() << ", " << this->cmag() << endl;

    // print angles
	vector<double> angs = this->angles();
	cout << "# alpha, beta, gamma = ";
	cout << angs[0] << ", " << angs[1] << ", " << angs[2] << endl;

    // print volume
	cout << "# volume = " << this->volume() << endl;
    cout << "#" << endl;

    // lattice vector components
    cout << "# ax, ay, az = ";
	cout << ax << ", " << ay << ", " << az << endl;
    cout << "# bx, by, bz = ";
	cout << bx << ", " << by << ", " << bz << endl;
    cout << "# cx, cy, cz = ";
	cout << cx << ", " << cy << ", " << cz << endl;
    cout << "#" << endl;

	cout << "# ====== end of lattice ======== #" << endl << "#" << endl;
}
//============================================================
/*
int main()
{
    // create artificial lattice and print it
	lattice lat(7.1, 8.2, 9.3);
	lat.print();

	return 0;
}
*/

