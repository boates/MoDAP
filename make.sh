#!/bin/bash

echo ""
echo "=== beginning make ==="

echo ""
echo "removing *.o"
rm -f *.o

echo ""
echo "=== g++ -c *.cpp ==="

echo "utils.cpp"
g++ -c utils.cpp
echo "atom.cpp"
g++ -c atom.cpp
#echo "bond.cpp"
#g++ -c bond.cpp
#echo "molecule.cpp"
#g++ -c molecule.cpp
echo "lattice.cpp"
g++ -c lattice.cpp
#echo "config.cpp"
#g++ -c config.cpp
#echo "physics.cpp"
#g++ -c physics.cpp

#echo "molecules.cpp"
#g++ -c molecules.cpp
#echo "pcf.cpp"
#g++ -c pcf.cpp
#echo "coordination.cpp"
#g++ -c coordination.cpp
#echo "trjvasp.cpp"
#g++ -c trjvasp.cpp
echo "=== finished all cpp ==="

echo ""
#echo "linking *.o for lattice"
#g++ -o lattice lattice.o utils.o
#echo "linking *.o for atom"
#g++ -o atom atom.o lattice.o utils.o

#echo "linking *.o for molecules"
#g++ -o molecules molecules.o atom.o bond.o molecule.o lattice.o config.o physics.o utils.o
#echo "linking *.o for pcf"
#g++ -o pcf pcf.o atom.o bond.o molecule.o lattice.o config.o physics.o utils.o
#echo "linking *.o for coordination"
#g++ -o coordination coordination.o atom.o bond.o molecule.o lattice.o config.o physics.o utils.o
#echo "linking *.o for trjvasp"
#g++ -o trjvasp trjvasp.o atom.o bond.o molecule.o lattice.o config.o physics.o utils.o


echo ""
echo "=== make successful ==="
echo ""