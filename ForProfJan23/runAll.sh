#!/bin/bash
DV=1

g++-7 -std=c++1z -O3 runBMv61_PIR_Power2.cpp -o runBMv61_PIR_Power2 -mavx2 -march=native -lbsd
g++-7 -std=c++1z -O3 runBMv62_Half.cpp -o runBMv62_Half -mavx2 -march=native -lbsd
g++-7 -std=c++1z -O3 runBMv63_Soundness.cpp -o runBMv63_Soundness -mavx2 -march=native -lbsd
g++-7 -std=c++1z -O3 runBMv64_PIR_arbitrary.cpp -o runBMv64_PIR_arbitrary -mavx2 -march=native -lbsd


./runBMv63_Soundness
./runBMv64_PIR_arbitrary
./runBMv61_PIR_Power2
./runBMv62_Half
