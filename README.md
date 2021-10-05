# mscheme
mscheme.cxx is ROOT macro, which calculates M and J values from all possible combination of spins
for N fermions in the same orbital
code by Matus Sedlak
****************************************************************************************

To run macro:
1. open ROOT in mscheme directory
2. type .L mscheme.cxx
3. type mscheme(spin,n)

spin is 2* spin of particles, e.g. when spin is 7/2, put 7 instead
n is number of particles, e.g. 2 

output is written in terminal window and histogram of M values is shown
