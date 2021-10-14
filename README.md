# mscheme
mscheme.cxx is ROOT macro, which calculates M and J values from all possible combination of spins for N fermions in the same orbital

written by Matus Sedlak
****************************************************************************************

To run macro:
1. open ROOT in mscheme directory
2. type .L mscheme.cxx
3. type mscheme(spin,n)

Input spin that is a double of an orbital spin, e.g. when spin is 7/2, put 7 instead

n is number of particles/holes, e.g. 2 

Output is written in terminal window and histogram of M values is shown
Please keep in mind, that resulting m-scheme will contain all seniorities. It is faster to calculate J configurations by specifying number of holes instead of particles after mid-shell.
