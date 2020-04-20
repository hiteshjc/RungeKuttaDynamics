#ifndef MC_HEADER
#define MC_HEADER

#include"global.h"
#include"math_utilities.h"

void mc(string lattice, int nstates, int L, double temp, double hx, double hy, double hz, double & eavg, 
	double &mxavg, double &myavg, double &mzavg, double &e2avg, 
	double &m2xavg, double &m2yavg, double &mz2avg);

#endif
