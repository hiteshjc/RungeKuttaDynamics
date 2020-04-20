#ifndef RUNGE_HEADER
#define RUNGE_HEADER

#include"global.h"
#include"math_utilities.h"

void time_evolve(   double spin, double deltat, 
		    double omega, double tottime, int L,
		    std::vector<double> &configx, 
		    std::vector<double> &configy, 
		    std::vector<double> &configz, 
		    std::vector< std::vector<int> > &neighbors, 
		    std::vector< std::vector<int> > &nneighbors, 
		    RMatrix &Jmat01, RMatrix &Jmat10,
		    RMatrix &Jmat02, RMatrix &Jmat20,
		    RMatrix &Jmat03, RMatrix &Jmat30,
		    RMatrix &Jmat12, RMatrix &Jmat21,
		    RMatrix &Jmat13, RMatrix &Jmat31,
		    RMatrix &Jmat23, RMatrix &Jmat32,
		    double &Jnnn, RMatrix &bond_disorder_matrix,
		    std::vector< std::vector<double> > &fullcoords,
		    std::vector< std::vector<int> > & ijkt, STensor &smunu);
		   
#endif
