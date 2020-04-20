#ifndef MC_PYRO_HEADER
#define MC_PYRO_HEADER

#include"global.h"
#include"math_utilities.h"

double total_j_energy(double &spin, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors, 
		      std::vector< std::vector<int> > &nneighbors, 
		      RMatrix &Jmat01, RMatrix &Jmat10,
		      RMatrix &Jmat02, RMatrix &Jmat20,
		      RMatrix &Jmat03, RMatrix &Jmat30,
		      RMatrix &Jmat12, RMatrix &Jmat21,
		      RMatrix &Jmat13, RMatrix &Jmat31,
		      RMatrix &Jmat23, RMatrix &Jmat32,
		      double &Jnnn, RMatrix &bond_disorder_matrix,
		      std::vector< std::vector<int> > & ijkt);
void make_qs(int L, std::vector<std::vector<double> > &qvals);

void make_phases( std::vector< std::vector<double> > &fullcoords,
		  std::vector< std::vector<double> > &qvals,
		  Matrix &phases);

void fourier_transforms(double spin,Matrix &phases, 
			std::vector<double> &configx,
			std::vector<double> &configy,
			std::vector<double> &configz,
			std::vector<complex<double> > &sxq,
			std::vector<complex<double> > &syq,
			std::vector<complex<double> > &szq);

void mc_pyrochlore_get_thermal_config_and_time_evolve(double spin, double deltat, double omega, double tottime, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, STensor &smunu); 

void mc_pyrochlore_ptall(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double tminKelvin, double tmaxKelvin, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void mc_pyrochlore_pt(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void iterative_pyrochlore(double spin, int L,int nsamples, int nburn, string start_config, 
		   	  double hx, double hy, double hz, 
		   	  double J1, double J2, double J3, double J4, double Jnnn,
		   	  double disorder_strength,
		   	  double gxy, double gz, 
		   	  double & eavg, 
		   	  double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   	  double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void mc_pyrochlore(double spin, int L,int nsamples, int nburn, string start_config, 
		   string mcmove, double temp, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn, double disorder,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

void mc_pyrochlore_groundstate(double spin, int L,int64_t nsamples, int64_t nburn, int nwait, 
			       string start_config, 
		   	       string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   	       double J1, double J2, double J3, double J4, double Jnnn,
		   	       double disorder_strength,
			       double gxy, double gz, 
			       double & eavg, 
			       double &mxavg, double &myavg, double &mzavg, double &e2avg, 
			       double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);



void mc_pyrochlore_slowcool(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);


void mc_pyrochlore_pt(double spin, int L,int nsamples, int nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs);

#endif
