#ifndef MC_D_HEADER
#define MC_D_HEADER

#include"global.h"
#include"math_utilities.h"

void mc_finite_D(int L, string lattice, 
		 int nsamples, int nburn, string start_config, 
		 string mcmove, double temp, double hx, double hy, 
		 double hz, double Jval, double D, double Dp,
		 double alphaL,int kcut,
		 double & eavg, double &mxavg, double &myavg, 
		 double &mzavg, double &e2avg, 
		 double &mx2avg, double &my2avg, double &mz2avg);
void normalize(std::vector<double> &vec);
void big_move_continuous_spin(double x, double y, double z, double &valx, double &valy, double &valz);
void infD_move_special_continuous_spin(double x, double y, double z, double &valx, double &valy, double &valz);
void random_move_continuous_spin(double &valx,double &valy, double &valz);
void conical_move_continuous_spin_rnds_provided(double r_max, double rnd1, double rnd2, double x, double y, double z, double &valx, double &valy, double &valz);
void conical_move_continuous_spin(double r_max, double x, double y, double z, double &valx, double &valy, double &valz);
void make_111_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz);
void make_x_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz);
void make_random_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz);
void get_two_normalized_orth_dirs(double x, double y, double z, 
		       double &x1, double &y1, double &z1, 
		       double &x2, double &y2, double &z2);

#endif
