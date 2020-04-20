#ifndef MATH_UTIL_HEADER
#define MATH_UTIL_HEADER

#include"global.h"

double uniform_rnd();
int uniform_rand_int(int lowest,int highest);
double mean(double sum_x, int nsamples);
double variance(double sum_x, double sum_x_2, int nsamples);
double sigma(double sum_x, double sum_x_2, int nsamples);
double sigma_over_root_n(double sum_x, double sum_x_2, int nsamples);
double error_in_mean(double sum_x, double sum_x_2, int nsamples);
int closest_int(double);
void convert_to_g6(int nsites,
                   std::vector< std::vector<int> > &ordered_pairs, 
                   std:: string &tmp_string);
void convert_all_to_g6(std::string read_file);


#endif
