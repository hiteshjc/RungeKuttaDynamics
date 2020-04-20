#ifndef PRINTINGFUNCTIONS_H
#define PRINTINGFUNCTIONS_H

#include"global.h"
#include"matrix.h"

///////////////////////////////////////////////////////////////////////////////
// Aesthetic/ Printing functions
///////////////////////////////////////////////////////////////////////////////

template<class T>
void print_vec(T &vec, bool vert=false);

void print_mat_int(std::vector< std::vector<int> > const &mat);

void print_vec_acc(std::vector<double> &vec,bool vert=false,int size=0);

void print_mat_double(std::vector< std::vector<double> > const &mat);

template<class T>
void print_mat(T &mat);

void print_real_mat(RMatrix &mat);

void print_mathematica_pairs(std::vector< std::vector<int> > const &pairs);

void print_mathematica_vector(std::vector<double> const &vec);

#endif

