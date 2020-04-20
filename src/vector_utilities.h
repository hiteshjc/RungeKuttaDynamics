#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include"global.h"
using namespace std;

void extract_subset(std::vector<int> const &config, 
                    std::vector<int> const &touched_sites, 
                    std::vector<int> &sub_config);

bool check_subset(std::vector<int> const &vec1, 
                  std::vector<int> const &vec2);

void sort_vecs(std::vector< std::vector<int> > &lists);

void remove_subsets(std::vector< std::vector<int> > &lists);

#endif
