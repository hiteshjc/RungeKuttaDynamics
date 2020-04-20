#ifndef NUMBER_FUNCTIONS_HEADER
#define NUMBER_FUNCTIONS_HEADER

#include"global.h"

int binary_search(int64_t det, std::vector<int64_t> &sorted_dets_list); 

//////////////////////////////////////////////////////////////////////////////
// Useful number functions  
//////////////////////////////////////////////////////////////////////////////
void load_wf(string mfile,
	     std::vector<double> &wf);
void dump_wfs(string mfile,
	      std::vector<double> &eigs,
	      const std::vector< std::vector<double> > &wfs);
int btest(int a,int pos);
int ibset(int a,int pos);
int ibclr(int a,int pos);
int btest64(int64_t a,int pos);
int64_t ibset64(int64_t a,int pos);
int64_t ibclr64(int64_t a,int pos);

int n_choose_k(int n,int k);

void constrained_dets(int num_sites,int num_ones,std::vector<int> &dets);
void constrained_dets_i64(int num_sites,int num_ones,std::vector<int64_t> &dets);

std::vector<int> convert_ind_to_vec(int index, 
				    std::vector<int> nstates);
int convert_vec_to_ind(std::vector<int> const &vec, 
		       std::vector<int> nstates);
void get_adj_list_from_pairs(std::vector< std::vector<int> > const &pairs,
			     std::vector< std::vector<int> > &adj_list);
void convert_num_to_vec(int num, int base, int num_bits, 
			std::vector<int> &config);
int convert_vec_to_num(std::vector<int> &config, int base);
void sort_by_energy(std::vector<double> &eigs,std::vector<double> &szs);
int str_to_int(std::string str);
double str_to_d(std::string str);
bool str_to_bool(std::string str);
template <class T>
std::string to_string (T const &t);
std::vector<int> convert_string_to_vec(std::string const &s);
std::vector<double> convert_string_to_vec_double(std::string const &s);
std::vector< std::vector<int> > convert_string_to_vec_of_vec(std::string const &s);
std::string dtos(double dbl);
int64_t representative(int64_t det, std::vector< std::vector<int> > maps);
void representativeud(int64_t &detu, int64_t &detd, std::vector< std::vector<int> > &maps, int64_t &repu, int64_t &repd);

std::vector< std::vector<int> > make_maps(std::vector<int> T1, std::vector<int> T2);
RMatrix translate(std::vector<int> T, int ntimes);
std::vector<int> make_vector_from_tmat(RMatrix Tmat);
//bool isrepud(int64_t &detu, int64_t &detd, std::vector< std::vector<int> > &maps);
bool isrepud(int64_t &detu, int64_t &detd, std::vector<int> &locsup,std::vector<int> &locsdn, std::vector< std::vector<int> > &maps);

#endif

