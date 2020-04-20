#ifndef HAMILTONIAN_SPIN_FUNCTIONS_HEADER
#define HAMILTONIAN_SPIN_FUNCTIONS_HEADER

#include"global.h"
using namespace std;
void compute_spsm(int num_sites,
                  std::vector<double> const &eigenvec,
		  RMatrix &si_sj);
void compute_szsz(int num_sites,
                  std::vector<double> const &eigenvec,
		  RMatrix &si_sj);

void calc_hints_xyz(double coupling, 
                   int first, int second, int third,
                   int spin_det,
                   std::vector<int> &new_spin_dets,
                   std::vector< complex<double> > &hints_list);

void calc_hints_xyzxyz(double coupling, 
                   int first, int second, int third, int fourth, int fifth, int sixth, 
                   int spin_det,
                   std::vector<int> &new_spin_dets,
                   std::vector< complex<double> > &hints_list);

void calc_hints_V_only(
                     double V,
		     std::vector< std::vector<int> > const &pairs_list, 
		     int const &uphole_det,
		     int const &dnhole_det,
                     std::vector<int> &new_uphole_dets,
                     std::vector<int> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list);

void calc_hints_V_only_ns(
                     double V,
		     std::vector< std::vector<int> > const &pairs_list, 
		     int64_t const &uphole_det,
		     int64_t const &dnhole_det,
                     std::vector<int64_t> &new_uphole_dets,
                     std::vector<int64_t> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list);


void calc_hints_U_plus_V(
                     double U, double V,
		     int nsites,
		     std::vector< std::vector<int> > const &pairs_list, 
		     std::vector< std::vector<int> > const &neighbors, 
		     int const &spin_det,
		     int const &uphole_det,
		     int const &dnhole_det,
                     std::vector<int> &new_spin_dets,
                     std::vector<int> &new_uphole_dets,
                     std::vector<int> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list);

void calc_hints_szsz_plus_U(
                     double J, double U,
		     int nsites,
		     std::vector< std::vector<int> > const &pairs_list, 
		     int const &spin_det,
		     int const &uphole_det,
		     int const &dnhole_det,
                     std::vector<int> &new_spin_dets,
                     std::vector<int> &new_uphole_dets,
                     std::vector<int> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list);

void calc_hints_sxsx_sysy(double coupling, 
                         int first, int second, 
                         int const &spin_det,
			 int const &uphole_det,
			 int const &dnhole_det,
                         std::vector<int> &new_spin_dets,
                         std::vector<int> &new_uphole_dets,
                         std::vector<int> &new_dnhole_dets,
                         std::vector< complex<double> > &hints_list);

void calc_hints_fermion_hop_ns(double t, 
                            int first, int second,
			    int sign_up, int sign_dn, 
			    int64_t const &uphole_det,
			    int64_t const &dnhole_det,
                            std::vector<int64_t> &new_uphole_dets,
                            std::vector<int64_t> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list);

void calc_hints_fermion_hop_up_spin_make_assumption(double t, 
                            int first, int second,
			    int sign_up, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list);
	
void calc_hints_fermion_hop_down_spin_make_assumption(double t, 
                            int first, int second,
			    int sign_down, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list);

void calc_hints_fermion_hop(double t, 
                            int first, int second,
			    int sign_up, int sign_dn, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list);

void calc_hints_U_plus_D(
                     double U, double D,
		     int nsites,
		     std::vector< std::vector<int> > const &pairs_list, 
		     std::vector< std::vector<int> > const &neighbors, 
		     std::vector< std::vector<int> > const &neighbors_within_rh, 
		     int const &spin_det,
		     int const &uphole_det,
		     int const &dnhole_det,
                     std::vector<int> &new_spin_dets,
                     std::vector<int> &new_uphole_dets,
                     std::vector<int> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list);

void compute_si_sj_spins(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &si_sj);

void compute_one_rdm_down_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &one_rdm);

void compute_one_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &one_rdm);

void compute_n_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &n_2);


void compute_nup_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &nup_2);

void compute_ndn_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &ndn_2);

void compute_nu_nu(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nu_nu);

void compute_nd_nd(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nd_nd);


void compute_nu_nd(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nu_nd);

void compute_explicit_two_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm);

void compute_explicit_two_rdm_dn_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm);

void compute_explicit_two_rdm_uddu(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm);


void compute_explicit_three_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &three_rdm);

#endif
