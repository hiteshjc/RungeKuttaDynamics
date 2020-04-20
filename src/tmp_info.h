#ifndef TMP_INFO_HEADER
#define TMP_INFO_HEADER

#include"global.h"
#include"matrix.h"

using namespace std;

class TI
{
	public:
	Matrix								couplings;
	std::vector< std::vector<int> > 				list_of_converted_vecs;
	std::vector<int> 		      				connected_site_locations;
	std::vector< std::vector<int> >      				connected_site_systems;
	std::vector< std::vector< std::vector<int> > >     		connected_system_system;
	std::vector< std::vector<int> >      				connected_site_sys_env;
	std::vector< std::vector<int> >      				system_pairs;
        std::vector<int>      						block_num_states;
	std::vector<double>   						tmp_eigs;
	std::vector< std::vector<double> >   				eigs;
	std::vector< std::vector<int> > 				map_for_states;
	std::vector< std::vector< std::vector<int> > > 			map_for_hints;
	std::vector< std::vector< std::vector<double> > >		hints;
	int								max_spin_index;
	int 								max_ham_spin_change;
	int 								max_trunc_states;
	int 								lanc_dav_it;
	int 								num_spaces;
	bool								diag;
	std::vector<Matrix>                                             tmp_opt_evecs;
	bool								subspace_check;
	std::vector<int>						target_states;
	std::vector<int>						list;
	std::vector<int>						inverse_sz_map,inverse_subspace_map;
	std::vector<double>						distinct_szs;
	std::vector<double>						all_szs;
	std::vector<int>						sites;

};

class Simulation_Params
{
	public:
	string 								loadwffile;
	string 								wffile;
	int 								iterations;
	int 								how_many_eigenvecs;
	bool		                                                store_ham;
	int 								num_cycles;
	bool								ipr;
	bool								rotate;
};

#endif
