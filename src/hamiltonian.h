#ifndef HAMILTONIAN_HEADER
#define HAMILTONIAN_HEADER

#include"global.h"
#include"config_classes.h"
using namespace std;

class Ham{

    public:
    int num_sites;
    std::string ham_name;
    std::vector< std::vector<int> > pairs_list;
    std::vector< std::vector<int> > hexagons;
    std::vector< std::vector<int> > up_triangles;
    std::vector<int> eta;
    virtual void operator()
                       (std::vector<int> const &config, 
                        std::vector< std::vector<int> > &touched_sites_list, 
                        std::vector< std::vector<int> > &vals_on_touched_list, 
                        std::vector< complex<double> > &hints_list) {};
    
    virtual void operator()
                       (int config, 
                        std::vector<int> &new_states, 
                        std::vector< complex<double> > &hints_list) {};
    
    virtual void operator()
                    (std::vector<int> const &config_spins,
                     std::vector<int> const &config_up_holes,
                     std::vector<int> const &config_down_holes,
                     std::vector< std::vector<int> > &touched_sites_list_spins, 
                     std::vector< std::vector<int> > &vals_on_touched_list_spins,
                     std::vector< std::vector<int> > &touched_sites_list_upholes, 
                     std::vector< std::vector<int> > &vals_on_touched_list_upholes,
                     std::vector< std::vector<int> > &touched_sites_list_dnholes, 
                     std::vector< std::vector<int> > &vals_on_touched_list_dnholes,
                     std::vector<double> &hints_list){};
    
    virtual void operator()
                    (std::vector<int> const &config_spins,
                     std::vector<int> const &config_up_holes,
                     std::vector<int> const &config_down_holes,
                     std::vector< std::vector<int> > &touched_sites_list_spins, 
                     std::vector< std::vector<int> > &vals_on_touched_list_spins,
                     std::vector< std::vector<int> > &touched_sites_list_upholes, 
                     std::vector< std::vector<int> > &vals_on_touched_list_upholes,
                     std::vector< std::vector<int> > &touched_sites_list_dnholes, 
                     std::vector< std::vector<int> > &vals_on_touched_list_dnholes,
                     std::vector< complex<double> >  &hints_list){};
    virtual void operator()
                   (int up_hole_det,
                    int dn_hole_det,
                    std::vector<int> &new_uphole_dets,
                    std::vector<int> &new_dnhole_dets,
                    std::vector< complex<double> > &hints_list){};


    virtual void operator()
                   (int spin_det,
                    int up_hole_det,
                    int dn_hole_det,
                    std::vector<int> &new_spin_dets,
                    std::vector<int> &new_uphole_dets,
                    std::vector<int> &new_dnhole_dets,
                    std::vector< complex<double> > &hints_list){};

    virtual void operator()(Spin_Config    &sc,
		    	    Fermion_Config &fc,
                    	    Hints          &hints){};

    virtual void operator()
                    (Spin_Config    		    &sc,
		     Fermion_Config 		    &fc,
                     std::vector<int>               &new_configs,
		     std::vector< complex<double> > &hints_list){};


    virtual void init(){};
    
    virtual Ham* clone() const=0;     
            
};

#endif
