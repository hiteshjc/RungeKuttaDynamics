#ifndef CLASSICAL_SPIN_HOLE_MODEL_HEADER
#define CLASSICAL_SPIN_HOLE_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 CLASSICAL SPIN HOLE MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////
#include"hamiltonian.h"
#include"global.h"
#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
#include"search_for.h"
#include"printing_functions.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////// 
class Classical_Spin_Hole_Model: public Ham
{
    public:
    int 			    l_x,l_y;
    int 			    nup_holes,ndn_holes,nup_spins;
    std::vector<int64_t>            uphole_dets;
    std::vector<int64_t>            dnhole_dets;
    int 			    sz;
    double 			    J,D,t,U,rh,J_h;
    std::vector< std::vector<int> > pairs;
    std::vector< std::vector<int> > neighbors;
    std::vector< std::vector<int> > neighbors_within_rh;
    bool			    pbc;
    int 			    hilbert_upholes;
    int 			    hilbert_dnholes;
    int 			    hilbert;
    RMatrix                         distance;

    public:
    void init(int l_x,int l_y, bool pbc, int nup_holes,int ndn_holes, 
	      double J, double D, double t, double U, double J_h,double rh)
    {
	this->l_x=l_x;
	this->l_y=l_y;
	this->pbc=pbc;
	this->rh=rh;
	this->num_sites=(this->l_x*this->l_y);
	this->nup_holes=nup_holes;
	this->ndn_holes=ndn_holes;
	this->sz=sz;
	this->nup_spins=(this->num_sites/2)+(this->sz);
	this->hilbert_upholes=n_choose_k(this->num_sites,this->nup_holes);
	this->hilbert_dnholes=n_choose_k(this->num_sites,this->ndn_holes);
	this->hilbert=this->hilbert_upholes*this->hilbert_dnholes;
        this->J=J;this->D=D;this->t=t;this->U=U; this->J_h=J_h;
	this->set_pairs_square_lattice();
        cout<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"l_x        ="<<this->l_x<<endl;
        cout<<"l_y        ="<<this->l_y<<endl;
        cout<<"pbc        ="<<this->pbc<<endl;
        cout<<"N_sites    ="<<this->num_sites<<endl;
        cout<<"N_bonds    ="<<this->pairs.size()<<endl;
        cout<<"nup_holes  ="<<this->nup_holes<<endl;
        cout<<"ndn_holes  ="<<this->ndn_holes<<endl;
        cout<<"J          ="<<this->J<<endl;
        cout<<"D          ="<<this->D<<endl;
        cout<<"J_h        ="<<this->J_h<<endl;
        cout<<"t          ="<<this->t<<endl;
        cout<<"U          ="<<this->U<<endl;
        cout<<"rh         ="<<this->rh<<endl;
        cout<<endl;
        cout<<"--------------------------------------------"<<endl;
	cout<<endl;
    }

/////////////////////////////////////////////////////////////////////////// 
   void operator()
                    (Spin_Config    		    &sc,
		     Fermion_Config 		    &fc,
                     std::vector<int>               &new_configs,
		     std::vector< complex<double> > &hints_list);

    //void operator()(Spin_Config    &sc,
//		    Fermion_Config &fc,
 //                   Hints          &hints);

    void set_pairs_square_lattice();

    Ham* clone() const
    {
	return new Classical_Spin_Hole_Model(*this);
    }    
};

void classical_spin_hole_setup(std::string filename, 
               Classical_Spin_Hole_Model &spin_hole);


#endif
