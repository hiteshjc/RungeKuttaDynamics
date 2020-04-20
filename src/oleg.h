#ifndef OLEG_MODEL_HEADER
#define OLEG_MODEL_HEADER

///////////////////////////////////////////////////////////////////////////////////////////
//
//                                 OLEG MODEL
//
//////////////////////////////////////////////////////////////////////////////////////////

#include"hamiltonian.h"
class Spin_Half_Oleg: public Ham
{
    public:
    int hilbert; 

    void init(std::vector< std::vector<int> > hexagons,
    	      std::vector< std::vector<int> > up_triangles)
    {
        this->hexagons=hexagons;
        this->up_triangles=up_triangles;
	int max=0;
        for (int i=0;i<hexagons.size();i++)
        {
            for (int j=0;j<6;j++)
            {
                if (hexagons[i][j]>max) {max=hexagons[i][j];}
            }
        }
        this->num_sites=max+1;
	this->hilbert=pow(2,this->num_sites);
    }
    
    void operator()(int spin_det,
		    std::vector<int> &new_spin_dets,
                    std::vector< complex<double> > &hints_list);
    
    Ham* clone() const
    {
        return new Spin_Half_Oleg(*this);
    }    
};

void oleg_setup(std::string filename, 
               Spin_Half_Oleg &oleg);


#endif
