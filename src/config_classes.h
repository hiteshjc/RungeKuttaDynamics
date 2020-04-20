#ifndef CONFIG_CLASSES_HEADER
#define CONFIG_CLASSES_HEADER

#include"global.h"
#include"math_utilities.h"

using namespace std;

class Spin_Config
{
	public:
	int    nsites;
	double spin_length;
	RMatrix spins;
        std::vector< std::vector<int> > neighbors;	
        std::vector< std::vector<int> > pairs;	
        double E_spin;
	double J;

	void init(int nsites, std::vector< std::vector<int> > pairs, 
		  std::vector< std::vector<int> > neighbors, 
		  double spin_length, double J)
	{
	    this->pairs=pairs;
	    this->neighbors=neighbors;
	    this->J = J;
	    this->nsites=nsites;
	    this->spin_length=spin_length;
	    this->spins.resize(this->nsites,3);
	    for (int i=0;i<this->nsites;i++)
	    {
		for (int j=0;j<3;j++) this->spins(i,j)=2.0*(uniform_rnd()-0.5);
	    };
	    this->normalize();
	    this->E_spin=this->energy_spins(J);
	}
	
	void init_afm(int nsites, std::vector< std::vector<int> > pairs, 
		  std::vector< std::vector<int> > neighbors, 
		  std::vector<int> eta, double spin_length, double J)
	{
	    this->pairs=pairs;
	    this->neighbors=neighbors;
	    this->J = J;
	    this->nsites=nsites;
	    this->spin_length=spin_length;
	    this->spins.resize(this->nsites,3);
	    for (int i=0;i<this->nsites;i++)
	    {
		this->spins(i,0)=0.0;
		this->spins(i,1)=0.0;
		this->spins(i,2)=double(eta[i]);
	    };
	    this->normalize();
	    this->E_spin=this->energy_spins(J);
	}

	void init_random_ising(int nsites, std::vector< std::vector<int> > pairs, 
		  std::vector< std::vector<int> > neighbors, 
		  std::vector<int> eta, double spin_length, double J)
	{
	    this->pairs=pairs;
	    this->neighbors=neighbors;
	    this->J = J;
	    this->nsites=nsites;
	    this->spin_length=spin_length;
	    this->spins.resize(this->nsites,3);
	    for (int i=0;i<this->nsites;i++)
	    {
		this->spins(i,0)=0.0;
		this->spins(i,1)=0.0;
		this->spins(i,2)=2.0*(uniform_rnd()-0.5);
	    };
	    this->normalize();
	    this->E_spin=this->energy_spins(J);
	}

	void normalize()
	{
	    	for (int i=0;i<this->nsites;i++)
	    	{
		      double norm=0.0;
		      for (int j=0;j<3;j++) norm=norm+(this->spins(i,j)*this->spins(i,j));
		      norm=sqrt(norm);
		      for (int j=0;j<3;j++) this->spins(i,j)=this->spins(i,j)*this->spin_length/norm;
	    	}
	}

    	double energy_spins(double J)
    	{
		this->E_spin=0.0;
		for (int i=0;i<this->pairs.size();i++)
		{
			double dot=0.0;
			int one=this->pairs[i][0];
			int two=this->pairs[i][1];
			for (int j=0;j<3;j++) {dot+=this->spins(one,j)*this->spins(two,j);}
			this->E_spin+=(J*dot);
		}
    	}    
	void move_random_spin(double delta)
	{
		cout<<"E_spin before move "<<this->E_spin<<endl;
		int site=uniform_rand_int(0,this->nsites);
		double part=0.0;
		for (int i=0;i<this->neighbors[site].size();i++)
		{
			for (int j=0;j<3;j++) part+=this->spins(site,j)*this->spins(this->neighbors[site][i],j);
		}
		part=part*this->J;
		this->E_spin=this->E_spin-part;	// subtract contribution	
		double rnd_x=2.0*(uniform_rnd()-0.5);
		double rnd_y=2.0*(uniform_rnd()-0.5);
		double rnd_z=2.0*(uniform_rnd()-0.5);
		this->spins(site,0)+=(rnd_x*delta);
		this->spins(site,1)+=(rnd_y*delta);
		this->spins(site,2)+=(rnd_z*delta);
		double norm=0.0;
		for (int j=0;j<3;j++) norm=norm+(this->spins(site,j)*this->spins(site,j));
		norm=sqrt(norm);
		for (int j=0;j<3;j++) this->spins(site,j)=this->spins(site,j)*this->spin_length/norm;
		part=0.0;
		for (int i=0;i<this->neighbors[site].size();i++)
		{
			for (int j=0;j<3;j++) part+=this->spins(site,j)*this->spins(neighbors[site][i],j);
		}
		part=part*this->J;
		this->E_spin=this->E_spin+part;	// add contribution	
		cout<<"E_spin after move "<<this->E_spin<<endl;
	}

	void move_random_ising_spin()
	{
		cout<<"E_spin before move "<<this->E_spin<<endl;
		int site=uniform_rand_int(0,this->nsites);
		double part=0.0;
		for (int i=0;i<this->neighbors[site].size();i++)
		{
			for (int j=0;j<3;j++) part+=this->spins(site,j)*this->spins(this->neighbors[site][i],j);
		}
		part=part*this->J;
		this->E_spin=this->E_spin-part;	// subtract contribution	
		this->spins(site,2)=-this->spins(site,2);
		part=0.0;
		for (int i=0;i<this->neighbors[site].size();i++)
		{
			for (int j=0;j<3;j++) part+=this->spins(site,j)*this->spins(neighbors[site][i],j);
		}
		part=part*this->J;
		this->E_spin=this->E_spin+part;	// add contribution	
		cout<<"E_spin after move "<<this->E_spin<<endl;
	}
};

class Fermion_Config
{
	public:
	int    nd;
	int    nsites;
        int    n_double;
	int    nup_holes;
        int    ndn_holes;
	std::vector<int> upsigns;
	std::vector<int> dnsigns;
        std::vector<int> uphole_locations;
        std::vector<int> dnhole_locations;
        int64_t uphole_det;
        int64_t dnhole_det;
	int i0,i1;
		
};

class Hints
{


};


#endif
