#include"classical_spin_hole_model.h"
#include"config_classes.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////
//                                 Classical Spin Hole Model
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Classical_Spin_Hole_Model::operator()
                    (Spin_Config    		    &sc,
		     Fermion_Config 		    &fc,
                     std::vector<int>               &new_configs,
		     std::vector< complex<double> > &hints_list)
{
	new_configs.clear();
	hints_list.clear();

	int ind;
	double U=this->U;
	double D=this->D;
	int i0=fc.i0;
	int i1=fc.i1;
	int nd=fc.nd;
	int nsites=this->num_sites;
	int64_t new_det;

	// Set fermion configs
	fc.n_double=0;
	fc.upsigns.resize(nsites);
	fc.dnsigns.resize(nsites);
	int ctru=0;
	for (int i=0;i<this->num_sites;i++)
	{
		fc.upsigns[i]=ctru;	
		if (btest64(fc.uphole_det,i)==1) ctru=ctru+1;
	}
	fc.nup_holes=ctru;
		
	int ctrd=0;
	for (int i=0;i<this->num_sites;i++)
	{
		fc.dnsigns[i]=ctrd;	
		if (btest64(fc.dnhole_det,i)==1) ctrd=ctrd+1;
	}
	fc.ndn_holes=ctrd;

	for (int i=0;i<nsites;i++) 
	{
		int a=btest64(fc.uphole_det,i);
		int b=btest64(fc.dnhole_det,i);
		if (a==1) fc.uphole_locations.push_back(i);
		if (b==1) fc.dnhole_locations.push_back(i);
		fc.n_double+=(a*b);
	}
	// J Si dot Sj - always saved as part of spin onfig
	double E_spin=sc.E_spin; // this is constant 
	//cout<<"E_spin="<<E_spin<<endl;
	
	// D f(S_hole dot Si, S_hole dot Sj) - This term is classical
	double spin_hole_energy=0.0;
	for (int k=0;k<fc.nup_holes;k++)
	{
		int i=fc.uphole_locations[k];
		for (int k=0;k<this->neighbors_within_rh[i].size();k++)
		{
		int l=this->neighbors_within_rh[i][k];
		double sh_dot_sl=0.5*sc.spins(l,2);
		if (sh_dot_sl<1 and sh_dot_sl>0)
		{
			for (int j=0;j<sc.neighbors[l].size();j++)
			{
				int m=sc.neighbors[l][j];
				if (std::count(this->neighbors_within_rh[i].begin(),this->neighbors_within_rh[i].end(),m)==1)
				{
					double sh_dot_sm=0.5*sc.spins(m,2);
					if (sh_dot_sm<1 and sh_dot_sm>0)
					{
				   	spin_hole_energy-=(D*sh_dot_sl*sh_dot_sm);
					}
				}
			}
		}
		}
	}
	
	for (int k=0;k<fc.ndn_holes;k++)
	{
		int i=fc.dnhole_locations[k];
		for (int k=0;k<this->neighbors_within_rh[i].size();k++)
		{
		int l=this->neighbors_within_rh[i][k];
		double sh_dot_sl=-0.5*sc.spins(l,2);
		if (sh_dot_sl<1 and sh_dot_sl>0)
		{
			for (int j=0;j<sc.neighbors[l].size();j++)
			{
				int m=sc.neighbors[l][j];
				if (std::count(this->neighbors_within_rh[i].begin(),this->neighbors_within_rh[i].end(),m)==1)
				{
					double sh_dot_sm=-0.5*sc.spins(m,2);
					if (sh_dot_sm<1 and sh_dot_sm>0)
					{
				   	spin_hole_energy-=(D*sh_dot_sl*sh_dot_sm);
					}
				}
			}
		}
		}
	}
	

	double hubbard_U=U*double(fc.n_double);                // Hubbard U
	double diag_term=hubbard_U+spin_hole_energy+E_spin;    // Diagonal term

	hints_list.push_back(diag_term);
	new_configs.push_back((i0*nd)+i1);

	// Hubbard t
	if (this->t!=0.0)
	{
		for (int i=0;i<this->pairs.size();i++)
		{
			int su=1;
			int sd=1;
			int first=this->pairs[i][0];
			int second=this->pairs[i][1];

			int mx=max(first,second);
			int mn=min(first,second);
			if (abs(fc.upsigns[mx]-fc.upsigns[mn+1])%2 !=0 ) su=-1;
			if (abs(fc.dnsigns[mx]-fc.dnsigns[mn+1])%2 !=0 ) sd=-1;
	
			int a=btest64(fc.uphole_det,first);
			int b=btest64(fc.uphole_det,second);	
			if (a!=b)
			{
			    if (a==0 and b==1)
			    { 
				new_det=fc.uphole_det;
				new_det=ibset64(new_det,first);
				new_det=ibclr64(new_det,second);
			    }
			    else
			    {
				new_det=fc.uphole_det;
				new_det=ibclr64(new_det,first);
				new_det=ibset64(new_det,second);
			    }
			    ind = binary_search(new_det,this->uphole_dets);
			    new_configs.push_back((ind*nd)+i1);
			    hints_list.push_back(-t*su);
			}
			a=btest64(fc.dnhole_det,first);
			b=btest64(fc.dnhole_det,second);	
			if (a!=b)
			{
			    if (a==0 and b==1)
			    { 
				new_det=fc.dnhole_det;
				new_det=ibset64(new_det,first);
				new_det=ibclr64(new_det,second);
			    }
			    else
			    {
				new_det=fc.dnhole_det;
				new_det=ibclr64(new_det,first);
				new_det=ibset64(new_det,second);
			    }
			    ind = binary_search(new_det,this->dnhole_dets);
			    new_configs.push_back((i0*nd)+ind);
			    hints_list.push_back(-t*sd);
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////////////////////////
void Classical_Spin_Hole_Model::set_pairs_square_lattice()
{
	std::vector<int> pair(2);
	int l_x=this->l_x;
	int l_y=this->l_y;
	bool pbc=this->pbc;
	this->pairs.clear();
	this->neighbors.clear();
	this->neighbors_within_rh.clear();
 	this->eta.clear();

	for(int i=0;i<this->num_sites;i++) 
	{
		this->neighbors.push_back(std::vector<int>());
		this->neighbors_within_rh.push_back(std::vector<int>());
	}
	
 	// Nearest neighbor pairs 
	for(int i=0;i<this->num_sites;i++)
	{
		int xi=i%l_x;
		int yi=(i-xi)/l_x;
		if ( (xi+yi) %2 ==0 ) {this->eta.push_back(1);}
		else                  {this->eta.push_back(-1);}

		int ind_l,ind_r,ind_u,ind_d;

		// Right
		if (xi==l_x-1 and pbc)
		{ 
			int xr=0; 
			ind_r=(yi*l_x) + xr;
			if (i<ind_r){pair[0]=i;pair[1]=ind_r;this->pairs.push_back(pair);}
		}
		
		if (xi<l_x-1)
		{
			int xr=xi+1;
			ind_r=(yi*l_x) + xr;
			if (i<ind_r){pair[0]=i;pair[1]=ind_r;this->pairs.push_back(pair);}
		}	
	
		// Left
		if (xi==0 and pbc)
		{ 
			int xl=l_x-1; 
			ind_l=(yi*l_x) + xl;
			if (i<ind_l and ind_l!=ind_r){pair[0]=i;pair[1]=ind_l;this->pairs.push_back(pair);}
		}
		
		if (xi>0)
		{
			int xl=xi-1;
			ind_l=(yi*l_x) + xl;
			if (i<ind_l and ind_l!=ind_r){pair[0]=i;pair[1]=ind_l;this->pairs.push_back(pair);}
		}	

		// Up 
		if (yi==l_y-1 and pbc)
		{ 
			int yu=0; 
			ind_u=(yu*l_x) + xi;
			if (i<ind_u){pair[0]=i;pair[1]=ind_u;this->pairs.push_back(pair);}
		}
		
		if (yi<l_y-1)
		{
			int yu=yi+1;
			ind_u=(yu*l_x) + xi;
			if (i<ind_u){pair[0]=i;pair[1]=ind_u;this->pairs.push_back(pair);}
		}	
	
		// Down	
		if (yi==0 and pbc)
		{ 
			int yd=l_y-1; 
			ind_d=(yd*l_x) + xi;
			if (i<ind_d and ind_d!=ind_u){pair[0]=i;pair[1]=ind_d;this->pairs.push_back(pair);}
		}
		
		if (yi>0)
		{
			int yd=yi-1;
			ind_d=(yd*l_x) + xi;
			if (i<ind_d and ind_d!=ind_u){pair[0]=i;pair[1]=ind_d;this->pairs.push_back(pair);}
		}	
	}

	// Nearest neighbors from pairs constructed
	for (int i=0;i<this->pairs.size();i++)
	{
		this->neighbors[pairs[i][0]].push_back(pairs[i][1]);
		this->neighbors[pairs[i][1]].push_back(pairs[i][0]);
	}	

	// Set distance matrix
	this->distance.resize(this->num_sites,this->num_sites);
	
	for (int i=0;i<this->num_sites;i++)
	{
		int xi=i%l_x;
		int yi=(i-xi)/l_x;
		for (int j=0;j<this->num_sites;j++)
		{
			int xj=j%l_x;
			int yj=(j-xj)/l_x;
			if (pbc)
			{
				int d=0,dmin=0;
				for (int a=-1;a<=1;a++)
				{
					for (int b=-1;b<=1;b++)
					{
					   d=((xi-xj-a*l_x)*(xi-xj-a*l_x)) + ((yi-yj-b*l_y)*(yi-yj-b*l_y));
					   if (d<dmin or (a==-1 and b==-1)) dmin=d; 
					}
				}
				this->distance(i,j)=sqrt(double(dmin));
			}
			else
			{
				this->distance(i,j)=double((xi-xj)*(xi-xj)) + double((yi-yj)*(yi-yj));
				this->distance(i,j)=sqrt(this->distance(i,j));
			}
		}
	}
	// Look at distance matrix
	for (int i=0;i<this->num_sites;i++)
	{
		for (int j=0;j<this->num_sites;j++)
		{
			if (this->distance(i,j)<=(rh+1.0e-6)) neighbors_within_rh[i].push_back(j); 
		}
	}

	// Print everything
	cout<<"Pairs"<<endl;
	print_mathematica_pairs(this->pairs);
	cout<<endl;
	cout<<"Neighbors"<<endl;
	print_mat_int(this->neighbors);
	cout<<endl;
	cout<<"Neighbors with rh="<<rh<<endl;
	print_mat_int(this->neighbors_within_rh);
	cout<<endl;
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////
void classical_spin_hole_setup(string filename, 
               Classical_Spin_Hole_Model &classical_spin_hole)

{    
    ///////////////////////////////////////////////////
    // Set the Classical spin quantum hole Hamiltonian
    ///////////////////////////////////////////////////

    bool found=true,pbc;
    string str,str_ret;
    int l_x,l_y,nup_holes,ndn_holes;
    double J,t,D,U,rh,J_h;

    // coupling parameters search
    search_for("l_x",filename,str_ret,found);
    if (found){l_x=str_to_int(str_ret);} else {l_x=1;}

    search_for("l_y",filename,str_ret,found);
    if (found){l_y=str_to_int(str_ret);} else {l_y=1;}
    
    search_for("nup_holes",filename,str_ret,found);
    if (found){nup_holes=str_to_int(str_ret);} else{nup_holes=0;}
    
    search_for("ndn_holes",filename,str_ret,found);
    if (found){ndn_holes=str_to_int(str_ret);} else{ndn_holes=0;}
    
    search_for("J",filename,str_ret,found);
    if (found){J=str_to_d(str_ret);} else{J=1.0;}
    
    search_for("D",filename,str_ret,found);
    if (found){D=str_to_d(str_ret);} else{D=0.0;}
    
    search_for("t",filename,str_ret,found);
    if (found){t=str_to_d(str_ret);} else{t=0.0;}
    
    search_for("pbc",filename,str_ret,found);
    if (found){pbc=str_to_bool(str_ret);} else{pbc=true;}
 
    search_for("U",filename,str_ret,found);
    if (found){U=str_to_d(str_ret);} else{U=0.0;}
    
    search_for("J_h",filename,str_ret,found);
    if (found){J_h=str_to_d(str_ret);} else{J_h=0.0;}
    
    search_for("rh",filename,str_ret,found);
    if (found){rh=str_to_d(str_ret);} else{rh=1.0;}

    if (nup_holes>l_x*l_y or ndn_holes>l_x*l_y) {cout<<"nholes (UP or DOWN ) can NOT exceed number of sites"<<endl;cout<<endl;exit(1);                                                }    
    
    classical_spin_hole.init(l_x,l_y,pbc,nup_holes,ndn_holes,J,D,t,U,J_h,rh);
	
}


