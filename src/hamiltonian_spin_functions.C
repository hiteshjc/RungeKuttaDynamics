#include"hamiltonian_spin_functions.h"
#include"number_functions.h"
using namespace std;

///////////////////////////////////////////////////////////////////////////////
// This is Oleg's W term... . the Hamiltonian will be - Sum W 
// Used sigma instead of S - so factor of 2 per site
void calc_hints_xyz(double coupling, 
                   int first, int second, int third,
                   int spin_det,
                   std::vector<int> &new_spin_dets,
                   std::vector< complex<double> > &hints_list)
{
       	int a=btest(spin_det,first); 
        int b=btest(spin_det,second); 
       	int c=btest(spin_det,third); 
	double sign1=1.0;
	if (b==1) sign1=-1.0;

        // XYZXYZ - only one out of 4 terms survive
	double signsz=(double(c)-0.5);
        hints_list.push_back( (8.0/4.0) * complex<double> (0,-1)* coupling*signsz*sign1);
	int new_spin_det=spin_det;
	new_spin_det=ibclr(new_spin_det,first);
	new_spin_det=ibclr(new_spin_det,second);
	if (a==0) new_spin_det=ibset(new_spin_det,first);
	if (b==0) new_spin_det=ibset(new_spin_det,second);
	new_spin_dets.push_back(new_spin_det);
}

void calc_hints_xyzxyz(double coupling, 
                   int first, int second, int third, int fourth, int fifth, int sixth, 
                   int spin_det,
                   std::vector<int> &new_spin_dets,
                   std::vector< complex<double> > &hints_list)
{
       	int a=btest(spin_det,first); 
        int b=btest(spin_det,second); 
       	int c=btest(spin_det,third); 
       	int d=btest(spin_det,fourth); 
       	int e=btest(spin_det,fifth); 
       	int f=btest(spin_det,sixth); 
	double sign1=1.0;
        double sign2=1.0; 

	if (b==1) sign1=-1.0;
	if (e==1) sign2=-1.0;

        // XYZXYZ - only one out of 16 terms survive
	double signsz=(double(c)-0.5)*(double(f)-0.5);
        hints_list.push_back( (-64.0/16.0) * coupling*signsz*sign1*sign2);

	int new_spin_det=spin_det;
	new_spin_det=ibclr(new_spin_det,first);
	new_spin_det=ibclr(new_spin_det,second);
	new_spin_det=ibclr(new_spin_det,fourth);
	new_spin_det=ibclr(new_spin_det,fifth);
	
	//if (a==1) new_spin_det=ibclr(new_spin_det,first);
	//if (b==1) new_spin_det=ibclr(new_spin_det,second);
	//if (d==1) new_spin_det=ibclr(new_spin_det,fourth);
	//if (e==1) new_spin_det=ibclr(new_spin_det,fifth);
	
	if (a==0) new_spin_det=ibset(new_spin_det,first);
	if (b==0) new_spin_det=ibset(new_spin_det,second);
	if (d==0) new_spin_det=ibset(new_spin_det,fourth);
	if (e==0) new_spin_det=ibset(new_spin_det,fifth);
	
	new_spin_dets.push_back(new_spin_det);
}

void set_signs(int det, std::vector<int> &signs)
{
            int ctru=0;
	    int num_sites=signs.size();
	    for (int i=0;i<num_sites;i++)
	    {
		signs[i]=ctru;	
		if (btest(det,i)==1) ctru=ctru+1;
	    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_sxsx_sysy(double coupling, 
                         int first, int second, 
                         int const &spin_det,
			 int const &uphole_det,
			 int const &dnhole_det,
                         std::vector<int> &new_spin_dets,
                         std::vector<int> &new_uphole_dets,
                         std::vector<int> &new_dnhole_dets,
                         std::vector< complex<double> > &hints_list)
{
	int new_state;
        // Apply S+ S- + S-S_+
        int a,b;
        a=btest(spin_det,first); b=btest(spin_det,second);

        if (a!=b)
	{ 
		if (a==0 and b==1)
		{
		    new_state=ibset(spin_det,first);
		    new_state=ibclr(new_state,second);
		}
		else
		{
		    new_state=ibset(spin_det,second);
		    new_state=ibclr(new_state,first);
		}
	        new_spin_dets.push_back(new_state);
	        new_uphole_dets.push_back(uphole_det);
	        new_dnhole_dets.push_back(dnhole_det);
	        hints_list.push_back(0.5*coupling);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_fermion_hop_ns(double t, 
                            int first, int second,
			    int sign_up, int sign_dn, 
			    int64_t const &uphole_det,
			    int64_t const &dnhole_det,
                            std::vector<int64_t> &new_uphole_dets,
                            std::vector<int64_t> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list)
{
	int64_t new_state;
        // Apply t a_i+ aj+ ai
        int a,b;
	//double mt=t;
        
	a=btest64(uphole_det,first); b=btest64(uphole_det,second);
        if (a!=b)
	{ 
		if (a==1 and b==0)
		{
		    new_state=ibclr(uphole_det,first);
		    new_state=ibset(new_state,second);
		}
		else
		{
		    new_state=ibclr(uphole_det,second);
		    new_state=ibset(new_state,first);
		}
		new_uphole_dets.push_back(new_state);
		new_dnhole_dets.push_back(dnhole_det);
		hints_list.push_back(t*sign_up);
		//hints_list.push_back(t);
	}
	
	a=btest64(dnhole_det,first); b=btest64(dnhole_det,second);
        if (a!=b)
	{ 
		if (a==1 and b==0)
		{
		    new_state=ibclr64(dnhole_det,first);
		    new_state=ibset64(new_state,second);
		}
		else
		{
		    new_state=ibclr64(dnhole_det,second);
		    new_state=ibset64(new_state,first);
		}
		new_dnhole_dets.push_back(new_state);
		new_uphole_dets.push_back(uphole_det);
		hints_list.push_back(t*sign_dn);
		//hints_list.push_back(t);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_fermion_hop_up_spin_make_assumption(double t, 
                            int first, int second,
			    int sign_up, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list)
{
	int new_state;
	int b=btest(uphole_det,second);
	if (b==0)
	{
	    new_state=ibclr(uphole_det,first);
	    new_state=ibset(new_state,second);
	    new_uphole_dets.push_back(new_state);
	    new_dnhole_dets.push_back(dnhole_det);
	    hints_list.push_back(t*sign_up);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_fermion_hop_down_spin_make_assumption(double t, 
                            int first, int second,
			    int sign_down, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list)
{
	int new_state;
	int b=btest(dnhole_det,second);
	if (b==0)
	{
	    new_state=ibclr(dnhole_det,first);
	    new_state=ibset(new_state,second);
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(new_state);
	    hints_list.push_back(t*sign_down);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_fermion_hop(double t, 
                            int first, int second,
			    int sign_up, int sign_dn, 
			    int const &uphole_det,
			    int const &dnhole_det,
                            std::vector<int> &new_uphole_dets,
                            std::vector<int> &new_dnhole_dets,
                            std::vector< complex<double> > &hints_list)
{
	int new_state;
        // Apply t a_i+ aj+ ai
        int a,b;
	//double mt=t;
        
	a=btest(uphole_det,first); b=btest(uphole_det,second);
        if (a!=b)
	{ 
		if (a==1 and b==0)
		{
		    new_state=ibclr(uphole_det,first);
		    new_state=ibset(new_state,second);
		}
		else
		{
		    new_state=ibclr(uphole_det,second);
		    new_state=ibset(new_state,first);
		}
		new_uphole_dets.push_back(new_state);
		new_dnhole_dets.push_back(dnhole_det);
		hints_list.push_back(t*sign_up);
	}
	
	a=btest(dnhole_det,first); b=btest(dnhole_det,second);
        if (a!=b)
	{ 
		if (a==1 and b==0)
		{
		    new_state=ibclr(dnhole_det,first);
		    new_state=ibset(new_state,second);
		}
		else
		{
		    new_state=ibclr(dnhole_det,second);
		    new_state=ibset(new_state,first);
		}
		new_dnhole_dets.push_back(new_state);
		new_uphole_dets.push_back(uphole_det);
		hints_list.push_back(t*sign_dn);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////
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
                     std::vector< complex<double> > &hints_list)
{
            double hint=0.0;

	    if (J!=0)
	    {
		    for (int i=0;i<pairs_list.size();i++) hint+=(J*((double(btest(spin_det,pairs_list[i][0]))-0.5)*(double(btest(spin_det,pairs_list[i][1]))-0.5)));
	    }
           
	    if (U!=0)
	    { 
		    for (int i=0;i<nsites;i++) hint+=U*double(btest(uphole_det,i))*double(btest(dnhole_det,i));
	    }
	    new_spin_dets.push_back(spin_det);
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(dnhole_det);
            hints_list.push_back( complex<double> (hint) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_V_only_ns(
                     double V,
		     std::vector< std::vector<int> > const &pairs_list, 
		     int64_t const &uphole_det,
		     int64_t const &dnhole_det,
                     std::vector<int64_t> &new_uphole_dets,
                     std::vector<int64_t> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list)
{
            double hint=0.0;
	    for (int i=0;i<pairs_list.size();i++) hint+=V*double(btest64(uphole_det,pairs_list[i][0])+btest64(dnhole_det,pairs_list[i][0]))*double(btest64(uphole_det,pairs_list[i][1])+btest64(dnhole_det,pairs_list[i][1]));
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(dnhole_det);
            hints_list.push_back( complex<double> (hint) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void calc_hints_V_only(
                     double V,
		     std::vector< std::vector<int> > const &pairs_list, 
		     int const &uphole_det,
		     int const &dnhole_det,
                     std::vector<int> &new_uphole_dets,
                     std::vector<int> &new_dnhole_dets,
                     std::vector< complex<double> > &hints_list)
{
            double hint=0.0;
	    for (int i=0;i<pairs_list.size();i++) hint+=V*double(btest(uphole_det,pairs_list[i][0])+btest(dnhole_det,pairs_list[i][0]))*double(btest(uphole_det,pairs_list[i][1])+btest(dnhole_det,pairs_list[i][1]));
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(dnhole_det);
            hints_list.push_back( complex<double> (hint) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
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
                     std::vector< complex<double> > &hints_list)
{
            double hint=0.0;

	    for (int i=0;i<nsites;i++) hint+=U*double(btest(uphole_det,i))*double(btest(dnhole_det,i));

	    for (int i=0;i<pairs_list.size();i++) hint+=V*double(btest(uphole_det,pairs_list[i][0]))*double(btest(uphole_det,pairs_list[i][1]));
	    for (int i=0;i<pairs_list.size();i++) hint+=V*double(btest(dnhole_det,pairs_list[i][0]))*double(btest(dnhole_det,pairs_list[i][1]));
	    new_spin_dets.push_back(spin_det);
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(dnhole_det);
            hints_list.push_back( complex<double> (hint) );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
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
                     std::vector< complex<double> > &hints_list)
{
            double hint=0.0;

	    for (int i=0;i<nsites;i++) hint+=U*double(btest(uphole_det,i))*double(btest(dnhole_det,i));

	    if (D!=0)
	    { 
		for (int i=0;i<nsites;i++)
		{
			if (btest(uphole_det,i)==1)
			{
				for (int k=0;k<neighbors_within_rh[i].size();k++)
				{
				int l=neighbors_within_rh[i][k];
				double sh_dot_sl=(btest(spin_det,l)-0.5);
				if (sh_dot_sl<1 and sh_dot_sl>0)
				{
					for (int j=0;j<neighbors[l].size();j++)
					{
						int m=neighbors[l][j];
						if (l<m)
						{
						if (std::count(neighbors_within_rh[i].begin(),neighbors_within_rh[i].end(),m)==1)
						{
						double sh_dot_sm=(btest(spin_det,m)-0.5);
						if (sh_dot_sm<1 and sh_dot_sm>0)
						{
						   hint=hint-(D*sh_dot_sl*sh_dot_sm);
						}
						}
						}
					}
				}
				}
			}
			if (btest(dnhole_det,i)==1)
			{
				for (int k=0;k<neighbors_within_rh[i].size();k++)
				{
				int l=neighbors_within_rh[i][k];
				double sh_dot_sl=-(btest(spin_det,l)-0.5);
				if (sh_dot_sl<1 and sh_dot_sl>0)
				{
					for (int j=0;j<neighbors[l].size();j++)
					{
						int m=neighbors[l][j];
						if(l<m)
						{
						if (std::count(neighbors_within_rh[i].begin(),neighbors_within_rh[i].end(),m)==1)
						{
						double sh_dot_sm=-(btest(spin_det,m)-0.5);
						if (sh_dot_sm<1 and sh_dot_sm>0)
						{
						   hint=hint-(D*sh_dot_sl*sh_dot_sm);
						}
						}
						}
					}
				}
				}
			}
			
		}	
	    }

	    new_spin_dets.push_back(spin_det);
	    new_uphole_dets.push_back(uphole_det);
	    new_dnhole_dets.push_back(dnhole_det);
            hints_list.push_back( complex<double> (hint) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_szsz(int num_sites,
                  std::vector<double> const &eigenvec,
		  RMatrix &si_sj)
{
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
      cout<<"hilbert ="<<hilbert<<endl;  
      //# pragma omp parallel for 
      for (int j=0;j<hilbert;j++)
      {
	    int spin_det=j;
	    for (int m=0;m<num_sites;m++)
	    {
		    int a=btest(spin_det,m); 
		    for (int n=0;n<num_sites;n++)
		    {
			    int b=btest(spin_det,n); 
			    //# pragma omp atomic 
				si_sj(m,n)+=((double(a)-0.5)*(double(b)-0.5)*(eigenvec[j]*eigenvec[j]));
		     }
	    }
     }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_spsm(int num_sites,
                  std::vector<double> const &eigenvec,
		  RMatrix &si_sj)
{
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
    
      cout<<"hilbert ="<<hilbert<<endl;  
      //# pragma omp parallel for 
      for (int j=0;j<hilbert;j++)
      {
	    int spin_det=j;
	    for (int m=0;m<num_sites;m++)
	    {
		    int a=btest(spin_det,m); 
		    for (int n=0;n<num_sites;n++)
		    {
			    int b=btest(spin_det,n); 
			    if (m!=n)
			    {
				    if (a==0 and b==1)
				    {
					int new_spin_det=ibset(spin_det,m);
					new_spin_det=ibclr(new_spin_det,n);
					//# pragma omp atomic 
						si_sj(m,n)+=(0.5*(eigenvec[new_spin_det]*eigenvec[j]));
				    }
				    if (b==0 and a==1)
				    {
					int new_spin_det=ibset(spin_det,n);
					new_spin_det=ibclr(new_spin_det,m);
					//# pragma omp atomic 
						si_sj(m,n)+=(0.5*(eigenvec[new_spin_det]*eigenvec[j]));
				    }
			    }
			    else
			    {
					//# pragma omp atomic 
						si_sj(m,n)+=(0.5*eigenvec[j]*eigenvec[j]);
			    }
		     }
	    }
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_si_sj_spins(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &si_sj)
{
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nd=n_dnhole_dets;
      int 		nud=n_uphole_dets*n_dnhole_dets;
      std::vector<int>  num_states;

      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      si_sj.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) si_sj[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int spin_det=spin_dets[i0];

	    for (int m=0;m<num_sites;m++)
	    {
	    int a=btest(spin_det,m); 
	    for (int n=m;n<num_sites;n++)
	    {
	    int b=btest(spin_det,n); 
	    if (m!=n)
	    {
		    if (a==0 and b==1)
		    {
			int new_spin_det=ibset(spin_det,m);
			new_spin_det=ibclr(new_spin_det,n);
			int loc=(inverse_map_spin[new_spin_det]*(nud))+(i1*(nd))+i2;
			si_sj(m,n)+=(0.5*(eigenvec[loc]*eigenvec[j]));
		    }

		    if (b==0 and a==1)
		    {
			int new_spin_det=ibset(spin_det,n);
			new_spin_det=ibclr(new_spin_det,m);
			int loc=(inverse_map_spin[new_spin_det]*(nud))+(i1*(nd))+i2;
			si_sj(m,n)+=(0.5*(eigenvec[loc]*eigenvec[j]));
		    }
	    }
	    else
	    {si_sj(m,n)+=(0.5*eigenvec[j]*eigenvec[j]);}

	    si_sj(m,n)+=(eigenvec[j]*eigenvec[j]*(double(a)-0.5)*(double(b)-0.5));
	    si_sj(n,m)=si_sj(m,n);
	    }
	    }
     }
}

/////////////////////////////////////////////////////////////////////////////////////////
void compute_one_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &one_rdm)
{
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      one_rdm.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) one_rdm[m]=0.0;

      //#pragma omp parallel for     
      for (int j=0;j<hilbert;j++)
      {
      	    std::vector<int>  upsigns(num_sites);
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
            int ctru=0;
	    for (int i=0;i<num_sites;i++)
	    {
		upsigns[i]=ctru;	
		if (btest(uphole_det,i)==1) ctru=ctru+1;
	    }
		

	    for (int m=0;m<num_sites;m++)
	    {
		    int a=btest(uphole_det,m); 
		    for (int n=m;n<num_sites;n++)
		    {
			    int b=btest(uphole_det,n); 
			    if (m!=n)
			    {
				    int su=1;
				    if (abs(upsigns[n]-upsigns[m+1])%2 !=0 ) su=-1;
				    if (a==0 and b==1)
				    {
					int new_uphole_det=ibset(uphole_det,m);
					new_uphole_det=ibclr(new_uphole_det,n);
					int loc=(i0*nud)+(inverse_map_uphole[new_uphole_det]*(nd))+i2;
					//#pragma omp atomic
						one_rdm(m,n)+=(su*(eigenvec[loc]*eigenvec[j]));
				    }
			    }
			    else
			    {
					//#pragma omp atomic
						one_rdm(m,m)+=(double(a)*eigenvec[j]*eigenvec[j]);
			    }
		    }
	    }
     }
}
//////////////////////////////////////////////////////////////////////////////////////
void compute_n_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &n_2)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      n_2.resize(num_sites);
      for (int m=0;m<num_sites;m++) n_2[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
	    int dnhole_det=dnhole_dets[i2];
            for (int m=0;m<num_sites;m++)
	    {
	    int a=btest(uphole_det,m)+btest(dnhole_det,m); 
	    n_2[m]+=(double(a)*double(a)*eigenvec[j]*eigenvec[j]);
	    }
     }
}

//////////////////////////////////////////////////////////////////////////////////////
void compute_nu_nu(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nu_nu)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      nu_nu.resize(num_sites,num_sites);
      for (int m=0;m<num_sites*num_sites;m++) nu_nu[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
	    int dnhole_det=dnhole_dets[i2];
            for (int m=0;m<num_sites;m++)
	    {
		for (int n=0;n<num_sites;n++)
		{
	    		int a=btest(uphole_det,m); 
	    		int b=btest(uphole_det,n); 
	    		nu_nu(m,n)+=(double(a)*double(b)*eigenvec[j]*eigenvec[j]);
		}
	    }
     }
}

//////////////////////////////////////////////////////////////////////////////////////
void compute_nd_nd(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nd_nd)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      nd_nd.resize(num_sites,num_sites);
      for (int m=0;m<num_sites*num_sites;m++) nd_nd[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
	    int dnhole_det=dnhole_dets[i2];
            for (int m=0;m<num_sites;m++)
	    {
		for (int n=0;n<num_sites;n++)
		{
	    		int a=btest(dnhole_det,m); 
	    		int b=btest(dnhole_det,n); 
	    		nd_nd(m,n)+=(double(a)*double(b)*eigenvec[j]*eigenvec[j]);
		}
	    }
     }
}


//////////////////////////////////////////////////////////////////////////////////////
void compute_nu_nd(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &nu_nd)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      nu_nd.resize(num_sites,num_sites);
      for (int m=0;m<num_sites*num_sites;m++) nu_nd[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
	    int dnhole_det=dnhole_dets[i2];
            for (int m=0;m<num_sites;m++)
	    {
		for (int n=0;n<num_sites;n++)
		{
	    		int a=btest(uphole_det,m); 
	    		int b=btest(dnhole_det,n); 
	    		nu_nd(m,n)+=(double(a)*double(b)*eigenvec[j]*eigenvec[j]);
		}
	    }
     }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_nup_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &nup_2)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      nup_2.resize(num_sites);
      for (int m=0;m<num_sites;m++) nup_2[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1];
            for (int m=0;m<num_sites;m++)
	    {
	    int a=btest(uphole_det,m); 
	    nup_2[m]+=(double(a)*double(a)*eigenvec[j]*eigenvec[j]);
	    }
     }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_ndn_2(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   std::vector<double> &ndn_2)
{
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      ndn_2.resize(num_sites);
      for (int m=0;m<num_sites;m++) ndn_2[m]=0.0;
       
      for (int j=0;j<hilbert;j++)
      {
	    std::vector<int> vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int dnhole_det=dnhole_dets[i2];
            for (int m=0;m<num_sites;m++)
	    {
	    int a=btest(dnhole_det,m); 
	    ndn_2[m]+=(double(a)*double(a)*eigenvec[j]*eigenvec[j]);
	    }
     }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_one_rdm_down_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &one_rdm)
{
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      one_rdm.resize(num_sites,num_sites);
      for (int m=0;m<num_pairs;m++) one_rdm[m]=0.0;
      
     // #pragma omp parallel for 
      for (int j=0;j<hilbert;j++)
      {
      	    std::vector<int>  dnsigns(num_sites);
	    std::vector<int>  vec=convert_ind_to_vec(j,num_states);
	    int i0=vec[0];
	    int i1=vec[1];
	    int i2=vec[2]; 
	    int dnhole_det=dnhole_dets[i2];
            int ctrd=0;
	    for (int i=0;i<num_sites;i++)
	    {
		dnsigns[i]=ctrd;	
		if (btest(dnhole_det,i)==1) ctrd=ctrd+1;
	    }
		
	    for (int m=0;m<num_sites;m++)
	    {
		    int a=btest(dnhole_det,m); 
		    for (int n=m;n<num_sites;n++)
		    {
			    int b=btest(dnhole_det,n); 
			    if (m!=n)
			    {
				    int sd=1;
				    if (abs(dnsigns[n]-dnsigns[m+1])%2 !=0 ) sd=-1;
				    if (a==0 and b==1)
				    {
					int new_dnhole_det=ibset(dnhole_det,m);
					new_dnhole_det=ibclr(new_dnhole_det,n);
					int loc=(i0*nud)+(i1*nd)+inverse_map_dnhole[new_dnhole_det];
					//#pragma omp atomic
						one_rdm(m,n)+=(sd*(eigenvec[loc]*eigenvec[j]));
				    }
			    }
			    else
			    {
					//#pragma omp atomic
						one_rdm(m,m)+=(double(a)*eigenvec[j]*eigenvec[j]);
			    }
	    	       }
     	    }
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_explicit_two_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm)
{
      //< C_i dagger C_j dagger C_l C_k >
      //  No restrictions on i,j,k,l 
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      int dim=num_sites*(num_sites);
      two_rdm.resize(dim,dim);
      for (int m=0;m<dim*dim;m++) two_rdm[m]=0.0;
      
      //#pragma omp parallel for 
      for (int h=0;h<hilbert;h++)
      {
	    std::vector<int>  upsigns(num_sites);
	    std::vector<int>  new_signs(num_sites);
	    std::vector<int>  new_signs2(num_sites);
	    std::vector<int>  new_signs3(num_sites);
	    std::vector<int>  new_signs4(num_sites);
	    std::vector<int>  new_signs5(num_sites);
	    std::vector<int>  vec=convert_ind_to_vec(h,num_states);
	    int i0=vec[0]; int i1=vec[1]; int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1]; set_signs(uphole_det,upsigns);

	    for (int k=0;k<num_sites;k++)
	    {
		int new_det=uphole_det;
	    	int a=btest(new_det,k);
		 
	    	if (a==1)
	    	{
		    int new_det2=ibclr(new_det,k); set_signs(new_det2,new_signs2);
		    int gamma1=1;if (abs(upsigns[k])%2 !=0 ) gamma1=-1;
		    
		    for (int l=0;l<num_sites;l++)
		    {
			int ind2=(num_sites*k)+l; 
			int b=btest(new_det2,l); 
		    	if (b==1)
		    	{
			    int gamma2=1; if (abs(new_signs2[l])%2 !=0 ) gamma2=-1;
			    int new_det3=ibclr(new_det2,l); set_signs(new_det3, new_signs3);

	    		    for (int j=0;j<num_sites;j++)
			    {
	    			int c=btest(new_det3,j);
				if (c==0)
				{
				    int gamma3=1; if (abs(new_signs3[j])%2 !=0 ) gamma3=-1;
				    int new_det4=ibset(new_det3,j); set_signs(new_det4, new_signs4);

				    for (int i=0;i<num_sites;i++)
				    { 
					int ind1=(num_sites*i)+j;
		    			int d=btest(new_det4,i); 
		    			if (d==0)
		    			{
					int gamma4=1; if (abs(new_signs4[i])%2 !=0 ) gamma4=-1;
					int new_det5=ibset(new_det4,i); set_signs(new_det5, new_signs5);
					int new_uphole_det=new_det5;
					int loc=(i0*nud)+(inverse_map_uphole[new_uphole_det]*(nd))+i2;
					//#pragma omp atomic 
						two_rdm(ind1,ind2)+=(double(gamma1*gamma2*gamma3*gamma4)*(eigenvec[loc]*eigenvec[h]));
					}
				    }
				}
			    }
		    	}
		    }
	     	}
	   }
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_explicit_two_rdm_dn_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm)
{
      //< C_i dagger C_j dagger C_l C_k >
      //  No restrictions on i,j,k,l 
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      int dim=num_sites*(num_sites);
      two_rdm.resize(dim,dim);
      for (int m=0;m<dim*dim;m++) two_rdm[m]=0.0;
      
      //#pragma omp parallel for 
      for (int h=0;h<hilbert;h++)
      {
	    std::vector<int>  dnsigns(num_sites);
	    std::vector<int>  new_signs(num_sites);
	    std::vector<int>  new_signs2(num_sites);
	    std::vector<int>  new_signs3(num_sites);
	    std::vector<int>  new_signs4(num_sites);
	    std::vector<int>  new_signs5(num_sites);
	    std::vector<int> vec=convert_ind_to_vec(h,num_states);
	    int i0=vec[0]; int i1=vec[1]; int i2=vec[2]; 
	    int dnhole_det=dnhole_dets[i2]; set_signs(dnhole_det,dnsigns);

	    for (int k=0;k<num_sites;k++)
	    {
		int new_det=dnhole_det;
	    	int a=btest(new_det,k);
		 
	    	if (a==1)
	    	{
		    int new_det2=ibclr(new_det,k); set_signs(new_det2,new_signs2);
		    int gamma1=1;if (abs(dnsigns[k])%2 !=0 ) gamma1=-1;
		    
		    for (int l=0;l<num_sites;l++)
		    {
			int ind2=(num_sites*k)+l; 
			int b=btest(new_det2,l); 
		    	if (b==1)
		    	{
			    int gamma2=1; if (abs(new_signs2[l])%2 !=0 ) gamma2=-1;
			    int new_det3=ibclr(new_det2,l); set_signs(new_det3, new_signs3);

	    		    for (int j=0;j<num_sites;j++)
			    {
	    			int c=btest(new_det3,j);
				if (c==0)
				{
				    int gamma3=1; if (abs(new_signs3[j])%2 !=0 ) gamma3=-1;
				    int new_det4=ibset(new_det3,j); set_signs(new_det4, new_signs4);

				    for (int i=0;i<num_sites;i++)
				    { 
					int ind1=(num_sites*i)+j;
		    			int d=btest(new_det4,i); 
		    			if (d==0)
		    			{
					int gamma4=1; if (abs(new_signs4[i])%2 !=0 ) gamma4=-1;
					int new_det5=ibset(new_det4,i); set_signs(new_det5, new_signs5);
					int new_dnhole_det=new_det5;
					int loc=(i0*nud)+(i1*nd)+(inverse_map_dnhole[new_dnhole_det]);
					//#pragma omp atomic 
						two_rdm(ind1,ind2)+=(double(gamma1*gamma2*gamma3*gamma4)*(eigenvec[loc]*eigenvec[h]));
					}
				    }
				}
			    }
		    	}
		    }
	     	}
	   }
     }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_explicit_two_rdm_uddu(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &two_rdm)
{
      //< C_i(u) dagger C_j(d) dagger C_l(d) C_k(u) >
      //  No restrictions on i,j,k,l 
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      int dim=num_sites*(num_sites);
      two_rdm.resize(dim,dim);
      for (int m=0;m<dim*dim;m++) two_rdm[m]=0.0;

      //# pragma omp parallel for       
      for (int h=0;h<hilbert;h++)
      {
	    std::vector<int>  upsigns(num_sites),dnsigns(num_sites);
	    std::vector<int>  new_signs(num_sites);
	    std::vector<int>  new_signs2(num_sites);
	    std::vector<int>  new_signs3(num_sites);
	    std::vector<int>  new_signs4(num_sites);
	    std::vector<int>  new_signs5(num_sites);
	    std::vector<int> vec=convert_ind_to_vec(h,num_states);
	    int i0=vec[0]; int i1=vec[1]; int i2=vec[2]; 
	    int uphole_det=uphole_dets[i1]; set_signs(uphole_det,upsigns);
	    int dnhole_det=dnhole_dets[i2]; set_signs(dnhole_det,dnsigns);

	    for (int k=0;k<num_sites;k++)
	    {
		int new_det=uphole_det;
	    	int a=btest(new_det,k);
		 
	    	if (a==1)
	    	{
		    int new_det2=ibclr(new_det,k); set_signs(new_det2,new_signs2);
		    int gamma1=1;if (abs(upsigns[k])%2 !=0 ) gamma1=-1;
		    
		    for (int l=0;l<num_sites;l++)
		    {
			int ind2=(num_sites*k)+l; 
			int b=btest(dnhole_det,l); 
		    	if (b==1)
		    	{
			    int gamma2=1; if (abs(dnsigns[l])%2 !=0 ) gamma2=-1;
			    int new_det3=ibclr(dnhole_det,l); set_signs(new_det3, new_signs3);

	    		    for (int j=0;j<num_sites;j++)
			    {
	    			int c=btest(new_det3,j);
				if (c==0)
				{
				    int gamma3=1; if (abs(new_signs3[j])%2 !=0 ) gamma3=-1;
				    int new_det4=ibset(new_det3,j); set_signs(new_det4, new_signs4);

				    for (int i=0;i<num_sites;i++)
				    { 
					int ind1=(num_sites*i)+j;
		    			int d=btest(new_det2,i); 
		    			if (d==0)
		    			{
					int gamma4=1; if (abs(new_signs2[i])%2 !=0 ) gamma4=-1;
					int new_det5=ibset(new_det2,i); set_signs(new_det5, new_signs5);
					int new_uphole_det=new_det5;
					int new_dnhole_det=new_det4;
					int loc=(i0*nud)+(inverse_map_uphole[new_uphole_det]*nd)+(inverse_map_dnhole[new_dnhole_det]);
					//#pragma omp atomic
						two_rdm(ind1,ind2)+=(double(gamma1*gamma2*gamma3*gamma4)*(eigenvec[loc]*eigenvec[h]));
					}
				    }
				}
			    }
		    	}
		    }
	     	}
	   }
     }
}



/////////////////////////////////////////////////////////////////////////////////////////
void compute_explicit_three_rdm_up_electrons(int num_sites,
                   std::vector<double> const &eigenvec,
   	 	   std::vector<int> const &spin_dets,
   	 	   std::vector<int> const &uphole_dets,
   	 	   std::vector<int> const &dnhole_dets,
   	 	   std::vector<int> const &inverse_map_spin,
   	 	   std::vector<int> const &inverse_map_uphole,
   	 	   std::vector<int> const &inverse_map_dnhole,
		   RMatrix &three_rdm)
{
      //< C_i dagger C_j dagger C_k dagger C_n C_m C_l >
      //  No restrictions on i,j,k,l,m,n 
      int 	   	num_pairs=num_sites*num_sites;
      int 		hilbert=eigenvec.size();
      int               n_uphole_dets=uphole_dets.size();
      int 		n_dnhole_dets=dnhole_dets.size();
      int 		nud=n_uphole_dets*n_dnhole_dets;
      int 		nd=n_dnhole_dets;
      std::vector<int>  num_states;
      std::vector<int>  upsigns(num_sites);
      std::vector<int>  new_signs(num_sites);
      std::vector<int>  new_signs2(num_sites);
      std::vector<int>  new_signs3(num_sites);
      std::vector<int>  new_signs4(num_sites);
      std::vector<int>  new_signs5(num_sites);
      std::vector<int>  new_signs6(num_sites);
      std::vector<int>  new_signs7(num_sites);
      
      num_states.push_back(spin_dets.size());
      num_states.push_back(n_uphole_dets);
      num_states.push_back(n_dnhole_dets);

      int dim=num_sites*num_sites*num_sites;
      three_rdm.resize(dim,dim);
      for (int m=0;m<dim*dim;m++) three_rdm[m]=0.0;
       
    for (int h=0;h<hilbert;h++)
    {
    std::vector<int> vec=convert_ind_to_vec(h,num_states);
    int i0=vec[0]; int i1=vec[1]; int i2=vec[2]; 
    int uphole_det=uphole_dets[i1]; set_signs(uphole_det,upsigns);

    for (int l=0;l<num_sites;l++)
    {
	int new_det=uphole_det;
	int a=btest(new_det,l);
	if (a==1)
	{
	    int new_det2=ibclr(new_det,l); set_signs(new_det2,new_signs2);
	    int gamma1=1; if (abs(upsigns[l])%2 !=0 ) gamma1=-1;
	    for (int m=0;m<num_sites;m++)
	    {
		int b=btest(new_det2,m); 
		if (b==1)
		{
		    int gamma2=1; if (abs(new_signs2[m])%2 !=0 ) gamma2=-1;
		    int new_det3=ibclr(new_det2,m); set_signs(new_det3,new_signs3);
		    for (int n=0;n<num_sites;n++)
		    {
			int ind2=(num_sites*num_sites*l)+(num_sites*m)+n;
			int c=btest(new_det3,n);
			if (c==1)
			{
			    int gamma3=1; if (abs(new_signs3[n])%2 !=0 ) gamma3=-1;
			    int new_det4=ibclr(new_det3,n); set_signs(new_det4,new_signs4);
			    for (int k=0;k<num_sites;k++)
			    { 
				int d=btest(new_det4,k); 
				if (d==0)
				{
				    int gamma4=1;if (abs(new_signs4[k])%2 !=0 ) gamma4=-1;
				    int new_det5=ibset(new_det4,k); set_signs(new_det5,new_signs5);
				    for (int j=0;j<num_sites;j++)
				    {
					int e=btest(new_det5,j); 
					if (e==0)
					{
					  int gamma5=1;if (abs(new_signs5[j])%2 !=0 ) gamma5=-1;
					  int new_det6=ibset(new_det5,j); set_signs(new_det6,new_signs6);
					  for (int i=0;i<num_sites;i++)
					  {
						int ind1=(num_sites*num_sites*i)+(num_sites*j)+k;
						int f=btest(new_det6,i);
						if (f==0)
						{
				    		 int gamma6=1;if (abs(new_signs6[i])%2 !=0 ) gamma6=-1;
				    		 int new_det7=ibset(new_det6,i); set_signs(new_det7,new_signs7);
						 int new_uphole_det=new_det7;
						 int loc=(i0*nud)+(inverse_map_uphole[new_uphole_det]*(nd))+i2;
						 three_rdm(ind1,ind2)+=(double(gamma1*gamma2*gamma3*gamma4*gamma5*gamma6)*(eigenvec[loc]*eigenvec[h]));
						} 
					  }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
   }
   }
}
