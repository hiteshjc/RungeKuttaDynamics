#include"mc.h"
#include"printing_functions.h"
using namespace std;

double local_j_energy(int i, int spinval, std::vector<int> &config, std::vector< std::vector<int> > &neighbors)
{
	// Ferro Potts
	// 0 - > +x
	// 1 - > -x
	// 2 - > +y
	// 3 - > -y
	// 4 - > +z
	// 5 - > -z

	double e=0;
	for (int j=0;j<neighbors[i].size();j++)
	{
		int k=neighbors[i][j];
		if (spinval==config[k]) e=e-1;
		if (abs(spinval-config[k])==1) 
		{
			if ( (spinval+config[k])==1 or (spinval+config[k])==5 or (spinval+config[k])==9)
			{
				e=e+1;
			}
		}
	}
	return e; // Pairs only counted once 
}

//double min_distance_diamond(int L,std::vector<double> &coord1, std::vector<double> &coord2)
double min_distance_diamond(int L, double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dmin=0.0;
	double deltax=(x2-x1)*(x2-x1);
	double deltay=(y2-y1)*(y2-y1);
	double deltaz=(z2-z1)*(z2-z1);
	dmin=sqrt(deltax + deltay + deltaz);
	//double n1=coord1[0]
	//double deltax=min(L-abs(n1-n1p),abs(n1-n1p))/2.0 + min(L-abs(n2-n2p),abs(n2-n2p))/2.0 + (shift);
	//double deltay=min(L-abs(n1-n1p),abs(n1-n1p))/2.0 + min(L-abs(n3-n3p),abs(n3-n3p))/2.0 ;
	//double deltaz=min(L-abs(n2-n2p),abs(n2-n2p))/2.0 + min(L-abs(n3-n3p),abs(n3-n3p))/2.0 ;
	//dmin=sqrt(dmin);
	return dmin;

}

double total_j_energy(std::vector<int> &config, std::vector< std::vector<int> > &neighbors)
{
	// Ferro Potts
	// 0 - > +x
	// 1 - > -x
	// 2 - > +y
	// 3 - > -y
	// 4 - > +z
	// 5 - > -z

	double e=0;
	for (int i=0;i<neighbors.size();i++)
	{
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			if (config[i]==config[k]) e=e-1;
			if (abs(config[i]-config[k])==1) 
			{
				if ((config[i]+config[k])==1 or (config[i]+config[k])==5 or (config[i]+config[k])==9)
				{
					e=e+1;
				}
			}
		}
	}
	return e/2;  // Pairs counted twice so divide by 2
}

std::vector<double> local_magnetization(int value)
{
	// Ferro Potts
	// 0 - > +x
	// 1 - > -x
	// 2 - > +y
	// 3 - > -y
	// 4 - > +z
	// 5 - > -z
	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	if (value==0) mx+=1;
	if (value==1) mx-=1;
	if (value==2) my+=1;
	if (value==3) my-=1;
	if (value==4) mz+=1;
	if (value==5) mz-=1;
	m.push_back(mx);
	m.push_back(my);
	m.push_back(mz);
	return m;
}

std::vector<double> total_magnetization(std::vector<int> &config)
{
	// 0 - > +x
	// 1 - > -x
	// 2 - > +y
	// 3 - > -y
	// 4 - > +z
	// 5 - > -z

	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	for (int i=0;i<config.size();i++)
	{
		if (config[i]==0) mx+=1;
		if (config[i]==1) mx-=1;
		if (config[i]==2) my+=1;
		if (config[i]==3) my-=1;
		if (config[i]==4) mz+=1;
		if (config[i]==5) mz-=1;
	}
	m.push_back(mx);
	m.push_back(my);
	m.push_back(mz);
	return m;
}

std::vector<int> make_random_config(string lattice, int Lx,int Ly,int Lz, int nstates)
{	
	int nsites=Lx*Ly*Lz;
	if (lattice=="cubic")   nsites=Lx*Ly*Lz;
	if (lattice=="diamond") nsites=2*Lx*Ly*Lz; // 2 FCC  // depends on how L is defined for diamond
	std::vector<int> config(nsites);
	
	//for (int i=0;i<nsites;i++) config[i]=uniform_rand_int(0,nstates);
	for (int i=0;i<nsites;i++) config[i]=0;
	return config;
}

void get_neighbor(int Lx,int Ly, int Lz, string direction, int x, int y, int z, int &xnew, int &ynew, int &znew, int &c)
{
	xnew=x;
	ynew=y;
	znew=z;

	if (direction=="+x") xnew=(x+1)%Lx;
	if (direction=="-x") xnew=(x-1+Lx)%Lx;
	if (direction=="+y") ynew=(y+1)%Ly;
	if (direction=="-y") ynew=(y-1+Ly)%Ly;
	if (direction=="+z") znew=(z+1)%Lz;
	if (direction=="-z") znew=(z-1+Lz)%Lz;

	c=(Ly*Lz*xnew) + (Lz*ynew) + znew;

}
	
void mc(string lattice, int nstates, int L, double temp, double hx, double hy, double hz, double & eavg, 
	double &mxavg, double &myavg, double &mzavg, double &e2avg, 
	double &mx2avg, double &my2avg, double &mz2avg)
{
	double beta=1.0/temp;
	// LxLxL pbc lattice
	int nsamples=1000000000;
	int nburn=nsamples/2;
	double nmeas=0;
	int nsites=L*L*L;
	if (lattice=="cubic") nsites=L*L*L;
	//if (lattice=="diamond") nsites=8*L*L*L; // 2 FCC // depends on how L is defined
	if (lattice=="diamond") nsites=2*L*L*L; // 2 FCC If L is defined as the max length of n1, n2, n3
						// The lattice vectors are n1a1+n2a2+n3a3
	std::vector<int> magnetizations(nsites);
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<double> > coords;	
	RMatrix dmat(nsites,nsites);
	
	if (lattice=="cubic")
	{
		// Make nearest neighbors
		for (int i=0;i<nsites;i++)
		{
			int x=i/(L*L);
			int y = (i - (x*L*L))/L;
			int z= (i-(x*L*L)-(y*L));
			int xnew1,xnew2,xnew3,xnew4,xnew5,xnew6;
			int ynew1,ynew2,ynew3,ynew4,ynew5,ynew6;
			int znew1,znew2,znew3,znew4,znew5,znew6;
			int c1,c2,c3,c4,c5,c6;

			get_neighbor(L,L,L,"+x",x,y,z,xnew1,ynew1,znew1,c1);
			//cout<<"+x "<<xnew1<<"  "<<ynew1<<"  "<<znew1<<"  "<<c1<<endl;	
			get_neighbor(L,L,L,"-x",x,y,z,xnew2,ynew2,znew2,c2);	
			//cout<<"-x "<<xnew2<<"  "<<ynew2<<"  "<<znew2<<"  "<<c2<<endl;	
			get_neighbor(L,L,L,"+y",x,y,z,xnew3,ynew3,znew3,c3);	
			//cout<<"+y "<<xnew3<<"  "<<ynew3<<"  "<<znew3<<"  "<<c3<<endl;	
			get_neighbor(L,L,L,"-y",x,y,z,xnew4,ynew4,znew4,c4);	
			//cout<<"-y "<<xnew4<<"  "<<ynew4<<"  "<<znew4<<"  "<<c4<<endl;	
			get_neighbor(L,L,L,"+z",x,y,z,xnew5,ynew5,znew5,c5);	
			//cout<<"+z "<<xnew5<<"  "<<ynew5<<"  "<<znew5<<"  "<<c5<<endl;	
			get_neighbor(L,L,L,"-z",x,y,z,xnew6,ynew6,znew6,c6);	
			//cout<<"-z "<<xnew6<<"  "<<ynew6<<"  "<<znew6<<"  "<<c6<<endl;	
			neighbors[i].push_back(c1);
			neighbors[i].push_back(c2);
			neighbors[i].push_back(c3);
			neighbors[i].push_back(c4);
			neighbors[i].push_back(c5);
			neighbors[i].push_back(c6);
		}
	}
	if (lattice=="diamond")
	{
		int site=0;
		// Make lattice sites
		for (int shift=0;shift<2;shift++)
		{
			for (int n1=0;n1<L;n1++)
			{
				for (int n2=0;n2<L;n2++)
				{
					for (int n3=0;n3<L;n3++)
					{
						double x=double(n1+n2)/2.0;
						double y=double(n1+n3)/2.0;
						double z=double(n2+n3)/2.0;	
						//if (x<L and y<L and z<L)
						{
							std::vector<double> coord;
							coord.push_back(x+double(shift)/4.0);
							coord.push_back(y+double(shift)/4.0);
							coord.push_back(z+double(shift)/4.0);
							coords.push_back(coord);
							site+=1;
						}
					}
				}	
			}
		}
		cout<<"Number of sites recorded is = "<<site<<endl;
		for (int i=0;i<nsites;i++)
		{
			for (int j=0;j<nsites;j++)
			{
				double x1=coords[i][0];
				double y1=coords[i][1];
				double z1=coords[i][2];
				double x2=coords[j][0];
				double y2=coords[j][1];
				double z2=coords[j][2];
				dmat(i,j)=min_distance_diamond(L,x1,y1,z1,x2,y2,z2);
				if (dmat(i,j)<0.44 and abs(dmat(i,j))>0.43) // Just nearest neighbor, do not include site itself as its neighbor
				{
					neighbors[i].push_back(j);
				}
			}
		}
		//print_real_mat(dmat);
	}
	int four=0;
	for (int i=0;i<neighbors.size();i++)
	{
		cout<<" i = "<<i<<", nsize ="<<neighbors[i].size()<<endl;
		if (neighbors[i].size()==4) four+=1;
	}
	cout<<" Number of sites with 4 neighbors = "<<four<<endl;
	cout<<endl;
	// Make random configuration 
	std::vector<int> config=make_random_config(lattice,L,L,L,nstates); 
	std::vector<double> magnetization=total_magnetization(config);
	double energy=total_j_energy(config,neighbors);
	cout<<"Total j energy ="<<energy<<endl;
	energy=energy-(hx*magnetization[0])-(hy*magnetization[1])-(hz*magnetization[2]); // Mag energy
	cout<<"Total j energy + total h energy ="<<energy<<endl;
	for (int i=0;i<config.size();i++) cout<<config[i]<<" ";
	cout<<endl;
	// energy
	// magnetization 
	double mx=magnetization[0];
	double my=magnetization[1];
	double mz=magnetization[2];
	double etot=0.0;
	double e2tot=0.0;
	double mxtot=0.0;
	double mytot=0.0;
	double mztot=0.0;
	double mx2tot=0.0;
	double my2tot=0.0;
	double mz2tot=0.0;

	// Accept reject Metropolis
	for (int i=0; i<(nsamples+nburn);i++)
	{
		//std::vector<int> newconfig=config;
		int site=uniform_rand_int(0,nsites);
		int val=config[site];
		int initial=val;
		while (val==config[site] or val%2!=0) // This latter condition imposes an infinite h model
		{	
			val=uniform_rand_int(0,nstates);	
		}
		// Change spin on site
		config[site]=val;
		//cout<<"Site    = "<<site<<endl;
		//cout<<"Initial = "<<initial<<endl;
		//cout<<"Val     = "<<val<<endl;
		std::vector<double> mag1=local_magnetization(initial);
		std::vector<double> mag2=local_magnetization(val);
		double local_energy1=local_j_energy(site,initial,config,neighbors);
		double local_energy2=local_j_energy(site,val,config,neighbors);

		double mxdiff=mag2[0]-mag1[0];
		double mydiff=mag2[1]-mag1[1];
		double mzdiff=mag2[2]-mag1[2];
		double mnewx=mx+mxdiff;
		double mnewy=my+mydiff;
		double mnewz=mz+mzdiff;
		double ediff=(-hx*mxdiff)+(-hy*mydiff)+(-hz*mzdiff)+(local_energy2-local_energy1);
		double prob=exp(-beta*ediff);
		double rand=uniform_rnd();
		if (rand<prob) // Accept
		{
			config[site]=val;
			//config=newconfig;
			double enew=energy+ediff;
			if (i>nburn and i%nsites==0)
			{
				etot=etot+enew;
				e2tot=e2tot+pow(enew,2.0);
				mxtot=mxtot+mnewx;
				mytot=mytot+mnewy;
				mztot=mztot+mnewz;
				mx2tot=mx2tot+pow(mnewx,2.0);
				my2tot=my2tot+pow(mnewy,2.0);
				mz2tot=mz2tot+pow(mnewz,2.0);
		       	        nmeas=nmeas+1;
			}
			energy=enew;
			mx=mnewx;
			my=mnewy;
			mz=mnewz;
		}
		else
		{
		       config[site]=initial;
		       if (i>nburn and i%nsites==0)
		       {
			       etot=etot+(energy);
			       e2tot=e2tot+pow((energy),2.0);
			       mx2tot=mx2tot+pow(mx,2.0);
			       my2tot=my2tot+pow(my,2.0);
			       mz2tot=mz2tot+pow(mz,2.0);
			       mxtot=mxtot+mx;
			       mytot=mytot+my;
			       mztot=mztot+mz;
			       nmeas=nmeas+1;
		       }
		}
		
	}
	eavg=etot/double(nmeas);
	e2avg=e2tot/double(nmeas);
	mxavg=mxtot/double(nmeas);
	myavg=mytot/double(nmeas);
	mzavg=mztot/double(nmeas);
	mx2avg=mx2tot/double(nmeas);
	my2avg=my2tot/double(nmeas);
	mz2avg=mz2tot/double(nmeas);
	double spheat=(e2avg-(eavg*eavg))/(temp*temp);
	double spheatpersite=(e2avg-(eavg*eavg))/(temp*temp*nsites);

	cout<<"mxavg = "<<boost::format("%+ .15f") %mxavg<<endl;
	cout<<"myavg = "<<boost::format("%+ .15f") %myavg<<endl;
	cout<<"mzavg = "<<boost::format("%+ .15f") %mzavg<<endl;
	cout<<"mx2avg = "<<boost::format("%+ .15f") %mx2avg<<endl;
	cout<<"my2avg = "<<boost::format("%+ .15f") %my2avg<<endl;
	cout<<"mz2avg = "<<boost::format("%+ .15f") %mz2avg<<endl;
	cout<<"eavg = "<<boost::format("%+ .15f") %eavg<<endl;
	cout<<"e2avg = "<<boost::format("%+ .15f") %e2avg<<endl;
	cout<<"Cv     = "<<boost::format("%+ .15f") %spheat<<endl;
	cout<<"Cvps   = "<<boost::format("%+ .15f") %spheatpersite<<endl;
}

