#include"mc_finite_D.h"
#include"printing_functions.h"
using namespace std;

double bfn(double r, double alpha)
{
	double invrootpi=0.56418958354;
	double alphar=alpha*r;
	double alphar2=alphar*alphar;
	double val1=erfc(alphar);
	double val2=(2.0*alphar*invrootpi)*exp(-alphar2);
	double val=(val1+val2)/(r*r*r);
	return val;
}

double cfn(double r, double alpha)
{
	double invrootpi=0.56418958354;
	double alphar=alpha*r;
	double alphar2=alphar*alphar;
	double val1=3.0*erfc(alphar);
	double val2=(2*alphar*invrootpi)*(3+(2.0*alphar2))*exp(-alphar2);
	double val=(val1+val2)/(r*r*r*r*r);
	return val;
}

void get_b_c(double r, double alpha, double &b, double &c)
{
	double invrootpi=0.56418958354;
	double alphar=alpha*r;
	double alphar2=alphar*alphar;
	double val1=erfc(alphar);
	double val2=(2.0*alphar*invrootpi);
	double val3=exp(-alphar2);
	double val4=(3.0+(2.0*alphar2));
	double r3=(r*r*r);
	double r5=(r3*r*r);
	b = (val1 + (val2*val3))/r3;
	c = ((3.0*val1) + (val2*val3*val4))/r5;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_bmat_cmat(RMatrix &dmat, double alpha, RMatrix &bmat, RMatrix &cmat, 
		    std::vector<int> &n1vec, 
		    std::vector<int> &n2vec, 
		    std::vector<int> &n3vec)
{
	int nsites=dmat.NRows();
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++)
		{
			if (i!=j)
			{
				double b,c;
				//if (n1vec[i]!=n1vec[j] or n2vec[i]!=n2vec[j] or n3vec[i]!=n3vec[j]) // Different units
				{
					get_b_c(dmat(i,j),alpha,b,c);
					bmat(i,j)=b;
					cmat(i,j)=c;
				}
				/*else // If 2 atoms in same "unit" then treat dipolar exactly 
				{
					get_b_c(dmat(i,j),0.0,b,c);
					bmat(i,j)=b;
					cmat(i,j)=c;
				}*/
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void normalize(std::vector<double> &vec)
{

	double norm=0;
	int size=vec.size();
	for (int i=0;i<size;i++) norm=norm+(vec[i]*vec[i]);
	norm=1.0/sqrt(norm);
	for (int i=0;i<size;i++) vec[i]=(vec[i]*norm);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void big_move_continuous_spin(double x, double y, double z, double &valx, double &valy, double &valz)
{
	int temp=uniform_rand_int(0,16);
	if(temp==0) {valx=y;valy=x;valz=z;}
	if(temp==1) {valx=-y;valy=x;valz=z;}
	if(temp==2) {valx=-y;valy=-x;valz=z;}
	if(temp==3) {valx=y;valy=-x;valz=z;}
	if(temp==4) {valx=-x;valy=y;valz=z;}
	if(temp==5) {valx=x;valy=-y;valz=z;}
	if(temp==6) {valx=x;valy=z;valz=y;}
	if(temp==7) {valx=x;valy=z;valz=-y;}
	if(temp==8) {valx=x;valy=-z;valz=-y;}
	if(temp==9) {valx=x;valy=-z;valz=y;}
	if(temp==10) {valx=x;valy=y;valz=-z;}
	if(temp==11) {valx=z;valy=y;valz=x;}
	if(temp==12) {valx=z;valy=y;valz=-x;}
	if(temp==13) {valx=-z;valy=y;valz=-x;}
	if(temp==14) {valx=-z;valy=y;valz=x;}
	if(temp==15) {valx=-x;valy=-y;valz=-z;}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void infD_move_special_continuous_spin(double x, double y, double z, double &valx, double &valy, double &valz)
{
	int oldspin,newspin;
	valx=0;
	valy=0;
	valz=0;

	if (x>1.0e-10)        oldspin=0;
	if (x<-1.0e-10)       oldspin=1;
	if (y>1.0e-10)        oldspin=2;
	if (y<-1.0e-10)       oldspin=3;
	if (z>1.0e-10)        oldspin=4;
	if (z<-1.0e-10)       oldspin=5;

	newspin=oldspin;

	while (newspin==oldspin)
	{
		newspin=uniform_rand_int(0,6);
	}

	if(newspin==0) valx=1;
	if(newspin==1) valx=-1;
	if(newspin==2) valy=1;
	if(newspin==3) valy=-1;
	if(newspin==4) valz=1;
	if(newspin==5) valz=-1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void random_move_continuous_spin(double &valx,double &valy, double &valz)
{
	bool cond=false;
	// Choose a completely random direction on the sphere - This is INEFFICENT at low temps
	while (cond==false)
	{
		valx=2*uniform_rnd() - 1;	
		valy=2*uniform_rnd() - 1;	
		valz=2*uniform_rnd() - 1;
		if (valx*valx + valy*valy + valz*valz <=1) cond=true; 
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_two_normalized_orth_dirs(double x, double y, double z, 
		       double &x1, double &y1, double &z1, 
		       double &x2, double &y2, double &z2)
{
	std::vector<double> d1vec(3),d2vec(3);
	if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=-y;
		d1vec[1]=x;
		d1vec[2]=0.0;
		d2vec[0]=1.0;
		d2vec[1]=(y/x);
		d2vec[2]=(-(x/z) - (y*y/(x*z)));
	}
	else if (abs(y)>1.0e-16 and abs(z)>1.0e-16 and abs(x)<1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=z;
		d2vec[2]=-y;
	}
	else if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=z;
		d2vec[1]=0.0;
		d2vec[2]=-x;
	}
	else if (abs(x)>1.0e-16 and abs(y)>1.0e-16 and abs(z)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=0.0;
		d1vec[2]=1.0;
		d2vec[0]=y;
		d2vec[1]=-x;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(y)<1.0e-16 and abs(z)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=1.0;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(z)<1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	else if (abs(y)<1.0e-16 and abs(z)<1.0e-16 and abs(x)>1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	normalize(d1vec);
	normalize(d2vec);

	x1=d1vec[0];
	y1=d1vec[1];
	z1=d1vec[2];
	
	x2=d2vec[0];
	y2=d2vec[1];
	z2=d2vec[2];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Problematiccccccccccccccc    routineeeeeeee
void conical_move_continuous_spin_rnds_provided(double r_max, double rnd1, double rnd2, double x, double y, double z, double &valx, double &valy, double &valz)
{
	// Instead we want to choose something CLOSE to the direction we have
	// It should not be too close either 
	bool cond=true;
	std::vector<double> d1vec(3),d2vec(3);
	double d1,d2;

	//while (cond==false) // Problem if fixed random numbers given !!!!!!!!!
	//{
		d1=r_max*(2.0*rnd1 - 1);	
		d2=r_max*(2.0*rnd2 - 1);	
		//ASSUMED if (d1*d1 + d2*d2 <=(r_max*r_max)) cond=true; 
	//}

	if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=-y;
		d1vec[1]=x;
		d1vec[2]=0.0;
		d2vec[0]=1.0;
		d2vec[1]=(y/x);
		d2vec[2]=(-(x/z) - (y*y/(x*z)));
	}
	else if (abs(y)>1.0e-16 and abs(z)>1.0e-16 and abs(x)<1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=z;
		d2vec[2]=-y;
	}
	else if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=z;
		d2vec[1]=0.0;
		d2vec[2]=-x;
	}
	else if (abs(x)>1.0e-16 and abs(y)>1.0e-16 and abs(z)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=0.0;
		d1vec[2]=1.0;
		d2vec[0]=y;
		d2vec[1]=-x;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(y)<1.0e-16 and abs(z)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=1.0;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(z)<1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	else if (abs(y)<1.0e-16 and abs(z)<1.0e-16 and abs(x)>1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	normalize(d1vec);
	normalize(d2vec);

	valx=x+(d1*d1vec[0])+(d2*d2vec[0]);	
	valy=y+(d1*d1vec[1])+(d2*d2vec[1]);	
	valz=z+(d1*d1vec[2])+(d2*d2vec[2]);	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void conical_move_continuous_spin(double r_max, double x, double y, double z, double &valx, double &valy, double &valz)
{
	// Instead we want to choose something CLOSE to the direction we have
	// It should not be too close either 
	bool cond=false;
	std::vector<double> d1vec(3),d2vec(3);
	double d1,d2;
	while (cond==false)
	{
		d1=r_max*(2*uniform_rnd() - 1);	
		d2=r_max*(2*uniform_rnd() - 1);	
		if (d1*d1 + d2*d2 <=(r_max*r_max)) cond=true; 
	}

	if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=-y;
		d1vec[1]=x;
		d1vec[2]=0.0;
		d2vec[0]=1.0;
		d2vec[1]=(y/x);
		d2vec[2]=(-(x/z) - (y*y/(x*z)));
	}
	else if (abs(y)>1.0e-16 and abs(z)>1.0e-16 and abs(x)<1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=z;
		d2vec[2]=-y;
	}
	else if (abs(x)>1.0e-16 and abs(z)>1.0e-16 and abs(y)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=z;
		d2vec[1]=0.0;
		d2vec[2]=-x;
	}
	else if (abs(x)>1.0e-16 and abs(y)>1.0e-16 and abs(z)<1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=0.0;
		d1vec[2]=1.0;
		d2vec[0]=y;
		d2vec[1]=-x;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(y)<1.0e-16 and abs(z)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=1.0;
		d2vec[2]=0.0;
	}
	else if (abs(x)<1.0e-16 and abs(z)<1.0e-16 and abs(y)>1.0e-16)
	{
		d1vec[0]=1.0;
		d1vec[1]=0.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	else if (abs(y)<1.0e-16 and abs(z)<1.0e-16 and abs(x)>1.0e-16)
	{
		d1vec[0]=0.0;
		d1vec[1]=1.0;
		d1vec[2]=0.0;
		d2vec[0]=0.0;
		d2vec[1]=0.0;
		d2vec[2]=1.0;
	}
	normalize(d1vec);
	normalize(d2vec);

	valx=x+(d1*d1vec[0])+(d2*d2vec[0]);	
	valy=y+(d1*d1vec[1])+(d2*d2vec[1]);	
	valz=z+(d1*d1vec[2])+(d2*d2vec[2]);	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double j_energy_diff(double coupling,int i, double x, double y, double z, double valx, double valy, double valz,
		      std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors)
{
	double e=0;
	for (int j=0;j<neighbors[i].size();j++)
	{
		int k=neighbors[i][j];
		e=e-((valx-x)*configx[k]);
		e=e-((valy-y)*configy[k]);
		e=e-((valz-z)*configz[k]);
	}
	return e*coupling; // Pairs only counted once 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double local_j_energy(int i, double x, double y, double z,
		      std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors)
{
	double e=0;
	for (int j=0;j<neighbors[i].size();j++)
	{
		int k=neighbors[i][j];
		e=e-(x*configx[k]);
		e=e-(y*configy[k]);
		e=e-(z*configz[k]);
	}
	return e; // Pairs only counted once 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_j_energy(double coupling, 
		      std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors)
{
	// Sum over nearest neighbors
	// coupling is > 0 but term is -coupling *spin i * spin j - so ferromagnetic
	double e=0;
	for (int i=0;i<neighbors.size();i++)
	{
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			e=e-(configx[i]*configx[k]);
			e=e-(configy[i]*configy[k]);
			e=e-(configz[i]*configz[k]);
		}
	}
	return (e/2.0)*coupling;  // Pairs counted twice so divide by 2
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void total_dipolar_energy(  double 				coupling,
			    int 				L,
			    std::vector<double> 		&configx, 
			    std::vector<double> 		&configy, 
			    std::vector<double> 		&configz, 
		      	    std::vector< std::vector<double> >  &coords, 
		      	    std::vector< std::vector<int> >     &groups, 
		      	    RMatrix 				&rijxmat, 
		      	    RMatrix 				&rijymat, 
		      	    RMatrix 				&rijzmat, 
		      	    RMatrix 				&bmat, 
		      	    RMatrix 				&cmat, 
		      	    std::vector< std::vector<double> >  &kvecs,
		      	    Matrix 				&expmikr,
			    std::vector<double> 		&kfactors,
			    double 				alpha,
			    double 				&udr,
			    double 				&udk,
			    double 				&usurf,
			    double 				&uself,
			    double 				&total_e,
			    std::vector<complex<double> >	&dmoments,
			    std::vector<double>			&magnetization)
{
	// coupling > 0 and is AFM 
	double invrootpi=0.56418958354;
	int numkvecs=kvecs.size();
	//cout<<"Number of k vectors is = "<<numkvecs<<endl;
	int nsites=coords.size();
	//double twopibv = (2.0*3.14159265359)/double(nsites); - cubic
	double twopibv = (2.0*3.14159265359)/double(L*L*L); //- diamond or cubic
	///////////////////////////////////////////////////////////////
	// Real space part of the energy
	//////////////////////////////////////////////////////////////
	// rij corresponds to the closest image (MINIMUM image convention)
	// All other images contribute nothing (zero) to the real space sum
	// alpha is chosen so that this condition holds

	double udrtemp=0.0;
	#pragma omp parallel for default(shared) reduction(+:udrtemp)
	for (int i=0;i<nsites;i++)
	{
		double x1=coords[i][0];
		double y1=coords[i][1];
		double z1=coords[i][2];
		
		double sx1=configx[i];
		double sy1=configy[i];
		double sz1=configz[i];
		for (int j=0;j<nsites;j++)
		{
			if (j!=i)
			{
				double x2=coords[j][0];
				double y2=coords[j][1];
				double z2=coords[j][2];
				double b=bmat(i,j);
				double c=cmat(i,j);
				double sx2=configx[j];
				double sy2=configy[j];
				double sz2=configz[j];
				double sisj=((sx1*sx2)+(sy1*sy2)+(sz1*sz2))*b;
				double rijx=rijxmat(i,j);	
				double rijy=rijymat(i,j);	
				double rijz=rijzmat(i,j);
				double sirsjr=((sx1*rijx+ sy1*rijy + sz1*rijz) * (sx2*rijx + sy2*rijy + sz2*rijz))*c;
				udrtemp+=(0.5*(sisj - sirsjr));
			}
		}
	}
	udr=udrtemp;
	///////////////////////////////////////////////////////////////
	// k space part of the energy
	//////////////////////////////////////////////////////////////
	double udktemp=0.0;
	double uself2=0.0;
	double uselffactor=1.0;
	#pragma omp parallel for default(shared) reduction(+:udktemp)
	for (int kint=0;kint<numkvecs;kint++)
	{
		double kx=kvecs[kint][0];	
		double ky=kvecs[kint][1];	
		double kz=kvecs[kint][2];	
		complex<double> sum=0.0;
		double sumr=0.0;
	
		/*		
		//Doulbe sum method - sum over i ! = j
		//OR sum method - sum over i  = j but then account for self energy
		// Are the 2 equivalent?? They should be, but looks like some error?
		for (int i=0;i<nsites;i++)
		{
			double x1=coords[i][0];
			double y1=coords[i][1];
			double z1=coords[i][2];
			double sx1=configx[i];
			double sy1=configy[i];
			double sz1=configz[i];
			double sidotk=(sx1*kx)+(sy1*ky)+(sz1*kz);
			for (int j=0;j<nsites;j++)
			{
				double x2=coords[j][0];
				double y2=coords[j][1];
				double z2=coords[j][2];
				double sx2=configx[j];
				double sy2=configy[j];
				double sz2=configz[j];
				double sjdotk=(sx2*kx)+(sy2*ky)+(sz2*kz);
				double kdotr=(kx*(x2-x1))+(ky*(y2-y1))+(kz*(z2-z1));
				//if (j!=i)
				{
					//sumr+=real(expmikr(kint,i)*conj(expmikr(kint,j)))*sidotk*sjdotk;
					sumr+=(cos(kdotr)*sidotk*sjdotk);
				}

			}
		}
		udk+=(sumr*kfactors[kint]);*/

		
		for (int i=0;i<nsites;i++)
		{
			double sx1=configx[i];
			double sy1=configy[i];
			double sz1=configz[i];
			double sidotk=(sx1*kx)+(sy1*ky)+(sz1*kz);
			sum+=((sidotk)*expmikr(kint,i));
		}
		
		double sum2=0.0;
		/*//cout<<"Groups size = "<<groups.size()<<endl;
		for (int g=0;g<groups.size();g++)
		{
			complex<double> sum3=0.0;
			//cout<<"Groups g . size() ="<<groups[g].size()<<endl;
			for (int i=0;i<groups[g].size();i++) // Sum for a single group
			{
				int site1=groups[g][i];
				double sx1=configx[site1];
				double sy1=configy[site1];
				double sz1=configz[site1];
				double sidotk=(sx1*kx)+(sy1*ky)+(sz1*kz);
				sum3+=(expmikr(kint,site1)*sidotk);
			}
			sum2+=(abs(sum3)*abs(sum3));
		}*/	
		dmoments[kint]=sum;
		udktemp+=( ((abs(sum)*abs(sum)) - sum2)*kfactors[kint]);
	}
	udk=udktemp;

	udk=udk*twopibv; //cubic
	//udk=udk*twopibv*64.0; // diamond
	///////////////////////////////////////////////////////////////
	// Surface part of the energy
	//////////////////////////////////////////////////////////////
	usurf=0.0;
	double mx=0.0;
	double my=0.0;
	double mz=0.0;

	for (int i=0;i<nsites;i++)
	{
		double sx=configx[i];
		double sy=configy[i];
		double sz=configz[i];
		mx+=sx;
		my+=sy;
		mz+=sz;
	}
	magnetization[0]=mx;
	magnetization[1]=my;
	magnetization[2]=mz;
	//usurf=(mx*mx)+(my*my)+(mz*mz);
	//usurf=usurf*twopibv/3.0; // ELiminate Usurf for speherical sample

	//////////////////////////////////////////////////////////////
	//Uself - self energy part of the energy
	//////////////////////////////////////////////////////////////

	uself=0.0;
	for (int i=0;i<nsites;i++)
	{
		double sx=configx[i];
		double sy=configy[i];
		double sz=configz[i];
		uself+=((sx*sx)+(sy*sy)+(sz*sz));
	}
	uself=uself*pow(alpha,3.0)*(-2.0/3.0)*invrootpi; /*Accounted for or not depends on summation idea for k space ... */

	/*cout<<"udr   = "<<udr<<endl;
	cout<<"udk   = "<<udk<<endl;
	cout<<"usurf = "<<usurf<<endl;
	cout<<"uself = "<<uself<<endl;*/
	total_e=udr+udk+usurf+uself;      
	total_e=total_e*coupling;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void difference_dipolar_energy(double 				coupling,
			    int 				L,
			    int 				site,
			    double 				x,
			    double 				y,
			    double 				z,
			    double				valx,
			    double				valy,
			    double				valz,
			    std::vector<double> 		&configx, 
			    std::vector<double> 		&configy, 
			    std::vector<double> 		&configz, 
		      	    std::vector< std::vector<double> >  &coords, 
		      	    RMatrix 				&rijxmat, 
		      	    RMatrix 				&rijymat, 
		      	    RMatrix 				&rijzmat, 
		      	    RMatrix 				&bmat, 
		      	    RMatrix 				&cmat, 
		      	    std::vector< std::vector<double> >  &kvecs,
		      	    Matrix 				&expmikr,
			    std::vector<double> 		&kfactors,
			    double 				alpha,
			    double 				&udrdiff,
			    double 				&udkdiff,
			    double 				&usurfdiff,
			    double 				&total_e_diff,
			    std::vector<complex<double> >	&dmoments,
			    std::vector<double>			&magnetization,
			    std::vector<complex<double> >	&dmomentsnew,
			    std::vector<double>			&magnetizationnew)
{
	// coupling > 0 and is AFM 
	double invrootpi=0.56418958354;
	int numkvecs=kvecs.size();
	//cout<<"Number of k vectors is = "<<numkvecs<<endl;
	int nsites=coords.size();
	double twopibv = (2.0*3.14159265359)/double(L*L*L);
	///////////////////////////////////////////////////////////////
	// Real space part of the energy
	//////////////////////////////////////////////////////////////
	// rij corresponds to the closest image (MINIMUM image convention)
	// All other images contribute nothing (zero) to the real space sum
	// alpha is chosen so that this condition holds

	double udrdifftemp=0.0;
	#pragma omp parallel for default(shared) reduction (+ : udrdifftemp)
	for (int j=0;j<nsites;j++)
	{
		if (j!=site)
		{
			double b=bmat(site,j);
			double c=cmat(site,j);
			double sx2=configx[j];
			double sy2=configy[j];
			double sz2=configz[j];
			double sisjdiff=(((valx-x)*sx2)+((valy-y)*sy2)+((valz-z)*sz2))*b;
			double rijx=rijxmat(site,j);	
			double rijy=rijymat(site,j);	
			double rijz=rijzmat(site,j);
			double sirsjrdiff=(((valx-x)*rijx+ (valy-y)*rijy + (valz-z)*rijz) * (sx2*rijx + sy2*rijy + sz2*rijz))*c;
			udrdifftemp+=((sisjdiff - sirsjrdiff));
		}
	}
	udrdiff=udrdifftemp;
	///////////////////////////////////////////////////////////////
	// k space part of the energy
	//////////////////////////////////////////////////////////////
	double udkdifftemp=0.0;
	#pragma omp parallel for default(shared) reduction (+ : udkdifftemp)
	for (int kint=0;kint<numkvecs;kint++)
	{
		double kx=kvecs[kint][0];	
		double ky=kvecs[kint][1];	
		double kz=kvecs[kint][2];	
		double sidotkdiff=((valx-x)*kx)+((valy-y)*ky)+((valz-z)*kz);
		dmomentsnew[kint]=dmoments[kint]+((sidotkdiff)*expmikr(kint,site));
		udkdifftemp+=(((abs(dmomentsnew[kint])*abs(dmomentsnew[kint])) -(abs(dmoments[kint])*abs(dmoments[kint]))) *kfactors[kint]);
	}
	udkdiff=udkdifftemp;

	udkdiff=udkdiff*twopibv;
	///////////////////////////////////////////////////////////////
	// Difference of the Surface part of the energy
	//////////////////////////////////////////////////////////////
	usurfdiff=0.0;
	double mx=magnetization[0];
	double my=magnetization[1];
	double mz=magnetization[2];
	double mxnew=mx+valx-x;
	double mynew=my+valy-y;
	double mznew=mz+valz-z;
	
	magnetizationnew[0]=mxnew;
	magnetizationnew[1]=mynew;
	magnetizationnew[2]=mznew;
	
	usurfdiff=((mxnew*mxnew)+(mynew*mynew)+(mznew*mznew)-(mx*mx)-(my*my)-(mz*mz));
	usurfdiff=usurfdiff*twopibv/3.0;
	
	//////////////////////////////////////////////////////////////
	//Uself - self energy part of the energy stays the same
	//////////////////////////////////////////////////////////////
	/*cout<<"udr   = "<<udr<<endl;
	cout<<"udk   = "<<udk<<endl;
	cout<<"usurf = "<<usurf<<endl;
	cout<<"uself = "<<uself<<endl;*/
	total_e_diff=udrdiff+udkdiff+usurfdiff;      
	total_e_diff=total_e_diff*coupling;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> total_magnetization(std::vector<double> &configx, 
					std::vector<double> &configy, 
					std::vector<double> &configz)
{
	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	for (int i=0;i<configx.size();i++)
	{
		mx+=configx[i];
		my+=configy[i];
		mz+=configz[i];
	}
	m.push_back(mx);
	m.push_back(my);
	m.push_back(mz);
	return m;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_111_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz)
{	
	configx.resize(nsites);
	configy.resize(nsites);
	configz.resize(nsites);
	for (int i=0;i<nsites;i++) 
	{
		double a=1.0;
		double b=1.0;
		double c=1.0;
		double norm=sqrt((a*a)+(b*b)+(c*c));
		configx[i]=a/norm;
		configy[i]=b/norm;
		configz[i]=c/norm;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_x_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz)
{	
	configx.resize(nsites);
	configy.resize(nsites);
	configz.resize(nsites);
	for (int i=0;i<nsites;i++) 
	{
		double a=1.0;
		double b=0.0;
		double c=0.0;
		double norm=sqrt((a*a)+(b*b)+(c*c));
		configx[i]=a/norm;
		configy[i]=b/norm;
		configz[i]=c/norm;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_random_config(int nsites, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz)
{	
	configx.resize(nsites);
	configy.resize(nsites);
	configz.resize(nsites);
	for (int i=0;i<nsites;i++) 
	{
		double a=2*uniform_rnd() - 1;
		double b=2*uniform_rnd() - 1;
		double c=2*uniform_rnd() - 1;
		double norm=sqrt((a*a)+(b*b)+(c*c));
		configx[i]=a/norm;
		configy[i]=b/norm;
		configz[i]=c/norm;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_neighbor2(int Lx,int Ly, int Lz, string direction, int x, int y, int z, int &xnew, int &ynew, int &znew, int &c)
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
void mc_finite_D(int L, string lattice, int nsamples, int nburn, string start_config, 
		string mcmove, double temperature, double hx, double hy, double hz, 
		double Jval, double D, double Dp, double alphaL, int kcut,
		double & eavg, 
		double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		double &mx2avg, double &my2avg, double &mz2avg)
{
	double alpha=alphaL/double(L);
	double mevbtesla=5.7883818012*(1.e-5)*(1.e3);
	double moment=5.7/2.0; // From Suzy/ Expt - the factor of 2.0 was added by HJC to half the coupling 
			       // to the magnetic field as spins are counted twice
	hx=hx*mevbtesla*moment;
	hy=hy*mevbtesla*moment;
	hz=hz*mevbtesla*moment;

	nburn=nsamples/10;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	temperature=temperature*0.08617; // Convert to meV
	double beta_final=1.0/temperature;
	double beta_initial=0.0;
	// LxLxL pbc lattice
	double beta_steps=(beta_final-beta_initial)/nburn;
	int nsites=L*L*L;
	if (lattice=="cubic") nsites=L*L*L;
	if (lattice=="diamond") nsites=8*L*L*L; // view it as a lattice with 8 sublattices
	RMatrix 				dmat(nsites,nsites); 
	RMatrix 				rijxmat(nsites,nsites); 
	RMatrix 				rijymat(nsites,nsites); 
	RMatrix 				rijzmat(nsites,nsites); 
	std::vector< std::vector<int> > 	neighbors(nsites);	
	std::vector< std::vector<int> > 	groups;	
	std::vector<std::vector<double> > 	coords;
	std::vector<std::vector<double> > 	kvecs;
	Matrix 					expmikr;
	std::vector<double> 			kfactors;
	std::vector<int> 			n1vec;
	std::vector<int> 			n2vec;
	std::vector<int> 			n3vec;
	std::vector<complex<double> > 		dmoments;
	std::vector<complex<double> > 		dmomentsnew;
	
	if (lattice=="cubic")
	{
		//////////////////////////////////////////////
		// Cubic lattice ........
		//////////////////////////////////////////////
		// Make sites on box of dimensions L,L,L
		int site=0;
		for (int n1=0;n1<L;n1++)
		{
			double x=double(n1);
			for (int n2=0;n2<L;n2++)
			{
				double y=double(n2);
				for (int n3=0;n3<L;n3++)
				{
					double z=double(n3);	
					std::vector<double> coord;
					coord.push_back(x);
					coord.push_back(y);
					coord.push_back(z);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord);
					std::vector<int> group;
					group.push_back(site);
					groups.push_back(group);
					site+=1;
				}
			}	
		}
		//////////////////////////////////////////////
		// Cubic lattice k vectors
		//////////////////////////////////////////////
		// Make sites on box of dimensions L,L,L
		double twopibL=(2.0*3.14159265359)/double(L);
		for (int n1=-kcut*L;n1<=kcut*L;n1++)
		{
			double kx=double(n1)*twopibL;
			for (int n2=-kcut*L;n2<=kcut*L;n2++)
			{
				double ky=double(n2)*twopibL;
				for (int n3=-kcut*L;n3<=kcut*L;n3++)
				{
					double kz=double(n3)*twopibL;
					if ((n1==0 and n2==0 and n3==0) or ( (n1*n1) + (n2*n2) + (n3*n3)> double(kcut*kcut*L*L)+1.0e-6) )
					{
						//cout<<"Disregarding k vec = 0,0,0 or k outside sphere"<<endl;
					}
					else
					{	
						std::vector<double> kvec;
						kvec.push_back(kx);
						kvec.push_back(ky);
						kvec.push_back(kz);
						kvecs.push_back(kvec);
					}
				}
			}	
		}
		int numkvecs=kvecs.size();
		expmikr.resize(numkvecs,nsites);
		kfactors.resize(numkvecs);
		dmoments.resize(numkvecs);
		dmomentsnew.resize(numkvecs);
		for (int k=0;k<numkvecs;k++)
		{
			double kx=kvecs[k][0];
			double ky=kvecs[k][1];
			double kz=kvecs[k][2];
			double k2=((kx*kx) + (ky*ky) + (kz*kz));
			kfactors[k]=exp(-k2/(4.0*alpha*alpha))/k2;
			for (int i=0;i<nsites;i++)
			{
				double x=coords[i][0];
				double y=coords[i][1];
				double z=coords[i][2];

				double kdotr=(kx*x + ky*y + kz*z);
				expmikr(k,i)=complex<double>(cos(kdotr),-sin(kdotr));
			}
		}
		//////////////////////////////////////////////
		// Make distance matrix ........
		//////////////////////////////////////////////
		cout<<"Number of sites recorded is = "<<site<<endl;
		for (int i=0;i<nsites;i++)
		{
			double x1=coords[i][0];
			double y1=coords[i][1];
			double z1=coords[i][2];
			for (int j=0;j<nsites;j++)
			{
				double x2=coords[j][0];
				double y2=coords[j][1];
				double z2=coords[j][2];
				
				double rijx=x2-x1;	
				double rijy=y2-y1;	
				double rijz=z2-z1;
				
				if (abs(rijx)>L/2.0 and x2>x1) rijx=(x2-x1)-L; // Translate x	
				if (abs(rijx)>L/2.0 and x2<x1) rijx=(x2-x1)+L; // Translate x 
				if (abs(rijy)>L/2.0 and y2>y1) rijy=(y2-y1)-L; // Translate y	
				if (abs(rijy)>L/2.0 and y2<y1) rijy=(y2-y1)+L; // Translate y
				if (abs(rijz)>L/2.0 and z2>z1) rijz=(z2-z1)-L; // Translate z	
				if (abs(rijz)>L/2.0 and z2<z1) rijz=(z2-z1)+L; // Translate z

				rijxmat(i,j)=rijx;
				rijymat(i,j)=rijy;
				rijzmat(i,j)=rijz;

				dmat(i,j)=sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));  // min distance between points in pbc
				//cout<<
				//cout<<"x1,y1,z1 ="<<x1<<","<<y1<<","<<z1<<"  --> x2,y2,z2 =  "<<x2<<","<<y2<<","<<z2<<" = "<<dmat(i,j)<<endl;
				if (dmat(i,j)>0.99 and dmat(i,j)<1.01) neighbors[i].push_back(j); // nearest neighbor on cube
			}
		}
		int six=0;
		for (int i=0;i<neighbors.size();i++)
		{
			cout<<" i = "<<i<<", nsize ="<<neighbors[i].size()<<endl;
			if (neighbors[i].size()==6) six+=1;
		}
		cout<<" Number of sites with 6 neighbors = "<<six<<endl;
	}

	if (lattice=="diamond")
	{
		//////////////////////////////////////////////
		// Diamond lattice ........
		//////////////////////////////////////////////
		// Make sites on box of dimensions L,L,L
		// Associate 8 sites with every n1,n2,n3
		int site=0;
		for (int n1=0;n1<L;n1++)
		{
			double x=double(n1);
			for (int n2=0;n2<L;n2++)
			{
				double y=double(n2);
				for (int n3=0;n3<L;n3++)
				{
					double z=double(n3);	
					std::vector<int> group;
					
					std::vector<double> coord1;
					coord1.push_back(x);
					coord1.push_back(y);
					coord1.push_back(z);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord1);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord2;
					coord2.push_back(x+0.5);
					coord2.push_back(y+0.5);
					coord2.push_back(z+0.0);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord2);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord3;
					coord3.push_back(x+0.5);
					coord3.push_back(y+0.0);
					coord3.push_back(z+0.5);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord3);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord4;
					coord4.push_back(x+0.0);
					coord4.push_back(y+0.5);
					coord4.push_back(z+0.5);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord4);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord5;
					coord5.push_back(x+0.25);
					coord5.push_back(y+0.25);
					coord5.push_back(z+0.25);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord5);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord6;
					coord6.push_back(x+0.75);
					coord6.push_back(y+0.75);
					coord6.push_back(z+0.25);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord6);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord7;
					coord7.push_back(x+0.25);
					coord7.push_back(y+0.75);
					coord7.push_back(z+0.75);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord7);
					group.push_back(site);
					site+=1;
					
					std::vector<double> coord8;
					coord8.push_back(x+0.75);
					coord8.push_back(y+0.25);
					coord8.push_back(z+0.75);
					n1vec.push_back(n1);	
					n2vec.push_back(n2);	
					n3vec.push_back(n3);	
					coords.push_back(coord8);
					group.push_back(site);
					site+=1;
					
					groups.push_back(group);
				}
			}	
		}
		//////////////////////////////////////////////
		// Cubic lattice k vectors used for diamond
		//////////////////////////////////////////////
		// Make sites on box of dimensions L,L,L
		double twopibL=(2.0*3.14159265359)/double(L);
		for (int n1=-kcut*L;n1<=kcut*L;n1++)
		{
			double kx=double(n1)*twopibL;
			for (int n2=-kcut*L;n2<=kcut*L;n2++)
			{
				double ky=double(n2)*twopibL;
				for (int n3=-kcut*L;n3<=kcut*L;n3++)
				{
					double kz=double(n3)*twopibL;
					if ((n1==0 and n2==0 and n3==0) or ( (n1*n1) + (n2*n2) + (n3*n3)> double(kcut*kcut*L*L)+1.0e-6) )
					{
						//cout<<"Disregarding k vec = 0,0,0 or k outside sphere"<<endl;
					}
					else
					{	
						std::vector<double> kvec;
						kvec.push_back(kx);
						kvec.push_back(ky);
						kvec.push_back(kz);
						kvecs.push_back(kvec);
					}
				}
			}	
		}
		int numkvecs=kvecs.size();
		expmikr.resize(numkvecs,nsites);
		kfactors.resize(numkvecs);
		dmoments.resize(numkvecs);
		dmomentsnew.resize(numkvecs);
		for (int k=0;k<numkvecs;k++)
		{
			double kx=kvecs[k][0];
			double ky=kvecs[k][1];
			double kz=kvecs[k][2];
			double k2=(kx*kx + ky*ky + kz*kz);
			kfactors[k]=exp(-k2/(4.0*alpha*alpha))/k2;
			for (int i=0;i<nsites;i++)
			{
				double x=coords[i][0];
				double y=coords[i][1];
				double z=coords[i][2];

				double kdotr=(kx*x + ky*y + kz*z);
				expmikr(k,i)=complex<double>(cos(kdotr),-sin(kdotr));
			}
		}
		//////////////////////////////////////////////
		// Make distance matrix ........
		//////////////////////////////////////////////
		cout<<"Number of sites recorded is = "<<site<<endl;
		for (int i=0;i<nsites;i++)
		{
			double x1=coords[i][0];
			double y1=coords[i][1];
			double z1=coords[i][2];
			
			int n1x=n1vec[i];
			int n1y=n2vec[i];
			int n1z=n3vec[i];
			for (int j=0;j<nsites;j++)
			{
				int n2x=n1vec[j];
				int n2y=n2vec[j];
				int n2z=n3vec[j];
				
				double x2=coords[j][0];
				double y2=coords[j][1];
				double z2=coords[j][2];
				
				int nijx=n2x-n1x;	
				int nijy=n2y-n1y;	
				int nijz=n2z-n1z;
				
				double rijx=x2-x1;	
				double rijy=y2-y1;	
				double rijz=z2-z1;
				
				if (abs(nijx)>L/2.0 and n2x>n1x) rijx=(x2-x1)-L; // Translate x	
				if (abs(nijx)>L/2.0 and n2x<n1x) rijx=(x2-x1)+L; // Translate x 
				if (abs(nijy)>L/2.0 and n2y>n1y) rijy=(y2-y1)-L; // Translate y	
				if (abs(nijy)>L/2.0 and n2y<n1y) rijy=(y2-y1)+L; // Translate y
				if (abs(nijz)>L/2.0 and n2z>n1z) rijz=(z2-z1)-L; // Translate z	
				if (abs(nijz)>L/2.0 and n2z<n1z) rijz=(z2-z1)+L; // Translate z

				rijxmat(i,j)=rijx;
				rijymat(i,j)=rijy;
				rijzmat(i,j)=rijz;

				dmat(i,j)=sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));  // min distance between points in pbc
				//cout<<
				//cout<<"x1,y1,z1 ="<<x1<<","<<y1<<","<<z1<<"  --> x2,y2,z2 =  "<<x2<<","<<y2<<","<<z2<<" = "<<dmat(i,j)<<endl;
				if (dmat(i,j)>0.43 and dmat(i,j)<0.44) neighbors[i].push_back(j); // nearest neighbor on diamond
			}
		}
		int four=0;
		for (int i=0;i<neighbors.size();i++)
		{
			cout<<" i = "<<i<<", nsize ="<<neighbors[i].size()<<endl;
			if (neighbors[i].size()==4) four+=1;
		}
		cout<<" Number of sites with 4 neighbors = "<<four<<endl;
	}
	
	cout<<endl;
	// Make random configuration 
	std::vector<double> configx,configy,configz;	
	if (start_config=="random") make_random_config(nsites,configx,configy,configz); 
	if (start_config=="111")    make_111_config(nsites,configx,configy,configz); 
	if (start_config=="x")      make_x_config(nsites,configx,configy,configz); 
	RMatrix bmat(nsites,nsites),cmat(nsites,nsites);
	//double alpha;
	/////////////////////////////////////////////////////////////////////////////
        // Exchange (Nearest neighbor)
	///////////////////////////////////////////////////////////////////////////
	double jenergy=total_j_energy(Jval,configx,configy,configz,neighbors);
	
	/*alpha=5.0/double(L);
	#make_bmat_cmat(dmat,alpha,bmat,cmat);
	#double dipolar1=total_dipolar_energy(L,configx,configy,configz,coords,dmat,bmat,cmat,kvecs,alpha);
	
	alpha=10.0/double(L);
	make_bmat_cmat(dmat,alpha,bmat,cmat);
	double dipolar2=total_dipolar_energy(L,configx,configy,configz,coords,dmat,bmat,cmat,kvecs,alpha);

	alpha=15.0/double(L);
	make_bmat_cmat(dmat,alpha,bmat,cmat);
	double dipolar3=total_dipolar_energy(L,configx,configy,configz,coords,dmat,bmat,cmat,kvecs,alpha);
	
	
	alpha=20.0/double(L);
	make_bmat_cmat(dmat,alpha,bmat,cmat);
	double dipolar4=total_dipolar_energy(L,configx,configy,configz,coords,dmat,bmat,cmat,kvecs,alpha);*/
	
	/*cout<<"dipolar 1 ="<<dipolar1<<endl;
	cout<<"dipolar 2 ="<<dipolar2<<endl;
	cout<<"dipolar 3 ="<<dipolar3<<endl;
	cout<<"dipolar 4 ="<<dipolar4<<endl;*/

	/////////////////////////////////////////////////////////////////////////////
        // Dipolar (Long range Ewald)
	///////////////////////////////////////////////////////////////////////////
	
	make_bmat_cmat(dmat,alpha,bmat,cmat,n1vec,n2vec,n3vec);
	std::vector<double> magnetization(3);
	std::vector<double> magnetizationnew(3);
	double udr,udk,usurf,uself,dipolar;
	
	total_dipolar_energy(Dp, L,configx,configy,configz,coords, groups, rijxmat,rijymat,rijzmat,bmat,cmat,kvecs,expmikr, kfactors, alpha, udr, udk, uself, usurf, dipolar, dmoments, magnetization);
		
	/////////////////////////////////////////////////////////////////////////////
        // Magnetic field coupling (on site)
	///////////////////////////////////////////////////////////////////////////
	double henergy=-(hx*magnetization[0])-(hy*magnetization[1])-(hz*magnetization[2]); // Mag energy
		
	/////////////////////////////////////////////////////////////////////////////
        // D (cubic anisotropy) energy (on site)
	///////////////////////////////////////////////////////////////////////////
	double denergy=0.0;
	for (int i=0;i<nsites;i++) denergy=denergy-D*pow(configx[i],4.0)-D*pow(configy[i],4.0)-D*pow(configz[i],4.0);

	cout<<"Total j (exchange)       energy ="<<jenergy<<endl;
	cout<<"Total dipolar            energy ="<<dipolar<<endl;
	cout<<"Total magnetic field     energy ="<<henergy<<endl;
	cout<<"Total D (anisotropy)     energy ="<<denergy<<endl;

	// Total energy
	double energy=jenergy+dipolar+henergy+denergy;
	cout<<"Total j energy + dipolar + total h + D energy ="<<energy<<endl;
	

	// magnetization 
	double mx=magnetization[0];
	double my=magnetization[1];
	double mz=magnetization[2];
	double etot=0.0;
	double e2tot=0.0;
	double e4tot=0.0;
	double mxtot=0.0;
	double mytot=0.0;
	double mztot=0.0;
	double mx2tot=0.0;
	double my2tot=0.0;
	double mz2tot=0.0;
	double mx4tot=0.0;
	double my4tot=0.0;
	double mz4tot=0.0;
	double accept=0.0;
	double reject=0.0;
        double nmeas=0;
	// Accept reject Metropolis
	for (int i=0; i<(nsamples+nburn);i++)
	{
		if (i%100000==0) cout<<"Move number "<<i<<endl;
		// Choose random site
		int site=uniform_rand_int(0,nsites);
		double valx,valy,valz;
		// Initialize the new configuration
		double x=configx[site];		
		double y=configy[site];		
		double z=configz[site];		
	
		// Choose a completely random direction - This is INEFFICENT at low temps
		if (mcmove=="random")  random_move_continuous_spin(valx,valy,valz);
		// Instead we want to choose something CLOSE to the direction we have
		// It should not be too close either 
		if (mcmove=="conical") conical_move_continuous_spin(0.3,x,y,z,valx,valy,valz);
		if (mcmove=="infDspecial")    infD_move_special_continuous_spin(x,y,z,valx,valy,valz);
		if (mcmove=="largeD")
		{
			double temp=uniform_rnd();
			if (temp>0.5)
			{
				conical_move_continuous_spin(0.3,x,y,z,valx,valy,valz);
			}
			else
			{
				big_move_continuous_spin(x,y,z,valx,valy,valz);
			}

		}	
		//if (i%2==0)  random_move_continuous_spin(valx,valy,valz);
		//random_move_continuous_spin(valx,valy,valz);
		// Instead we want to choose something CLOSE to the direction we have
		// It should not be too close either 
		//else conical_move_continuous_spin(0.3,x,y,z,valx,valy,valz);
		//conical_move_continuous_spin(0.3,x,y,z,valx,valy,valz);


		// Normalize new direction
		double norm=sqrt(valx*valx + valy*valy + valz*valz); 
		valx=valx/norm;
		valy=valy/norm;
		valz=valz/norm;

		// Calculate local J energy of old and new configs
		//double local_jenergy1=local_j_energy(site,x,y,z,configx,configy,configz,neighbors);
		//double local_jenergy2=local_j_energy(site,valx,valy,valz,configx,configy,configz,neighbors);
		double jen_diff=j_energy_diff(Jval,site,x,y,z,valx,valy,valz,configx,configy,configz,neighbors);
		
		// Calculate dipolar energy difference
		//Temporary change
		configx[site]=valx;		
		configy[site]=valy;		
		configz[site]=valz;	

		double udrnew,udknew,usurfnew,uselfnew;
		////// will be sped up....
		double dipolarnew;
		/*total_dipolar_energy(Dp, L,configx,configy,configz,coords,groups,rijxmat,rijymat,rijzmat,bmat,cmat,kvecs,expmikr,kfactors,alpha,udrnew,udknew,usurfnew,uselfnew, dipolarnew, dmomentsnew,magnetizationnew);
		double dipen_diff1=dipolarnew-dipolar;	*/
		double dipen_diff=0.0;
		if (abs(Dp)>1.0e-6)
		{
		difference_dipolar_energy(Dp, L, site, x,y,z, valx, valy, valz, configx,configy,configz,coords,rijxmat,rijymat,rijzmat,bmat,cmat,kvecs,expmikr,kfactors,alpha,udrnew,udknew,usurfnew, dipen_diff, dmoments, magnetization, dmomentsnew,magnetizationnew);
		}
		//cout<<"Dipolar energy difference method 1 = "<<dipen_diff1<<endl;
		//cout<<"Dipolar energy difference method 2 = "<<dipen_diff<<endl;
		// Calculate difference of magnetization 
		dipolarnew=dipolar+dipen_diff;
		double mxdiff=valx-x;
		double mydiff=valy-y;
		double mzdiff=valz-z;
		double mx4diff=pow(valx,4)-pow(x,4);
		double my4diff=pow(valy,4)-pow(y,4);
		double mz4diff=pow(valz,4)-pow(z,4);
		double mnewx=mx+mxdiff;
		double mnewy=my+mydiff;
		double mnewz=mz+mzdiff;
				// h on site                         // D on site              // local j energy  // Dipole energy diff
		double ediff=(-hx*mxdiff)+(-hy*mydiff)+(-hz*mzdiff)-D*(mx4diff+my4diff+mz4diff)+(jen_diff)+(dipen_diff);
		double beta;
		if (start_config=="random") 
		{
			beta=min(beta_final,i*beta_steps); // Gradually cool from infinite temp if start config is random
		}
		else       // if ordered phase then dont cool but let the temperature play its role
		{
			beta=beta_final;                  
		}
		//cout<<"Beta = "<<beta<<endl;
		double prob=exp(-beta*ediff);
		double rand=uniform_rnd();
		if (rand<prob) // Accept
		{
			accept=accept+1;
			// New values
			configx[site]=valx;
			configy[site]=valy;
			configz[site]=valz;
			double enew=energy+ediff;
						     // 1 sweep = Nsite moves
			if (i>nburn and i%nsites==0) // Update averages - after 1 sweep 
			{
				etot=etot+enew;
				e2tot=e2tot+pow(enew,2.0);
				e4tot=e4tot+pow(enew,4.0);
				mxtot=mxtot+mnewx;
				mytot=mytot+mnewy;
				mztot=mztot+mnewz;
				mx2tot=mx2tot+pow(mnewx,2.0);
				my2tot=my2tot+pow(mnewy,2.0);
				mz2tot=mz2tot+pow(mnewz,2.0);
				mx4tot=mx4tot+pow(mnewx,4.0);
				my4tot=my4tot+pow(mnewy,4.0);
				mz4tot=mz4tot+pow(mnewz,4.0);
		       	        nmeas=nmeas+1;
			}
			energy=enew;  // Update for the next move
			mx=mnewx;
			my=mnewy;
			mz=mnewz;
			dipolar=dipolarnew;
			magnetization=magnetizationnew;
			#pragma omp parallel for
			for (int kint=0;kint<dmoments.size();kint++) dmoments[kint]=dmomentsnew[kint];
		}
		else // reject and Update averages
		{
		       reject=reject+1;
		       if (i>nburn and i%nsites==0)
		       {
			       etot=etot+(energy);
			       e2tot=e2tot+pow((energy),2.0);
			       e4tot=e4tot+pow((energy),4.0);
			       mx2tot=mx2tot+pow(mx,2.0);
			       my2tot=my2tot+pow(my,2.0);
			       mz2tot=mz2tot+pow(mz,2.0);
			       mx4tot=mx4tot+pow(mx,4.0);
			       my4tot=my4tot+pow(my,4.0);
			       mz4tot=mz4tot+pow(mz,4.0);
			       mxtot=mxtot+mx;
			       mytot=mytot+my;
			       mztot=mztot+mz;
			       nmeas=nmeas+1;
		       }
		       configx[site]=x; // Old values		
		       configy[site]=y;		
		       configz[site]=z;		
		}
	}
	accept=accept/double(nsamples+nburn);
	eavg=etot/double(nmeas);
	e2avg=e2tot/double(nmeas);
	e4avg=e4tot/double(nmeas);
	mxavg=mxtot/double(nmeas);
	myavg=mytot/double(nmeas);
	mzavg=mztot/double(nmeas);
	mx2avg=mx2tot/double(nmeas);
	my2avg=my2tot/double(nmeas);
	mz2avg=mz2tot/double(nmeas);
	mx4avg=mx4tot/double(nmeas);
	my4avg=my4tot/double(nmeas);
	mz4avg=mz4tot/double(nmeas);
	double spheat=(e2avg-(eavg*eavg))/(temperature*temperature);
	double spheatpersite=spheat/double(nsites);
       
	cout<<"accept = "<<boost::format("%+ .15f") %accept<<endl;
	cout<<"mxavg  = "<<boost::format("%+ .15f") %mxavg<<endl;
	cout<<"myavg  = "<<boost::format("%+ .15f") %myavg<<endl;
	cout<<"mzavg  = "<<boost::format("%+ .15f") %mzavg<<endl;
	cout<<"mx2avg = "<<boost::format("%+ .15f") %mx2avg<<endl;
	cout<<"my2avg = "<<boost::format("%+ .15f") %my2avg<<endl;
	cout<<"mz2avg = "<<boost::format("%+ .15f") %mz2avg<<endl;
	cout<<"mx4avg = "<<boost::format("%+ .15f") %mx4avg<<endl;
	cout<<"my4avg = "<<boost::format("%+ .15f") %my4avg<<endl;
	cout<<"mz4avg = "<<boost::format("%+ .15f") %mz4avg<<endl;
	cout<<"eavg   = "<<boost::format("%+ .15f") %eavg<<endl;
	cout<<"e2avg  = "<<boost::format("%+ .15f") %e2avg<<endl;
	cout<<"e4avg  = "<<boost::format("%+ .15f") %e4avg<<endl;
	cout<<"Cv     = "<<boost::format("%+ .15f") %spheat<<endl;
	cout<<"Cvps   = "<<boost::format("%+ .15f") %spheatpersite<<endl;
}
