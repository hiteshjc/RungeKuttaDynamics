#include"mc_pyrochlore.h"
#include"mc_finite_D.h"
#include"printing_functions.h"
#include"matrix_functions.h"
#include"number_functions.h"
#include"stensor.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_j_tensor( RMatrix Jmat01, RMatrix Jmat10,
		  RMatrix Jmat02, RMatrix Jmat20,
		  RMatrix Jmat03, RMatrix Jmat30,
		  RMatrix Jmat12, RMatrix Jmat21,
		  RMatrix Jmat13, RMatrix Jmat31,
		  RMatrix Jmat23, RMatrix Jmat32,	
		  std::vector<RMatrix> &Jtensor)
{
	RMatrix empty(0,0);
	Jtensor.clear();
	Jtensor.push_back(empty);
	Jtensor.push_back(Jmat01);
	Jtensor.push_back(Jmat02);
	Jtensor.push_back(Jmat03);
	Jtensor.push_back(Jmat10);
	Jtensor.push_back(empty);
	Jtensor.push_back(Jmat12);
	Jtensor.push_back(Jmat13);
	Jtensor.push_back(Jmat20);
	Jtensor.push_back(Jmat21);
	Jtensor.push_back(empty);
	Jtensor.push_back(Jmat23);
	Jtensor.push_back(Jmat30);
	Jtensor.push_back(Jmat31);
	Jtensor.push_back(Jmat32);
	Jtensor.push_back(empty);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void local_fields(double &spin, int &site, int &t,
		      std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors,
		      std::vector< std::vector<int> > &nneighbors,
		      std::vector<RMatrix> &Jtensor, 
		      double &Jnnn, RMatrix &bond_disorder_matrix,
		      std::vector< std::vector<int> > &ijkt,
		      double &eff_field_x,
		      double &eff_field_y,
		      double &eff_field_z)
{
	eff_field_x=0.0; eff_field_y=0.0; eff_field_z=0.0;
	RMatrix Jmat;
	// Disorder free nearest neighbor terms
	for (int j=0;j<neighbors[site].size();j++)
	{
		int k=neighbors[site][j];
		int t2=ijkt[k][3];
		int cind=(4*t)+t2;
		eff_field_x+=(Jtensor[cind](0,0)*configx[k])+(Jtensor[cind](0,1)*configy[k])+(Jtensor[cind](0,2)*configz[k]);
		eff_field_y+=(Jtensor[cind](1,0)*configx[k])+(Jtensor[cind](1,1)*configy[k])+(Jtensor[cind](1,2)*configz[k]);
		eff_field_z+=(Jtensor[cind](2,0)*configx[k])+(Jtensor[cind](2,1)*configy[k])+(Jtensor[cind](2,2)*configz[k]);
	}
	// Bond disorder nearest neighbor Heisenberg
	for (int j=0;j<neighbors[site].size();j++)
	{
		int k=neighbors[site][j];
		double J=bond_disorder_matrix(site,k);
		eff_field_x+=(J*configx[k]);
		eff_field_y+=(J*configy[k]);
		eff_field_z+=(J*configz[k]);
	}
	// Next nn Heisenberg
	for (int j=0;j<nneighbors[site].size();j++)
	{
		int k=nneighbors[site][j];
		eff_field_x+=(Jnnn*configx[k]);
		eff_field_y+=(Jnnn*configy[k]);
		eff_field_z+=(Jnnn*configz[k]);
	}
	eff_field_x=eff_field_x*spin; eff_field_y=eff_field_y*spin; eff_field_z=eff_field_z*spin;
}
			      
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void cross_product(double spin, double sx, double sy, double sz, double eff_field_x, double eff_field_y, double eff_field_z, 
		   double &crossx, double &crossy, double &crossz)
{

	crossx=sy*eff_field_z - sz*eff_field_y;
	crossy=sz*eff_field_x - sx*eff_field_z;
	crossz=sx*eff_field_y - sy*eff_field_x;

	crossx=crossx*spin; crossy=crossy*spin; crossz=crossz*spin;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Time evolution of initial configuration
// This initial configuration can be anything - but we will typically choose something from a thermal distribution
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void time_evolve(   double spin, double deltat,
		    double tottime, int L, 
		    std::vector<double> &configx, 
		    std::vector<double> &configy, 
		    std::vector<double> &configz, 
		    std::vector< std::vector<int> > &neighbors, 
		    std::vector< std::vector<int> > &nneighbors, 
		    RMatrix &Jmat01, RMatrix &Jmat10,
		    RMatrix &Jmat02, RMatrix &Jmat20,
		    RMatrix &Jmat03, RMatrix &Jmat30,
		    RMatrix &Jmat12, RMatrix &Jmat21,
		    RMatrix &Jmat13, RMatrix &Jmat31,
		    RMatrix &Jmat23, RMatrix &Jmat32,
		    double &Jnnn, RMatrix &bond_disorder_matrix,
		    std::vector< std::vector<double> > &fullcoords,
		    std::vector< std::vector<int> >   & ijkt, 
		    std::vector<std::vector<double> > & qvals,
		    Matrix 			      & phases,
		    std::vector<double> 	      & omegas, 
		    STensor 			      & smunu)
{
	int numqs=qvals.size();
	int nsites=configx.size();
	int numomegas=omegas.size();
      
	smunu.init(numqs,numomegas);
      	smunu.qvals=qvals;
      	smunu.omegas=omegas;

	std::vector< std::vector<complex<double> > > sxtot, sytot, sztot;
	std::vector<complex<double> > 		     sxtq,  sytq,  sztq;

	sxtq.resize(numqs);sytq.resize(numqs);sztq.resize(numqs);	
	for (int om=0; om<numomegas; om++)
	{
		sxtot.push_back(std::vector< complex<double> >());
		sytot.push_back(std::vector< complex<double> >());
		sztot.push_back(std::vector< complex<double> >());
		sxtot[om].resize(numqs);sytot[om].resize(numqs);sztot[om].resize(numqs);	
		#pragma omp parallel for
		for (int nq=0;nq<numqs;nq++) 
		{
			sxtot[om][nq]=0.0;sytot[om][nq]=0.0;sztot[om][nq]=0.0;
		}
	}

	// 4th order Runge-Kutta intermediate vectors
	std::vector<double> config0x(nsites),config0y(nsites),config0z(nsites);
	std::vector<double> config1x(nsites),config1y(nsites),config1z(nsites);
	std::vector<double> config2x(nsites),config2y(nsites),config2z(nsites);
	std::vector<double> config3x(nsites),config3y(nsites),config3z(nsites);
	
	std::vector<double> k1x(nsites),k1y(nsites),k1z(nsites);
	std::vector<double> k2x(nsites),k2y(nsites),k2z(nsites);
	std::vector<double> k3x(nsites),k3y(nsites),k3z(nsites);
	std::vector<double> k4x(nsites),k4y(nsites),k4z(nsites);
	
	// Initialize configuration to configuration provided 
	config0x=configx; config0y=configy; config0z=configz;
	
	// How many time steps should one evolve for
	int numtimes=int(tottime/deltat);

	std::vector<RMatrix> Jtensor; 
	make_j_tensor(Jmat01, Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32, Jtensor);

	// Loop over time steps
	tottime=double(numtimes)*deltat;
	for (int intt=0;intt<numtimes;intt++)
	{
		double energyj=total_j_energy(spin, config0x,config0y,config0z,
							 neighbors,nneighbors,
							 Jmat01, Jmat10, 
							 Jmat02, Jmat20, 
							 Jmat03, Jmat30, 
							 Jmat12, Jmat21, 
							 Jmat13, Jmat31, 
							 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
	
		double time=double(intt)*deltat;
		//cout<<"Time = "<<boost::format("%+ .5f") %time<<" S(0)S(t) = "<<boost::format("%+ .10f") %((config0x[0]*configx[0])+(config0y[0]*configy[0])+(config0z[0]*configz[0]))<<"  E = "<<boost::format("%+ .10f") %energyj<<endl;
		// Print configuration of spins along with coordinates
		if (intt==0)
		{
			for (int i=0;i<nsites;i++)
			{
				double x=fullcoords[i][0]; double y=fullcoords[i][1]; double z=fullcoords[i][2];
				double sx=config0x[i]; double sy=config0y[i]; double sz=config0z[i];
				double norm=(sx*sx) + (sy*sy) + (sz*sz);
				cout<<boost::format("%+ .10f  %+ .10f  %+ .10f  %+ .10f   %+ .10f   %+ .10f  %+ .10f") %x %y %z %sx %sy %sz %norm<<endl;
			}
		}
		cout<<"Time = "<<boost::format("%+ .5f") %time<<" S(0)S(t) = "<<boost::format("%+ .10f") %((config0x[0]*configx[0])+(config0y[0]*configy[0])+(config0z[0]*configz[0]))<<"  E = "<<boost::format("%+ .10f") %energyj<<endl;
		//cout<<endl;

		fourier_transforms(spin,phases,config0x,config0y,config0z,sxtq,sytq,sztq);
		
		for (int om=0; om<numomegas; om++)
		{
			double omegat=omegas[om]*time;
			complex<double> timephase=complex<double>(cos(omegat),sin(omegat))*(deltat/tottime);	
			#pragma omp parallel for
			for (int nq=0;nq<numqs;nq++)
			{
				sxtot[om][nq]+=(sxtq[nq]*timephase);
				sytot[om][nq]+=(sytq[nq]*timephase);
				sztot[om][nq]+=(sztq[nq]*timephase);
			}
	        }

		///////////////////////////////////////////////////////////
		// 1 st step of Runge Kutta - get k1
		double sign=-1.0;

		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      int t=ijkt[i][3];
			      double crossx,crossy,crossz;
			      double eff_field_x,eff_field_y,eff_field_z;
			      double sx_i=config0x[i]; double sy_i=config0y[i]; double sz_i=config0z[i];
			      local_fields(spin,i,t,config0x,config0y,config0z,neighbors, nneighbors, Jtensor,
					     Jnnn, bond_disorder_matrix, ijkt, eff_field_x, eff_field_y, eff_field_z);
			      cross_product(spin, sign*sx_i, sign*sy_i, sign*sz_i, eff_field_x, eff_field_y, eff_field_z, crossx, crossy, crossz);
			      k1x[i]=crossx;
			      k1y[i]=crossy;
			      k1z[i]=crossz;
		}

		///////////////////////////////////////////////////////////
		// 2 nd step of Runge Kutta - get k2
		// yn + hk1/2 configuration
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      config1x[i]=config0x[i]+(deltat*0.5*k1x[i]); 
			      config1y[i]=config0y[i]+(deltat*0.5*k1y[i]); 
			      config1z[i]=config0z[i]+(deltat*0.5*k1z[i]); 
		}
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      int t=ijkt[i][3];
			      double crossx,crossy,crossz;
			      double eff_field_x,eff_field_y,eff_field_z;
			      double sx_i=config1x[i]; double sy_i=config1y[i]; double sz_i=config1z[i];
			      local_fields(spin,i,t, config1x,config1y,config1z,neighbors, nneighbors, Jtensor, 
					     Jnnn, bond_disorder_matrix, ijkt, eff_field_x, eff_field_y, eff_field_z);
			      cross_product(spin, sign*sx_i, sign*sy_i, sign*sz_i, eff_field_x, eff_field_y, eff_field_z, crossx, crossy, crossz);
			      k2x[i]=crossx;
			      k2y[i]=crossy;
			      k2z[i]=crossz;
		}

		///////////////////////////////////////////////////////////
		// 3 rd step of Runge Kutta - get k3
		// yn+hk2/2 configuration
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      config2x[i]=config0x[i]+(deltat*0.5*k2x[i]); 
			      config2y[i]=config0y[i]+(deltat*0.5*k2y[i]); 
			      config2z[i]=config0z[i]+(deltat*0.5*k2z[i]); 
		}
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      int t=ijkt[i][3];
			      double crossx,crossy,crossz;
			      double eff_field_x,eff_field_y,eff_field_z;
			      double sx_i=config2x[i]; double sy_i=config2y[i]; double sz_i=config2z[i];
			      local_fields(spin,i,t, config2x,config2y,config2z,neighbors, nneighbors, Jtensor, 
					     Jnnn, bond_disorder_matrix, ijkt, eff_field_x, eff_field_y, eff_field_z);
			      cross_product(spin, sign*sx_i, sign*sy_i, sign*sz_i, eff_field_x, eff_field_y, eff_field_z, crossx, crossy, crossz);
			      k3x[i]=crossx;
			      k3y[i]=crossy;
			      k3z[i]=crossz;
		}
		
		///////////////////////////////////////////////////////////
		// 4th step of Runge Kutta - get k4
		// yn + h k3
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      config3x[i]=config0x[i]+(deltat*1.0*k3x[i]); 
			      config3y[i]=config0y[i]+(deltat*1.0*k3y[i]); 
			      config3z[i]=config0z[i]+(deltat*1.0*k3z[i]); 
		}
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			      int t=ijkt[i][3];
			      double crossx,crossy,crossz;
			      double eff_field_x,eff_field_y,eff_field_z;
			      double sx_i=config3x[i]; double sy_i=config3y[i]; double sz_i=config3z[i];
			      local_fields(spin,i,t, config3x,config3y,config3z,neighbors, nneighbors, Jtensor, 
					     Jnnn, bond_disorder_matrix, ijkt, eff_field_x, eff_field_y, eff_field_z);
			      cross_product(spin, sign*sx_i, sign*sy_i, sign*sz_i, eff_field_x, eff_field_y, eff_field_z, crossx, crossy, crossz);
			      k4x[i]=crossx;
			      k4y[i]=crossy;
			      k4z[i]=crossz;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// NEW Configuration for the next step given by Runge Kutta formula
		// y_n+1 = y_n + h/6 (k1 + 2 k2 + 2 k3 + k4 )
		#pragma omp parallel for
		for (int i=0;i<nsites;i++)
		{
			config0x[i]+=((deltat/6.0)*(k1x[i]+(2.0*k2x[i])+(2.0*k3x[i])+k4x[i]));
			config0y[i]+=((deltat/6.0)*(k1y[i]+(2.0*k2y[i])+(2.0*k3y[i])+k4y[i]));
			config0z[i]+=((deltat/6.0)*(k1z[i]+(2.0*k2z[i])+(2.0*k3z[i])+k4z[i]));
		}
      }

      for (int om=0;om<numomegas;om++)
      {	
	      #pragma omp parallel for
	      for (int nq=0;nq<numqs;nq++) 
	      {
			smunu.sxx[om][nq]=sxtot[om][nq]*conj(sxtot[om][nq])/double(nsites);
			smunu.sxy[om][nq]=sxtot[om][nq]*conj(sytot[om][nq])/double(nsites);
			smunu.sxz[om][nq]=sxtot[om][nq]*conj(sztot[om][nq])/double(nsites);

			smunu.syx[om][nq]=sytot[om][nq]*conj(sxtot[om][nq])/double(nsites);
			smunu.syy[om][nq]=sytot[om][nq]*conj(sytot[om][nq])/double(nsites);
			smunu.syz[om][nq]=sytot[om][nq]*conj(sztot[om][nq])/double(nsites);

			smunu.szx[om][nq]=sztot[om][nq]*conj(sxtot[om][nq])/double(nsites);
			smunu.szy[om][nq]=sztot[om][nq]*conj(sytot[om][nq])/double(nsites);
			smunu.szz[om][nq]=sztot[om][nq]*conj(sztot[om][nq])/double(nsites);
	      }
      }

      /*cout<<"======================================================================================================================================="<<endl;
      cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)          Sperp(Q)     "<<endl;
      cout<<"======================================================================================================================================="<<endl;*/
      for (int om=0; om<numomegas; om++)
      {
      	      #pragma omp parallel for
	      for (int i=0;i<numqs;i++)
	      {
		double qx=qvals[i][0]; double qy=qvals[i][1]; double qz=qvals[i][2];
		complex<double> sxx=smunu.sxx[om][i]; complex<double> syy=smunu.syy[om][i]; complex<double> szz=smunu.szz[om][i];
		complex<double> sxy=smunu.sxy[om][i]; complex<double> syx=smunu.syx[om][i];
		complex<double> sxz=smunu.sxz[om][i]; complex<double> szx=smunu.szx[om][i];
		complex<double> syz=smunu.syz[om][i]; complex<double> szy=smunu.szy[om][i];

		double q2=((qx*qx)+(qy*qy)+(qz*qz));
		complex<double> sperp=0.0;
		if (abs(q2)>1.0e-6)
		{
			sperp=(1.0 - (qx*qx/q2))*sxx + (1.0 - (qy*qy/q2))*syy + (1.0 - (qz*qz/q2))*szz + (-sxy*qx*qy/q2) + (-syx*qy*qx/q2) + (-sxz*qx*qz/q2) + (-szx*qx*qz/q2) + (-syz*qz*qy/q2) + (-szy*qz*qy/q2); 
		}
		else
		{
			sperp=sxx + syy + szz; 
		}
		//cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f  %+ .8f  %+ .8f %+ .8f  %+ .8f") %qx %qy %qz %sxx %sxy %sxz %syx %syy %syz %szx %szy %szz %sperp<<endl;
	       smunu.sperp[om][i]=sperp; 
	     }
     }
}
