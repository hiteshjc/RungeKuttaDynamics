#ifndef STENSOR_HEADER
#define STENSOR_HEADER

#include"global.h"

class STensor
{
	public:
	int 				     numqs;
	int 				     numomegas;
	std::vector< std::vector<double> >   qvals;
	std::vector< std::vector<complex<double> >  >	sxx, sxy, sxz, syx, syy, syz, szx, szy, szz, sperp;
	int 				     nmeas;
	std::vector<double> 		     omegas;
 
	void init(int numqs, int numomegas)
	{
		this->qvals.clear();
		this->nmeas=1;
		this->numqs=numqs;
		this->numomegas=numomegas; 
		this->sxx.clear(); this->sxy.clear(); this->sxz.clear();
		this->syx.clear(); this->syy.clear(); this->syz.clear();
		this->szx.clear(); this->szy.clear(); this->szz.clear();
		this->sperp.clear();
	
		for (int om=0; om<this->numomegas; om++)
		{
			this->sxx.push_back(std::vector<complex<double> > ());
			this->sxy.push_back(std::vector<complex<double> > ());
			this->sxz.push_back(std::vector<complex<double> > ());
			
			this->syx.push_back(std::vector<complex<double> > ());
			this->syy.push_back(std::vector<complex<double> > ());
			this->syz.push_back(std::vector<complex<double> > ());
			
			this->szx.push_back(std::vector<complex<double> > ());
			this->szy.push_back(std::vector<complex<double> > ());
			this->szz.push_back(std::vector<complex<double> > ());

			this->sperp.push_back(std::vector<complex<double> > ());
			
			this->sxx[om].resize(this->numqs);this->sxy[om].resize(this->numqs);this->sxz[om].resize(this->numqs);
			this->syx[om].resize(this->numqs);this->syy[om].resize(this->numqs);this->syz[om].resize(this->numqs);
			this->szx[om].resize(this->numqs);this->szy[om].resize(this->numqs);this->szz[om].resize(this->numqs);
			this->sperp[om].resize(this->numqs);

			#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++)
			{
				this->sxx[om][nq]=0.0; this->sxy[om][nq]=0.0;this->sxz[om][nq]=0.0;
				this->syx[om][nq]=0.0; this->syy[om][nq]=0.0;this->syz[om][nq]=0.0;
				this->szx[om][nq]=0.0; this->szy[om][nq]=0.0;this->szz[om][nq]=0.0;
				this->sperp[om][nq]=0.0;
			}
		}
	}
	
	void update_totals(STensor &stemp)
	{
		this->nmeas+=1;
		for (int om=0; om<this->numomegas; om++)
		{
			#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++)
			{
				this->sxx[om][nq]+=stemp.sxx[om][nq];
				this->sxy[om][nq]+=stemp.sxy[om][nq];
				this->sxz[om][nq]+=stemp.sxz[om][nq];
				
				this->syx[om][nq]+=stemp.syx[om][nq];
				this->syy[om][nq]+=stemp.syy[om][nq];
				this->syz[om][nq]+=stemp.syz[om][nq];
				
				this->szx[om][nq]+=stemp.szx[om][nq];
				this->szy[om][nq]+=stemp.szy[om][nq];
				this->szz[om][nq]+=stemp.szz[om][nq];
				
				this->sperp[om][nq]+=stemp.sperp[om][nq];
			}
		}
	}
	
	void copy(STensor &stemp)
	{
		this->nmeas=stemp.nmeas;
		this->qvals=stemp.qvals;
		for (int om=0; om<this->numomegas; om++)
		{
			#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++)
			{
				this->sxx[om][nq]=stemp.sxx[om][nq];
				this->sxy[om][nq]=stemp.sxy[om][nq];
				this->sxz[om][nq]=stemp.sxz[om][nq];
				
				this->syx[om][nq]=stemp.syx[om][nq];
				this->syy[om][nq]=stemp.syy[om][nq];
				this->syz[om][nq]=stemp.syz[om][nq];
				
				this->szx[om][nq]=stemp.szx[om][nq];
				this->szy[om][nq]=stemp.szy[om][nq];
				this->szz[om][nq]=stemp.szz[om][nq];
				
				this->sperp[om][nq]=stemp.sperp[om][nq];
			}
		}
	}

	void average()
	{
		for (int om=0; om<this->numomegas; om++)
		{
			//#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++) 
			{
				this->sxx[om][nq]=this->sxx[om][nq]/double(this->nmeas);
				this->sxy[om][nq]=this->sxy[om][nq]/double(this->nmeas);
				this->sxz[om][nq]=this->sxz[om][nq]/double(this->nmeas);
			
				this->syx[om][nq]=this->syx[om][nq]/double(this->nmeas);
				this->syy[om][nq]=this->syy[om][nq]/double(this->nmeas);
				this->syz[om][nq]=this->syz[om][nq]/double(this->nmeas);
			
				this->szx[om][nq]=this->szx[om][nq]/double(this->nmeas);
				this->szy[om][nq]=this->szy[om][nq]/double(this->nmeas);
				this->szz[om][nq]=this->szz[om][nq]/double(this->nmeas);
				
				this->sperp[om][nq]=this->sperp[om][nq]/double(this->nmeas);
			}
		}
	}
	
};

#endif
