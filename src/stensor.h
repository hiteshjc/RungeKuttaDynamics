#ifndef STENSOR_HEADER
#define STENSOR_HEADER

#include"global.h"

class STensor
{
	public:
	int 				     numqs;
	std::vector< std::vector<double> >   qvals;
	std::vector<complex<double> > 	     sxx;
	std::vector<complex<double> > 	     sxy;
	std::vector<complex<double> > 	     sxz;
	std::vector<complex<double> > 	     syx;
	std::vector<complex<double> > 	     syy;
	std::vector<complex<double> > 	     syz;
	std::vector<complex<double> > 	     szx;
	std::vector<complex<double> > 	     szy;
	std::vector<complex<double> > 	     szz;
	std::vector<complex<double> > 	     sperp;
	int 				     nmeas;
 
	void init(int numqs)
	{
		this->qvals.clear();
		this->nmeas=1;
		this->numqs=numqs;
		this->sxx.resize(this->numqs);this->sxy.resize(this->numqs);this->sxz.resize(this->numqs);
		this->syx.resize(this->numqs);this->syy.resize(this->numqs);this->syz.resize(this->numqs);
		this->szx.resize(this->numqs);this->szy.resize(this->numqs);this->szz.resize(this->numqs);
		this->sperp.resize(this->numqs);

		#pragma omp parallel for
		for (int nq=0;nq<this->numqs;nq++)
		{
			this->sxx[nq]=0.0; this->sxy[nq]=0.0;this->sxz[nq]=0.0;
			this->syx[nq]=0.0; this->syy[nq]=0.0;this->syz[nq]=0.0;
			this->szx[nq]=0.0; this->szy[nq]=0.0;this->szz[nq]=0.0;
			this->sperp[nq]=0.0;
		}
	}
	
	void update_totals(STensor &stemp)
	{
		this->nmeas+=1;
		#pragma omp parallel for
		for (int nq=0;nq<this->numqs;nq++)
		{
			this->sxx[nq]+=stemp.sxx[nq];
			this->sxy[nq]+=stemp.sxy[nq];
			this->sxz[nq]+=stemp.sxz[nq];
			
			this->syx[nq]+=stemp.syx[nq];
			this->syy[nq]+=stemp.syy[nq];
			this->syz[nq]+=stemp.syz[nq];
			
			this->szx[nq]+=stemp.szx[nq];
			this->szy[nq]+=stemp.szy[nq];
			this->szz[nq]+=stemp.szz[nq];
			
			this->sperp[nq]+=stemp.sperp[nq];
		}
	}
	
	void copy(STensor &stemp)
	{
		this->nmeas=stemp.nmeas;
		this->qvals=stemp.qvals;
		#pragma omp parallel for
		for (int nq=0;nq<this->numqs;nq++)
		{
			this->sxx[nq]=stemp.sxx[nq];
			this->sxy[nq]=stemp.sxy[nq];
			this->sxz[nq]=stemp.sxz[nq];
			
			this->syx[nq]=stemp.syx[nq];
			this->syy[nq]=stemp.syy[nq];
			this->syz[nq]=stemp.syz[nq];
			
			this->szx[nq]=stemp.szx[nq];
			this->szy[nq]=stemp.szy[nq];
			this->szz[nq]=stemp.szz[nq];
			
			this->sperp[nq]=stemp.sperp[nq];
		}
	}

	void average()
	{
		//#pragma omp parallel for
		for (int nq=0;nq<this->numqs;nq++) 
		{
			this->sxx[nq]=this->sxx[nq]/double(this->nmeas);
			this->sxy[nq]=this->sxy[nq]/double(this->nmeas);
			this->sxz[nq]=this->sxz[nq]/double(this->nmeas);
		
			this->syx[nq]=this->syx[nq]/double(this->nmeas);
			this->syy[nq]=this->syy[nq]/double(this->nmeas);
			this->syz[nq]=this->syz[nq]/double(this->nmeas);
		
			this->szx[nq]=this->szx[nq]/double(this->nmeas);
			this->szy[nq]=this->szy[nq]/double(this->nmeas);
			this->szz[nq]=this->szz[nq]/double(this->nmeas);
			
			this->sperp[nq]=this->sperp[nq]/double(this->nmeas);
		}
	}
	
};

#endif
