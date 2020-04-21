#include"mc_pyrochlore.h"
#include"mc_finite_D.h"
#include"printing_functions.h"
#include"matrix_functions.h"
#include"number_functions.h"
#include"runge_kutta.h"
//#include"spin_wave.h"
#include"stensor.h"

using namespace std;

class QMC_Info
{
	public:
	int 			nsites;
	int 			numqs;
	std::vector< std::vector<double> >    qvals;
	std::vector<double> 	configx, configy, configz;
	std::vector<complex<double> > sxq;
	std::vector<complex<double> > syq;
	std::vector<complex<double> > szq;
	std::vector<complex<double> > sxsxtot;
	std::vector<complex<double> > sxsytot;
	std::vector<complex<double> > sxsztot;
	std::vector<complex<double> > sysxtot;
	std::vector<complex<double> > sysytot;
	std::vector<complex<double> > sysztot;
	std::vector<complex<double> > szsxtot;
	std::vector<complex<double> > szsytot;
	std::vector<complex<double> > szsztot;
	double 			etot,e2tot,e4tot;
	double 			mxtot,mx2tot,mx4tot;
	double 			mytot,my2tot,my4tot;
	double 			mztot,mz2tot,mz4tot;
	double 			eavg,e2avg,e4avg;
	double 			mxavg,mx2avg,mx4avg;
	double 			myavg,my2avg,my4avg;
	double 			mzavg,mz2avg,mz4avg;
	double 			nmeas;
	double 			spheatpersitenodim; 
	double 			spheat;
	double			spheatpersite;
	double			temp,tempKelvin,beta;
	double 			energy,mx,my,mz;
	Matrix 			phases;
	bool                    measure_corrs;
        double                  accept;
        double                  reject;
        double                  nswaps;
 
	void init(int L, int nsites, double tempKelvin, bool measure_corrs)
	{
		this->accept=0.0;
		this->reject=0.0;
		this->nswaps=0.0;
	
		this->measure_corrs=measure_corrs;
		this->nsites=nsites;
		this->tempKelvin=tempKelvin;
        	this->temp=tempKelvin*0.08621738;
		this->beta=1.0/temp;
		this->etot=0.0;this->e2tot=0.0;this->e4tot=0.0;
		this->mxtot=0.0;this->mx2tot=0.0;this->mx4tot=0.0;
		this->mytot=0.0;this->my2tot=0.0;this->my4tot=0.0;
		this->mztot=0.0;this->mz2tot=0.0;this->mz4tot=0.0;
		this->nmeas=0;

		std::vector<double>qs;
		
		if (this->measure_corrs)
		{
			for (int i=0;i<=4*L;i++)
			{
				double q=double(i)/double(L);
				qs.push_back(q);
			} 
			int onedq=qs.size();	
			this->numqs=onedq*onedq*onedq;

			for (int nqx=0;nqx<onedq;nqx++)
			{
				for (int nqy=0;nqy<onedq;nqy++)
				{
					for (int nqz=0;nqz<onedq;nqz++)
					{
						std::vector<double> qvec;
						qvec.push_back(qs[nqx]);
						qvec.push_back(qs[nqy]);
						qvec.push_back(qs[nqz]);
						this->qvals.push_back(qvec);
					}
				}
			}
			this->sxq.resize(this->numqs);
			this->syq.resize(this->numqs);
			this->szq.resize(this->numqs);
			this->sxsxtot.resize(this->numqs);
			this->sxsytot.resize(this->numqs);
			this->sxsztot.resize(this->numqs);
			this->sysxtot.resize(this->numqs);
			this->sysytot.resize(this->numqs);
			this->sysztot.resize(this->numqs);
			this->szsxtot.resize(this->numqs);
			this->szsytot.resize(this->numqs);
			this->szsztot.resize(this->numqs);
		}

		if (this->measure_corrs)
		{
			#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++)
			{
				this->sxq[nq]=0.0;
				this->syq[nq]=0.0;
				this->szq[nq]=0.0;

				this->sxsxtot[nq]=0.0;
				this->sxsytot[nq]=0.0;
				this->sxsztot[nq]=0.0;
				
				this->sysxtot[nq]=0.0;
				this->sysytot[nq]=0.0;
				this->sysztot[nq]=0.0;
				
				this->szsxtot[nq]=0.0;
				this->szsytot[nq]=0.0;
				this->szsztot[nq]=0.0;
			}
		}
	}
	
	void update_totals()
	{
		this->etot+=this->energy;this->e2tot+=pow(this->energy,2.0);this->e4tot+=pow(this->energy,4.0);
		this->mxtot+=this->mx;this->mytot+=this->my;this->mztot+=this->mz;
		this->mx2tot+=pow(this->mx,2.0);this->my2tot+=pow(this->my,2.0);this->mz2tot+=pow(this->mz,2.0);
		this->mx4tot+=pow(this->mx,4.0);this->my4tot+=pow(this->my,4.0);this->mz4tot+=pow(this->mz,4.0);
		this->nmeas+=1;
		
		if (this->measure_corrs)
		{
			#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++)
			{
				this->sxsxtot[nq]+=(this->sxq[nq]*conj(this->sxq[nq]));
				this->sxsytot[nq]+=(this->sxq[nq]*conj(this->syq[nq]));
				this->sxsztot[nq]+=(this->sxq[nq]*conj(this->szq[nq]));
				
				this->sysxtot[nq]+=(this->syq[nq]*conj(this->sxq[nq]));
				this->sysytot[nq]+=(this->syq[nq]*conj(this->syq[nq]));
				this->sysztot[nq]+=(this->syq[nq]*conj(this->szq[nq]));
				
				this->szsxtot[nq]+=(this->szq[nq]*conj(this->sxq[nq]));
				this->szsytot[nq]+=(this->szq[nq]*conj(this->syq[nq]));
				this->szsztot[nq]+=(this->szq[nq]*conj(this->szq[nq]));
			}
		}
	}

	void average()
	{
		//cout<<"Nmeas = "<<this->nmeas<<endl;
		//cout<<"Numqs = "<<this->numqs<<endl;
		this->eavg=this->etot/this->nmeas;this->e2avg=this->e2tot/this->nmeas;this->e4avg=this->e4tot/this->nmeas;
		this->mxavg=this->mxtot/this->nmeas;this->myavg=this->mytot/this->nmeas;this->mzavg=this->mztot/this->nmeas;
		this->mx2avg=this->mx2tot/this->nmeas;this->my2avg=this->my2tot/this->nmeas;this->mz2avg=this->mz2tot/this->nmeas;
		this->mx4avg=this->mx4tot/this->nmeas;this->my4avg=this->my4tot/this->nmeas;this->mz4avg=this->mz4tot/this->nmeas;
		this->spheatpersitenodim=(this->e2avg-(this->eavg*this->eavg))/(this->temp*this->temp*double(this->nsites));
		this->spheat=(1119.67107046)*(this->e2avg-(this->eavg*this->eavg))/(this->tempKelvin*this->tempKelvin);
		this->spheatpersite=this->spheat/double(this->nsites);
		
		if (this->measure_corrs)
		{
			//#pragma omp parallel for
			for (int nq=0;nq<this->numqs;nq++) 
			{
				this->sxsxtot[nq]=this->sxsxtot[nq]/double(this->nmeas*this->nsites);
				this->sxsytot[nq]=this->sxsytot[nq]/double(this->nmeas*this->nsites);
				this->sxsztot[nq]=this->sxsztot[nq]/double(this->nmeas*this->nsites);
			
				this->sysxtot[nq]=this->sysxtot[nq]/double(this->nmeas*this->nsites);
				this->sysytot[nq]=this->sysytot[nq]/double(this->nmeas*this->nsites);
				this->sysztot[nq]=this->sysztot[nq]/double(this->nmeas*this->nsites);
			
				this->szsxtot[nq]=this->szsxtot[nq]/double(this->nmeas*this->nsites);
				this->szsytot[nq]=this->szsytot[nq]/double(this->nmeas*this->nsites);
				this->szsztot[nq]=this->szsztot[nq]/double(this->nmeas*this->nsites);
			}
		}
	}
	
};

////////////////////////////////////////////////////////////////////////
void make_bond_disorder_matrix(double disorder_strength,
			       std::vector< std::vector<int> > &neighbors,
			       RMatrix &bond_disorder_matrix)
{
	int nsites=neighbors.size();
	#pragma omp parallel for
	for (int i=0;i<nsites*nsites;i++) bond_disorder_matrix[i]=0.0;

	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<neighbors[i].size();j++)
		{
			int site=neighbors[i][j];
			if (i<site)
			{
				bond_disorder_matrix(i,site)=(2.0*uniform_rnd()-1.0)*disorder_strength;
				bond_disorder_matrix(site,i)=bond_disorder_matrix(i,site);
			}
		}
	}
	
	cout<<"========================================================================================================================="<<endl;
	cout<<" i       j      Bij				                                                                        "<<endl;
	cout<<"========================================================================================================================="<<endl;
	for (int i=0;i<nsites;i++)
	{
		for (int j=0;j<nsites;j++) if (abs(bond_disorder_matrix(i,j))>1.0e-6) cout<<boost::format("%5d    %5d    %+ .5f") %i %j %bond_disorder_matrix(i,j)<<endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
void fourier_transforms(double spin,Matrix &phases, 
			std::vector<double> &configx,
			std::vector<double> &configy,
			std::vector<double> &configz,
			std::vector<complex<double> > &sxq,
			std::vector<complex<double> > &syq,
			std::vector<complex<double> > &szq)
{
	int numsites=phases.NCols();
	int numqs=phases.NRows();
	#pragma omp parallel for
	for (int i=0;i<numqs;i++) 
	{
		sxq[i]=0.0;
		syq[i]=0.0;
		szq[i]=0.0;
		for (int j=0;j<numsites;j++)
		{
			sxq[i]+=(configx[j]*phases(i,j)*spin);
			syq[i]+=(configy[j]*phases(i,j)*spin);
			szq[i]+=(configz[j]*phases(i,j)*spin);
		}
	}


}
////////////////////////////////////////////////////////////////////////////////
void fourier_transforms_slow(double spin,
			std::vector<std::vector<double> > &fullcoords, 
			std::vector<std::vector<double> > &qvals, 
			std::vector<double> &configx,
			std::vector<double> &configy,
			std::vector<double> &configz,
			std::vector<complex<double> > &sxq,
			std::vector<complex<double> > &syq,
			std::vector<complex<double> > &szq)
{
	double tpi=2.0*3.1415926;
	int numsites=fullcoords.size();
	int numqs=qvals.size();
	#pragma omp parallel for
	for (int i=0;i<numqs;i++) 
	{
		sxq[i]=0.0;
		syq[i]=0.0;
		szq[i]=0.0;
		for (int j=0;j<numsites;j++)
		{
			double qdotr=(qvals[i][0]*fullcoords[j][0]);
			       qdotr+=(qvals[i][1]*fullcoords[j][1]);
			       qdotr+=(qvals[i][2]*fullcoords[j][2]);
			complex<double> phase=exp(tpi*complex<double>(0,1)*qdotr);
			sxq[i]+=(configx[j]*phase*spin);
			syq[i]+=(configy[j]*phase*spin);
			szq[i]+=(configz[j]*phase*spin);
		}
	}


}


////////////////////////////////////////////////////////////////////////////////

void make_phases( std::vector< std::vector<double> > &fullcoords,
		  std::vector< std::vector<double> > &qvals,
		  Matrix &phases)
{
	double tpi=2.0*3.1415926;
	int nsites=fullcoords.size();
	int numqs=qvals.size();
	phases.resize(numqs,nsites);
	for (int i=0;i<numqs;i++)
	{
		for (int j=0;j<nsites;j++)
		{
			double qdotr=(qvals[i][0]*fullcoords[j][0]);
			       qdotr+=(qvals[i][1]*fullcoords[j][1]);
			       qdotr+=(qvals[i][2]*fullcoords[j][2]);
			phases(i,j)=exp(tpi*complex<double>(0,1)*qdotr);
		}

	}

}
////////////////////////////////////////////////////////////////////////
void make_qs(int L, std::vector<std::vector<double> > &qvals)
{	
	qvals.clear();
	std::vector<double> qs; 
	for (int i=0;i<=4*L;i++)
	{
		double q=double(i)/double(L);
		qs.push_back(q);
	} 
	int onedq=qs.size();	
	int numqs=onedq*onedq*onedq;

	for (int nqx=0;nqx<onedq;nqx++)
	{
		for (int nqy=0;nqy<onedq;nqy++)
		{
			for (int nqz=0;nqz<onedq;nqz++)
			{
				std::vector<double> qvec;
				qvec.push_back(qs[nqx]);
				qvec.push_back(qs[nqy]);
				qvec.push_back(qs[nqz]);
				qvals.push_back(qvec);
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////
void make_single_tet(int L,
		     std::vector< std::vector<double> > &fullcoords,
		     std::vector< std::vector<int> > &ijkt,
		     std::vector< std::vector<int> > &neighbors,
		     std::vector< std::vector<int> > &nneighbors,
		     bool measure_corrs)
{
	cout<<"Making pyrochlore..."<<endl;
	fullcoords.clear();
	ijkt.clear();
	//////////////////////////////////////////////
	//             Pyrochlore lattice
	//////////////////////////////////////////////
	// Make sites on cube of dimensions L,L,L
	// Associate 16 sites with every n1,n2,n3
	//
	// Tetrahedron coords (used in Lucile/Ross paper)
	// r0 = 1/8 (+1,+1,+1)
	// r1 = 1/8 (+1,-1,-1)
	// r2 = 1/8 (-1,+1,-1)
	// r3 = 1/8 (-1,-1,+1)

	double u=1.0/8.0;
	double r0x=+u;double r0y=+u;double r0z=+u;
	double r1x=+u;double r1y=-u;double r1z=-u;
	double r2x=-u;double r2y=+u;double r2z=-u;
	double r3x=-u;double r3y=-u;double r3z=+u;

	std::vector<int>    ijktentry;
	std::vector<double> fullcoordsentry;
	int nx=0;
	int ny=0;
	int nz=0;	
	int site=0;
	
	///////////////////////////////////////////////////////////////////////////////////////////////
	ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(0);	
	fullcoordsentry.push_back(r0x);fullcoordsentry.push_back(r0y);fullcoordsentry.push_back(r0z);
	ijkt.push_back(ijktentry);
	fullcoords.push_back(fullcoordsentry);
	site+=1;
	ijktentry.clear();
	fullcoordsentry.clear();
	///////////////////////////////////////////////////////////////////////////////////////////////

	ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(1);	
	fullcoordsentry.push_back(r1x);fullcoordsentry.push_back(r1y);fullcoordsentry.push_back(r1z);
	ijkt.push_back(ijktentry);
	fullcoords.push_back(fullcoordsentry);
	site+=1;
	ijktentry.clear();
	fullcoordsentry.clear();
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(2);	
	fullcoordsentry.push_back(r2x);fullcoordsentry.push_back(r2y);fullcoordsentry.push_back(r2z);
	ijkt.push_back(ijktentry);
	fullcoords.push_back(fullcoordsentry);
	site+=1;
	ijktentry.clear();
	fullcoordsentry.clear();
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(3);	
	fullcoordsentry.push_back(r3x);fullcoordsentry.push_back(r3y);fullcoordsentry.push_back(r3z);
	ijkt.push_back(ijktentry);
	fullcoords.push_back(fullcoordsentry);
	site+=1;
	ijktentry.clear();
	fullcoordsentry.clear();
	///////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"Number of sites recorded is = "<<site<<endl;
	int nsites=site;
	RMatrix dmat(nsites,nsites);
	for (int i=0;i<nsites;i++)
	{
		double x1=fullcoords[i][0];
		double y1=fullcoords[i][1];
		double z1=fullcoords[i][2];
		for (int j=0;j<nsites;j++)
		{
			double x2=fullcoords[j][0];
			double y2=fullcoords[j][1];
			double z2=fullcoords[j][2];
			
			double nijx=ijkt[j][0]-ijkt[i][0];	
			double nijy=ijkt[j][1]-ijkt[i][1];	
			double nijz=ijkt[j][2]-ijkt[i][2];
			
                        double rijx=x2-x1;	
			double rijy=y2-y1;	
			double rijz=z2-z1;
			
			/*if ((x2-x1)>=L) rijx=(x2-x1)-L; // Translate x	
			if ((x1-x2)>=L) rijx=(x2-x1)+L; // Translate x 
			if ((y2-y1)>=L) rijy=(y2-y1)-L; // Translate y	
			if ((y1-y2)>=L) rijy=(y2-y1)+L; // Translate y	
			if ((z2-z1)>=L) rijz=(z2-z1)-L; // Translate z	
			if ((z1-z2)>=L) rijz=(z2-z1)+L; // Translate z
			*/
			
			/*if (abs(nijx)>=L/2.0 and nijx>0) rijx=(x2-x1)-L; // Translate x	
			if (abs(nijx)>=L/2.0 and nijx<0) rijx=(x2-x1)+L; // Translate x 
			if (abs(nijy)>=L/2.0 and nijy>0) rijy=(y2-y1)-L; // Translate y	
			if (abs(nijy)>=L/2.0 and nijy<0) rijy=(y2-y1)+L; // Translate y
			if (abs(nijz)>=L/2.0 and nijz>0) rijz=(z2-z1)-L; // Translate z	
			if (abs(nijz)>=L/2.0 and nijz<0) rijz=(z2-z1)+L; // Translate z*/

			dmat(i,j)=sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));  // min distance between points in pbc
			//cout<<"dmat(i,j) = "<<dmat(i,j)<<endl;
			if (dmat(i,j)>0.353 and dmat(i,j)<0.354)  neighbors[i].push_back(j);  //      nearest neighbor on pyrochlore
			if (dmat(i,j)>0.612 and dmat(i,j)<0.613)  nneighbors[i].push_back(j); // next nearest neighbor on pyrochlore
		}
	}

	if (measure_corrs)
	{	
		int six1=0;	
		int twelve1=0;	
		for (int i=0;i<nsites;i++)
		{
			cout<<" i = "<<i<<", nsize,nnsize ="<<neighbors[i].size()<<"  "<<nneighbors[i].size()<<endl;
			if (neighbors[i].size()==6)   six1+=1;
			if (nneighbors[i].size()==12) twelve1+=1;
		}

		cout<<"Number of sites with  6      nearest neighbors = "<<six1<<endl;
		cout<<"Number of sites with 12 next nearest neighbors = "<<twelve1<<endl;

		/*cout<<"========================================================================================================================="<<endl;
		cout<<" i       j      dij(min dist between i and j in pbc)                                                                     "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<nsites;i++)
		{
			for (int j=0;j<nsites;j++) cout<<boost::format("%3d    %3d    %+ .5f") %i %j %dmat(i,j)<<endl;
		}*/
	}
}

////////////////////////////////////////////////////////////////////////
void make_pyrochlore(int L,
		     std::vector< std::vector<double> > &fullcoords,
		     std::vector< std::vector<int> > &ijkt,
		     std::vector< std::vector<int> > &neighbors,
		     std::vector< std::vector<int> > &nneighbors,
		     bool measure_corrs)
{
	cout<<"Making pyrochlore..."<<endl;
	fullcoords.clear();
	ijkt.clear();
	//////////////////////////////////////////////
	//             Pyrochlore lattice
	//////////////////////////////////////////////
	// Make sites on cube of dimensions L,L,L
	// Associate 16 sites with every n1,n2,n3
	//
	// Tetrahedron coords (used in Lucile/Ross paper)
	// r0 = 1/8 (+1,+1,+1)
	// r1 = 1/8 (+1,-1,-1)
	// r2 = 1/8 (-1,+1,-1)
	// r3 = 1/8 (-1,-1,+1)

	double u=1.0/8.0;
	double r0x=+u;double r0y=+u;double r0z=+u;
	double r1x=+u;double r1y=-u;double r1z=-u;
	double r2x=-u;double r2y=+u;double r2z=-u;
	double r3x=-u;double r3y=-u;double r3z=+u;

        // FCC coords....
	// (0.0,0.0,0.0)
	// (0.5,0.5,0.0)
	// (0.5,0.0,0.5)
	// (0.0,0.5,0.5)
        	
        std::vector< std::vector<double> > fcc_coords;
        std::vector<double> coord;
       
	coord.push_back(0.0);
        coord.push_back(0.0);
        coord.push_back(0.0);
	fcc_coords.push_back(coord);
	coord.clear();
        
	coord.push_back(0.5);
        coord.push_back(0.5);
        coord.push_back(0.0);
	fcc_coords.push_back(coord);
	coord.clear();

	coord.push_back(0.5);
        coord.push_back(0.0);
        coord.push_back(0.5);
	fcc_coords.push_back(coord);
	coord.clear();
	
	coord.push_back(0.0);
        coord.push_back(0.5);
        coord.push_back(0.5);
	fcc_coords.push_back(coord);
	
        int site=0;
	for (int nx=0;nx<L;nx++)
	{
		double x=double(nx);
		for (int ny=0;ny<L;ny++)
		{
			double y=double(ny);
			for (int nz=0;nz<L;nz++)
			{
				double z=double(nz);
				for (int p=0; p<fcc_coords.size();p++)
				{
					double xb=fcc_coords[p][0];	
					double yb=fcc_coords[p][1];	
					double zb=fcc_coords[p][2];
	
					std::vector<int>    ijktentry;
					std::vector<double> fullcoordsentry;
					
					///////////////////////////////////////////////////////////////////////////////////////////////
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(0);	
					fullcoordsentry.push_back(x+xb+r0x);fullcoordsentry.push_back(y+yb+r0y);fullcoordsentry.push_back(z+zb+r0z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
	
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(1);	
					fullcoordsentry.push_back(x+xb+r1x);fullcoordsentry.push_back(y+yb+r1y);fullcoordsentry.push_back(z+zb+r1z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
					
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(2);	
					fullcoordsentry.push_back(x+xb+r2x);fullcoordsentry.push_back(y+yb+r2y);fullcoordsentry.push_back(z+zb+r2z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
					
					ijktentry.push_back(nx);ijktentry.push_back(ny);ijktentry.push_back(nz); ijktentry.push_back(3);	
					fullcoordsentry.push_back(x+xb+r3x);fullcoordsentry.push_back(y+yb+r3y);fullcoordsentry.push_back(z+zb+r3z);
					ijkt.push_back(ijktentry);
					fullcoords.push_back(fullcoordsentry);
					site+=1;
					ijktentry.clear();
					fullcoordsentry.clear();
					///////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}
	
	cout<<"Number of sites recorded is = "<<site<<endl;
	int nsites=site;
	RMatrix dmat(nsites,nsites);
	for (int i=0;i<nsites;i++)
	{
		double x1=fullcoords[i][0];
		double y1=fullcoords[i][1];
		double z1=fullcoords[i][2];
		for (int j=0;j<nsites;j++)
		{
			double x2=fullcoords[j][0];
			double y2=fullcoords[j][1];
			double z2=fullcoords[j][2];
			
			double nijx=ijkt[j][0]-ijkt[i][0];	
			double nijy=ijkt[j][1]-ijkt[i][1];	
			double nijz=ijkt[j][2]-ijkt[i][2];
			
                        double rijx=x2-x1;	
			double rijy=y2-y1;	
			double rijz=z2-z1;
	
			if (L==1)
			{
				rijx=min(min(abs((x2-x1)-L),abs((x2-x1)+L)),abs(x2-x1));		
				rijy=min(min(abs((y2-y1)-L),abs((y2-y1)+L)),abs(y2-y1));		
				rijz=min(min(abs((z2-z1)-L),abs((z2-z1)+L)),abs(z2-z1));		
				//if ((x2-x1)>=L) rijx=(x2-x1)-L; // Translate x	
				//if ((x1-x2)>=L) rijx=(x2-x1)+L; // Translate x 
				//if ((y2-y1)>=L) rijy=(y2-y1)-L; // Translate y	
				//if ((y1-y2)>=L) rijy=(y2-y1)+L; // Translate y	
				//if ((z2-z1)>=L) rijz=(z2-z1)-L; // Translate z	
				//if ((z1-z2)>=L) rijz=(z2-z1)+L; // Translate z
			}
			else
			{	
				rijx=min(min(abs((x2-x1)-L),abs((x2-x1)+L)),abs(x2-x1));		
				rijy=min(min(abs((y2-y1)-L),abs((y2-y1)+L)),abs(y2-y1));		
				rijz=min(min(abs((z2-z1)-L),abs((z2-z1)+L)),abs(z2-z1));		
				/*if (abs(nijx)>=L/2.0 and nijx>0) rijx=(x2-x1)-L; // Translate x	
				if (abs(nijx)>=L/2.0 and nijx<0) rijx=(x2-x1)+L; // Translate x 
				if (abs(nijy)>=L/2.0 and nijy>0) rijy=(y2-y1)-L; // Translate y	
				if (abs(nijy)>=L/2.0 and nijy<0) rijy=(y2-y1)+L; // Translate y
				if (abs(nijz)>=L/2.0 and nijz>0) rijz=(z2-z1)-L; // Translate z	
				if (abs(nijz)>=L/2.0 and nijz<0) rijz=(z2-z1)+L; // Translate z*/
			}
			dmat(i,j)=sqrt((rijx*rijx)+(rijy*rijy)+(rijz*rijz));  // min distance between points in pbc
			//cout<<"dmat(i,j) = "<<dmat(i,j)<<endl;
			if (dmat(i,j)>0.353 and dmat(i,j)<0.354)  neighbors[i].push_back(j);  //      nearest neighbor on pyrochlore
			if (dmat(i,j)>0.612 and dmat(i,j)<0.613)  nneighbors[i].push_back(j); // next nearest neighbor on pyrochlore
		}
	}

	//if (measure_corrs)
	{	
		int six1=0;	
		int twelve1=0;	
		for (int i=0;i<nsites;i++)
		{
			cout<<" i = "<<i<<", nsize,nnsize ="<<neighbors[i].size()<<"  "<<nneighbors[i].size()<<endl;
			if (neighbors[i].size()==6)   six1+=1;
			if (nneighbors[i].size()==12) twelve1+=1;
		}

		cout<<"Number of sites with  6      nearest neighbors = "<<six1<<endl;
		cout<<"Number of sites with 12 next nearest neighbors = "<<twelve1<<endl;

		/*cout<<"========================================================================================================================="<<endl;
		cout<<" i       j      dij(min dist between i and j in pbc)                                                                     "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<nsites;i++)
		{
			for (int j=0;j<nsites;j++) cout<<boost::format("%3d    %3d    %+ .5f") %i %j %dmat(i,j)<<endl;
		}*/
	}
}

////////////////////////////////////////////////////////////////////////
void make_J_mats(double J1,double J2,double J3,double J4, 
		 RMatrix &Jmat01, RMatrix &Jmat10,
		 RMatrix &Jmat02, RMatrix &Jmat20,
		 RMatrix &Jmat03, RMatrix &Jmat30,
		 RMatrix &Jmat12, RMatrix &Jmat21,
		 RMatrix &Jmat13, RMatrix &Jmat31,
		 RMatrix &Jmat23, RMatrix &Jmat32)
{
	Jmat01.resize(3,3);
	Jmat10.resize(3,3);
	Jmat02.resize(3,3);
	Jmat20.resize(3,3);
	Jmat03.resize(3,3);
	Jmat30.resize(3,3);
	Jmat12.resize(3,3);
	Jmat21.resize(3,3);
	Jmat13.resize(3,3);
	Jmat31.resize(3,3);
	Jmat23.resize(3,3);
	Jmat32.resize(3,3);

	Jmat01(0,0)=+J2; Jmat01(0,1)=+J4; Jmat01(0,2)=+J4;
	Jmat01(1,0)=-J4; Jmat01(1,1)=+J1; Jmat01(1,2)=+J3;
	Jmat01(2,0)=-J4; Jmat01(2,1)=+J3; Jmat01(2,2)=+J1;
	
	Jmat10(0,0)=+J2;  Jmat10(0,1)=-J4; Jmat10(0,2)=-J4;
	Jmat10(1,0)=+J4;  Jmat10(1,1)=+J1; Jmat10(1,2)=+J3;
	Jmat10(2,0)=+J4;  Jmat10(2,1)=+J3; Jmat10(2,2)=+J1;

	Jmat02(0,0)=+J1; Jmat02(0,1)=-J4; Jmat02(0,2)=+J3;
	Jmat02(1,0)=+J4; Jmat02(1,1)=+J2; Jmat02(1,2)=+J4;
	Jmat02(2,0)=+J3; Jmat02(2,1)=-J4; Jmat02(2,2)=+J1;

	Jmat20(0,0)=+J1;  Jmat20(0,1)=+J4; Jmat20(0,2)=+J3;
	Jmat20(1,0)=-J4;  Jmat20(1,1)=+J2; Jmat20(1,2)=-J4;
	Jmat20(2,0)=+J3;  Jmat20(2,1)=+J4; Jmat20(2,2)=+J1;

	Jmat03(0,0)=+J1; Jmat03(0,1)=+J3; Jmat03(0,2)=-J4;
	Jmat03(1,0)=+J3; Jmat03(1,1)=+J1; Jmat03(1,2)=-J4;
	Jmat03(2,0)=+J4; Jmat03(2,1)=+J4; Jmat03(2,2)=+J2;
	
        Jmat30(0,0)=+J1;  Jmat30(0,1)=+J3; Jmat30(0,2)=+J4;
	Jmat30(1,0)=+J3;  Jmat30(1,1)=+J1; Jmat30(1,2)=+J4;
	Jmat30(2,0)=-J4;  Jmat30(2,1)=-J4; Jmat30(2,2)=+J2;

	Jmat12(0,0)=+J1; Jmat12(0,1)=-J3; Jmat12(0,2)=+J4;
	Jmat12(1,0)=-J3; Jmat12(1,1)=+J1; Jmat12(1,2)=-J4;
	Jmat12(2,0)=-J4; Jmat12(2,1)=+J4; Jmat12(2,2)=+J2;
	
        Jmat21(0,0)=+J1;  Jmat21(0,1)=-J3; Jmat21(0,2)=-J4;
	Jmat21(1,0)=-J3;  Jmat21(1,1)=+J1; Jmat21(1,2)=+J4;
	Jmat21(2,0)=+J4;  Jmat21(2,1)=-J4; Jmat21(2,2)=+J2;

	Jmat13(0,0)=+J1; Jmat13(0,1)=+J4; Jmat13(0,2)=-J3;
	Jmat13(1,0)=-J4; Jmat13(1,1)=+J2; Jmat13(1,2)=+J4;
	Jmat13(2,0)=-J3; Jmat13(2,1)=-J4; Jmat13(2,2)=+J1;

	Jmat31(0,0)=+J1;  Jmat31(0,1)=-J4; Jmat31(0,2)=-J3;
	Jmat31(1,0)=+J4;  Jmat31(1,1)=+J2; Jmat31(1,2)=-J4;
	Jmat31(2,0)=-J3;  Jmat31(2,1)=+J4; Jmat31(2,2)=+J1;

	Jmat23(0,0)=+J2; Jmat23(0,1)=-J4; Jmat23(0,2)=+J4;
	Jmat23(1,0)=+J4; Jmat23(1,1)=+J1; Jmat23(1,2)=-J3;
	Jmat23(2,0)=-J4; Jmat23(2,1)=-J3; Jmat23(2,2)=+J1;
	
	Jmat32(0,0)=+J2;  Jmat32(0,1)=+J4; Jmat32(0,2)=-J4;
	Jmat32(1,0)=-J4;  Jmat32(1,1)=+J1; Jmat32(1,2)=-J3;
	Jmat32(2,0)=+J4;  Jmat32(2,1)=-J3; Jmat32(2,2)=+J1;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void make_g_mats(double gxy,double gz,RMatrix &gmat0, RMatrix &gmat1, RMatrix &gmat2, RMatrix &gmat3)
{

	gmat0.resize(3,3);
	gmat1.resize(3,3);
	gmat2.resize(3,3);
	gmat3.resize(3,3);

	double gp=(2.0*gxy+ gz)/3.0;
	double gm=(gxy - gz)/3.0;

	gmat0(0,0)=gp;gmat0(0,1)=-gm;gmat0(0,2)=-gm;
	gmat0(1,0)=-gm;gmat0(1,1)=gp;gmat0(1,2)=-gm;
	gmat0(2,0)=-gm;gmat0(2,1)=-gm;gmat0(2,2)=gp;

	gmat1(0,0)=gp;gmat1(0,1)=gm;gmat1(0,2)=gm;
	gmat1(1,0)=gm;gmat1(1,1)=gp;gmat1(1,2)=-gm;
	gmat1(2,0)=gm;gmat1(2,1)=-gm;gmat1(2,2)=gp;

	gmat2(0,0)=gp;gmat2(0,1)=gm;gmat2(0,2)=-gm;
	gmat2(1,0)=gm;gmat2(1,1)=gp;gmat2(1,2)=gm;
	gmat2(2,0)=-gm;gmat2(2,1)=gm;gmat2(2,2)=gp;

	gmat3(0,0)=gp;gmat3(0,1)=-gm;gmat3(0,2)=gm;
	gmat3(1,0)=-gm;gmat3(1,1)=gp;gmat3(1,2)=gm;
	gmat3(2,0)=gm;gmat3(2,1)=gm;gmat3(2,2)=gp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void local_j_energy(double &spin, int &site, int &t, double &sx, double &sy, double &sz,
		      std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors,
		      std::vector< std::vector<int> > &nneighbors,
		      RMatrix &Jmat01, RMatrix &Jmat10,
		      RMatrix &Jmat02, RMatrix &Jmat20,
		      RMatrix &Jmat03, RMatrix &Jmat30,
		      RMatrix &Jmat12, RMatrix &Jmat21,
		      RMatrix &Jmat13, RMatrix &Jmat31,
		      RMatrix &Jmat23, RMatrix &Jmat32,
		      double &Jnnn, RMatrix &bond_disorder_matrix,
		      std::vector< std::vector<int> > &ijkt,
		      double &eff_field_x,
		      double &eff_field_y,
		      double &eff_field_z)
{
	//cout<<"Calculating local j energy"<<endl;
	//local_energy=0.0;
	eff_field_x=0.0;
	eff_field_y=0.0;
	eff_field_z=0.0;
	RMatrix Jmat;
	// Disorder free nearest neighbor terms
	for (int j=0;j<neighbors[site].size();j++)
	{
		int k=neighbors[site][j];
		int t2=ijkt[k][3];
		if (t==0 and t2==1) Jmat=Jmat01;
		if (t==1 and t2==0) Jmat=Jmat10;
		if (t==0 and t2==2) Jmat=Jmat02;
		if (t==2 and t2==0) Jmat=Jmat20;
		if (t==0 and t2==3) Jmat=Jmat03;
		if (t==3 and t2==0) Jmat=Jmat30;
		if (t==1 and t2==2) Jmat=Jmat12;
		if (t==2 and t2==1) Jmat=Jmat21;
		if (t==1 and t2==3) Jmat=Jmat13;
		if (t==3 and t2==1) Jmat=Jmat31;
		if (t==2 and t2==3) Jmat=Jmat23;
		if (t==3 and t2==2) Jmat=Jmat32;
		eff_field_x+=(Jmat(0,0)*configx[k])+(Jmat(0,1)*configy[k])+(Jmat(0,2)*configz[k]);
		eff_field_y+=(Jmat(1,0)*configx[k])+(Jmat(1,1)*configy[k])+(Jmat(1,2)*configz[k]);
		eff_field_z+=(Jmat(2,0)*configx[k])+(Jmat(2,1)*configy[k])+(Jmat(2,2)*configz[k]);
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
	//local_energy=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; // Pairs only counted once 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double local_h_energy(double &spin, int &site, int &t, double &sx, double &sy, double &sz,
		      double &hx, double &hy, double &hz, 
		      std::vector<RMatrix> &gmats)
{
        double mu_b=5.7883818012*0.01; // meV/Tesla
	RMatrix gmat;
	RMatrix spinvec(3,1);
	RMatrix hvec(3,1);
	hvec(0,0)=hx;
	hvec(1,0)=hy;
	hvec(2,0)=hz;
	double e=0;
	spinvec(0,0)=sx;
	spinvec(1,0)=sy;
	spinvec(2,0)=sz;
	e=e+(yt_A_x(hvec,gmats[t],spinvec));
	return -e*mu_b*spin;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_h_energy(double &spin, 
		      std::vector<double>  &configx, 
		      std::vector<double>  &configy, 
		      std::vector<double>  &configz,
		      double &hx, double &hy, double &hz, 
		      std::vector<RMatrix> &gmats,
		      std::vector< std::vector<int> > & ijkt)
{
        double mu_b=5.7883818012*0.01; // meV/Tesla
	RMatrix gmat;
	RMatrix spinvec(3,1);
	RMatrix hvec(3,1);
	hvec(0,0)=hx;
	hvec(1,0)=hy;
	hvec(2,0)=hz;
	double e=0;
	int nsites=int(configx.size());
	for (int site=0;site<nsites;site++)
	{
		int t=ijkt[site][3];
		spinvec(0,0)=configx[site];
		spinvec(1,0)=configy[site];
		spinvec(2,0)=configz[site];
		e=e+(yt_A_x(hvec,gmats[t],spinvec));
	}
	return -e*mu_b*spin;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double total_j_energy(double &spin, 
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
		        std::vector< std::vector<int> > & ijkt)
{
	/*RMatrix Jmat;
	RMatrix vecl(3,1);
	RMatrix vecr(3,1);*/

	double e=0;
	// Nearest neighbor disorder-free contributions
	# pragma omp parallel for default(shared) reduction(+:e)
	for (int i=0;i<neighbors.size();i++)
	{
		RMatrix Jmat;
		RMatrix vecl(3,1);
		RMatrix vecr(3,1);
		int t1=ijkt[i][3];
		vecl(0,0)=configx[i];vecl(1,0)=configy[i];vecl(2,0)=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			int t2=ijkt[k][3];
			vecr(0,0)=configx[k];vecr(1,0)=configy[k];vecr(2,0)=configz[k];
			if (t1==0 and t2==1) Jmat=Jmat01;
			if (t1==1 and t2==0) Jmat=Jmat10; 
			if (t1==0 and t2==2) Jmat=Jmat02; 
			if (t1==2 and t2==0) Jmat=Jmat20;
			if (t1==0 and t2==3) Jmat=Jmat03;
			if (t1==3 and t2==0) Jmat=Jmat30;
			if (t1==1 and t2==2) Jmat=Jmat12;
			if (t1==2 and t2==1) Jmat=Jmat21;
			if (t1==1 and t2==3) Jmat=Jmat13;
			if (t1==3 and t2==1) Jmat=Jmat31;
			if (t1==2 and t2==3) Jmat=Jmat23;
			if (t1==3 and t2==2) Jmat=Jmat32;
			e+=(yt_A_x(vecl,Jmat,vecr));
		}
	}

	// Nearest neighbor bond disorder Heisenberg contributions
	# pragma omp parallel for default(shared) reduction(+:e)
	for (int i=0;i<neighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e+=(bond_disorder_matrix(i,k)*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		}
	}
	// Next nn Heisenberg 
	# pragma omp parallel for default(shared) reduction(+:e)
	for (int i=0;i<nneighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<nneighbors[i].size();j++)
		{
			int k=nneighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e+=(Jnnn*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		
		}
	}
	return (e/2.0)*spin*spin;  // Pairs counted twice so divide by 2
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*double total_j_energy(double &spin, std::vector<double> &configx, std::vector<double> &configy, std::vector<double> &configz, 
		      std::vector< std::vector<int> > &neighbors, 
		      std::vector< std::vector<int> > &nneighbors, 
		      RMatrix &Jmat01, RMatrix &Jmat10,
		      RMatrix &Jmat02, RMatrix &Jmat20,
		      RMatrix &Jmat03, RMatrix &Jmat30,
		      RMatrix &Jmat12, RMatrix &Jmat21,
		      RMatrix &Jmat13, RMatrix &Jmat31,
		      RMatrix &Jmat23, RMatrix &Jmat32,
		      double &Jnnn, RMatrix &bond_disorder_matrix,
		      std::vector< std::vector<int> > & ijkt)
{
	RMatrix Jmat;
	RMatrix vecl(3,1);
	RMatrix vecr(3,1);
	double e=0;

	// Nearest neighbor disorder-free contributions
	for (int i=0;i<neighbors.size();i++)
	{
		int t1=ijkt[i][3];
		vecl(0,0)=configx[i];vecl(1,0)=configy[i];vecl(2,0)=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			int t2=ijkt[k][3];
			vecr(0,0)=configx[k];vecr(1,0)=configy[k];vecr(2,0)=configz[k];
			if (t1==0 and t2==1) e=e+(yt_A_x(vecl,Jmat01,vecr));
			if (t1==1 and t2==0) e=e+(yt_A_x(vecl,Jmat10,vecr));
			if (t1==0 and t2==2) e=e+(yt_A_x(vecl,Jmat02,vecr));
			if (t1==2 and t2==0) e=e+(yt_A_x(vecl,Jmat20,vecr));
			if (t1==0 and t2==3) e=e+(yt_A_x(vecl,Jmat03,vecr));
			if (t1==3 and t2==0) e=e+(yt_A_x(vecl,Jmat30,vecr));
			if (t1==1 and t2==2) e=e+(yt_A_x(vecl,Jmat12,vecr));
			if (t1==2 and t2==1) e=e+(yt_A_x(vecl,Jmat21,vecr));
			if (t1==1 and t2==3) e=e+(yt_A_x(vecl,Jmat13,vecr));
			if (t1==3 and t2==1) e=e+(yt_A_x(vecl,Jmat31,vecr));
			if (t1==2 and t2==3) e=e+(yt_A_x(vecl,Jmat23,vecr));
			if (t1==3 and t2==2) e=e+(yt_A_x(vecl,Jmat32,vecr));
		}
	}

	// Nearest neighbor bond disorder Heisenberg contributions
	for (int i=0;i<neighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<neighbors[i].size();j++)
		{
			int k=neighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e=e+(bond_disorder_matrix(i,k)*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		}
	}
	// Next nn Heisenberg 
	for (int i=0;i<nneighbors.size();i++)
	{
		double sx1=configx[i];double sy1=configy[i];double sz1=configz[i];
		for (int j=0;j<nneighbors[i].size();j++)
		{
			int k=nneighbors[i][j];
			double sx2=configx[k];double sy2=configy[k];double sz2=configz[k];
			e=e+(Jnnn*((sx1*sx2) + (sy1*sy2) + (sz1*sz2)));
		
		}
	}
	return (e/2.0)*spin*spin;  // Pairs counted twice so divide by 2
}*/

/////////////////////////////////////////////////////////////////////////////////
std::vector<double> local_magnetization(double 				&spin,
					int 				&t,
					double 				&sx, 
					double 				&sy, 
					double 				&sz,
					std::vector< RMatrix> 		&gmats)
{
	double  mx=0;
	double  my=0;
	double  mz=0;
	std::vector<double> m;
	mx+=(gmats[t](0,0)*sx)+(gmats[t](0,1)*sy)+(gmats[t](0,2)*sz);
	my+=(gmats[t](1,0)*sx)+(gmats[t](1,1)*sy)+(gmats[t](1,2)*sz);
	mz+=(gmats[t](2,0)*sx)+(gmats[t](2,1)*sy)+(gmats[t](2,2)*sz);
	
	m.push_back(mx*spin);
	m.push_back(my*spin);
	m.push_back(mz*spin);
	return m;
}
/////////////////////////////////////////////////////////////////////////////////
std::vector<double> total_magnetization(double 				&spin,
					std::vector<double> 		&configx, 
					std::vector<double> 		&configy, 
					std::vector<double> 		&configz,
					std::vector<RMatrix>		&gmats,
					std::vector< std::vector<int> > &ijkt)
{
	double mx=0;
	double my=0;
	double mz=0;
	std::vector<double> m;
	int nsites=int(configx.size());
	for (int n=0;n<nsites;n++)
	{
		int t=ijkt[n][3];
		double sx=configx[n];
		double sy=configy[n];
		double sz=configz[n];
		mx+=(gmats[t](0,0)*sx)+(gmats[t](0,1)*sy)+(gmats[t](0,2)*sz);
		my+=(gmats[t](1,0)*sx)+(gmats[t](1,1)*sy)+(gmats[t](1,2)*sz);
		mz+=(gmats[t](2,0)*sx)+(gmats[t](2,1)*sy)+(gmats[t](2,2)*sz);
	}
	m.push_back(mx*spin);
	m.push_back(my*spin);
	m.push_back(mz*spin);
	return m;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void iterative_pyrochlore(double spin, int L,int nsamples, int nburn, string start_config, 
		   	  double hx, double hy, double hz, 
		   	  double J1, double J2, double J3, double J4, double Jnnn,
		   	  double disorder_strength,
		   	  double gxy, double gz, 
		   	  double & eavg, 
		   	  double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   	  double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,measure_corrs);
	//Matrix phases;
	//make_phases(fullcoords,qvals,phases);
	// Make disorder matrix to be added to the Heisenberg couplings 
	make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	std::vector<QMC_Info> infos;
	int b=0;
	QMC_Info qmc;
	infos.push_back(qmc);
	cout<<"Beta = "<<infos[b].beta<<endl;
	if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
	if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
	if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);

	double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
							 neighbors,nneighbors,
							 Jmat01, Jmat10, 
							 Jmat02, Jmat20, 
							 Jmat03, Jmat30, 
							 Jmat12, Jmat21, 
							 Jmat13, Jmat31, 
							 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
	
	infos[b].energy=energyj;
	cout<<"Total J energy                  ="<<energyj<<endl;
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	for (int n=0; n<(nsamples+nburn);n++)
	{
		int site=uniform_rand_int(0,nsites);
		double sx=infos[b].configx[site];
		double sy=infos[b].configy[site];
		double sz=infos[b].configz[site];
		int t=ijkt[site][3];
		// Calculate local energy of old and new configs
		double eff_field_x, eff_field_y, eff_field_z;
		local_j_energy(spin, site,t,sx,sy,sz,
		infos[b].configx,infos[b].configy,infos[b].configz, 
		neighbors, nneighbors,
		Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
		Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
		eff_field_x, eff_field_y,eff_field_z);
		double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
		norm=sqrt(norm);
		double factor=0.0;
		//if (n<nburn) factor=0.01;
		double r=0.5*uniform_rnd();
		double sxnew=((-eff_field_x/norm)*spin)*r + sx + factor*(2.0*uniform_rnd()-1);
		double synew=((-eff_field_y/norm)*spin)*r + sy + factor*(2.0*uniform_rnd()-1);
		double sznew=((-eff_field_z/norm)*spin)*r + sz + factor*(2.0*uniform_rnd()-1);
		norm=(sxnew*sxnew) + (synew*synew) + (sznew*sznew);
		norm=sqrt(norm);
		sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;
		double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
		double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
		if (local_energyj2<local_energyj1)
		{
			infos[b].energy+=(local_energyj2-local_energyj1);
			infos[b].configx[site]=sxnew;
			infos[b].configy[site]=synew;
			infos[b].configz[site]=sznew;
		}
		if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}	
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mc_pyrochlore_get_thermal_config_and_time_evolve(double spin, double deltat, 
		   double omegamin, double omegamax, double omegaspacing, 
		   double tottime, int L, int nstarts, int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz) 
{
	STensor smunu, savg;
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	double tempKelvin=temp;
        temp=tempKelvin*0.08621738;

	std::vector<double> omegas;
	int numomegas=int(((omegamax-omegamin)/omegaspacing) + 1 + 1.0e-6);
	for (int i=0;i<numomegas;i++) 
	{
		double omega=omegamin+(double(i)*omegaspacing);
		omegas.push_back(omega);
	}
	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
       
	//////////////////////////////////////
	// Make pyrochlore lattice 
	time_t start,end; 
	time (&start);	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,false);
	// Make disorder matrix to be added to the Heisenberg couplings 
	make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	time (&end);	
	double seconds=difftime(end,start);
    	cout<<"Time to make pyrochlore = "<<seconds<<" seconds"<<endl;
	time (&start);	
	//////////////////////////////////////
	
	
	/////////////////////////////////////
	// Make phases and q values	
	std::vector<std::vector<double> > qvals;
	make_qs(L,qvals);
	int numqs=qvals.size();
	Matrix phases;
	make_phases(fullcoords, qvals,phases);
	time (&end);	
	seconds=difftime(end,start);
    	cout<<"Time to make q, phases = "<<seconds<<" seconds"<<endl;
	/////////////////////////////////////

	for (int nstart=0;nstart<nstarts;nstart++) // Number of times i.e. starts to run MC
	{	
		time(&start);
		/////////////////////////////////////////////////////////////////////////
		// Make random configuration of spins or selected type 
		ntemps=1;
		std::vector<QMC_Info> infos;
		cout<<"Ntemps   = "<<ntemps<<endl;
		double exponent=pow(100.0,1.0/double(ntemps)); // Tmax/Tmin=30
		cout<<"exponent = "<<exponent<<endl;
		for (int b=0;b<ntemps;b++)
		{
			QMC_Info qmc;
			qmc.init(L,nsites,tempKelvin*pow(exponent,double(b)),false);
			infos.push_back(qmc);
			cout<<"Beta = "<<infos[b].beta<<endl;
			if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
			if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
			if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
		}
		for (int b=0;b<ntemps;b++)
		{
			//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
			double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
									 neighbors,nneighbors,
									 Jmat01, Jmat10, 
									 Jmat02, Jmat20, 
									 Jmat03, Jmat30, 
									 Jmat12, Jmat21, 
									 Jmat13, Jmat31, 
									 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
			double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
								   hx,hy,hz, gmats , ijkt);
			std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
									      gmats, ijkt);
			infos[b].energy=energyj+energyh;
			infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
			cout<<"Total J energy                  ="<<energyj<<endl;
			cout<<"Total h energy                  ="<<energyh<<endl;
			cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
			cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
		}

		cout<<endl;	
		cout<<endl;	
		cout<<"========================================================"<<endl;
		cout<<"Energy history "<<endl;	
		/////////////////////////////////////////////////////////////////////////
		// Accept reject Metropolis
		double accept=0.0;
		double reject=0.0;
		double nswaps=0.0;
		double nreplicatries=0.0;
		double elowest=0.0;
		for (int64_t n=0; n<(nsamples+nburn);n++) // Begin MC
		{
			///////////////////////////////////////////////////////////////////////
			// Usual Moves of a serial Metropolis Monte Carlo
			///////////////////////////////////////////////////////////////////////
			int b=0;
			// Very Small moves needed at low temperatures to increase acceptance rates
			double move_size=0.3;  
			// Choose random site
			int site=uniform_rand_int(0,nsites);
			int t=ijkt[site][3];
			double sxnew,synew,sznew;
			// Current sx,sy,sz on chosen site
			double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
			if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
			// Choose a completely random direction within a cone
			if (mcmove=="conical") conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
			if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
			if (mcmove=="largeD")
			{
				double tempnum=uniform_rnd();
				if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
				else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
			}	
			// Normalize new direction
			double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
			sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

			// Calculate local energy of old and new configs
			double eff_field_x, eff_field_y, eff_field_z;
			local_j_energy(spin, site,t,sx,sy,sz,
				       infos[b].configx,infos[b].configy,infos[b].configz, 
				       neighbors, nneighbors,
				       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
				       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
				       eff_field_x, eff_field_y,eff_field_z);
			double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
			double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
			double mxdiff=(sxnew-sx)*spin;	
			double mydiff=(synew-sy)*spin;
			double mzdiff=(sznew-sz)*spin;
					// Jterms                        // hterms - no field for now
			double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
			double beta=infos[b].beta; 
			if (n<nburn/2) beta=infos[b].beta*(2.0*double(n)/double(nburn));
			double prob=exp(-beta*ediff);
			double rand=uniform_rnd();
			if (rand<prob) // Metropolis Accept-reject for a given temperature
			{
				if (b==0 and n>(nburn)) accept=accept+1.0;
				// Reset configs to new configs, because accepted
				infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
				infos[b].energy+=ediff;
				if (infos[b].energy<elowest) elowest=infos[b].energy;
				infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
			}
			else
			{
				if (b==0 and n>(nburn)) reject=reject+1.0;
				// Still move spin if rejected, but this conserves energy!!
				double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
				norm=sqrt(norm);
				eff_field_x=eff_field_x/norm;
				eff_field_y=eff_field_y/norm;
				eff_field_z=eff_field_z/norm;
				double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
				double sintheta=sqrt(1.0-(costheta*costheta));	
				double phi=uniform_rnd()*2.0*3.14159; // random in range 0 to 2 pi
				double x1,y1,z1, x2,y2,z2;	
				get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
							     x1, y1, z1, 
							     x2, y2, z2);
				double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
				double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
				double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
				double mxdiff=(sxnew-sx)*spin;	
				double mydiff=(synew-sy)*spin;
				double mzdiff=(sznew-sz)*spin;
				infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
				infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
				// ediff=0 by construction!!
			}
			if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .10f") %infos[0].energy<<endl;
		} // End nsamples and i.e. MC run
		cout<<"========================================================"<<endl;
		cout<<endl;
		time_evolve(spin, deltat, tottime, L, infos[0].configx,infos[0].configy,infos[0].configz,neighbors, 
			    nneighbors, Jmat01, Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
			    Jmat13, Jmat31,Jmat23, Jmat32, Jnnn, bond_disorder_matrix,fullcoords,ijkt, qvals, phases, omegas, smunu);
	     
                if (nstart>0) savg.update_totals(smunu);
		else	 { savg.init(int(smunu.qvals.size()),numomegas); savg.copy(smunu);}
		time(&end);
		double seconds=difftime(end,start);
		cout<<"Done with run "<<nstart<<" in  "<<seconds<<" seconds"<<endl;
		cout<<"======================================================================================================================================="<<endl;
		cout<<endl;
        }  // End nstarts
	
        savg.average();
	for (int om=0; om<numomegas;om++)
	{   
		cout<<"======================================================================================================================================="<<endl;
		cout<<"                                              AVERAGED DATA for omega = "<<omegas[om]<<"                                               "<<endl;
		cout<<"======================================================================================================================================="<<endl;
		cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)          Sperp(Q)     "<<endl;
		cout<<"======================================================================================================================================="<<endl;
		for (int i=0;i<numqs;i++)
		{
			double qx=savg.qvals[i][0];double qy=savg.qvals[i][1];double qz=savg.qvals[i][2];
			complex<double> sxx=savg.sxx[om][i];complex<double> syy=savg.syy[om][i];complex<double> szz=savg.szz[om][i];
			complex<double> sxy=savg.sxy[om][i];complex<double> syx=savg.syx[om][i];
			complex<double> sxz=savg.sxz[om][i];complex<double> szx=savg.szx[om][i];
			complex<double> syz=savg.syz[om][i];complex<double> szy=savg.szy[om][i];
			complex<double> sperp=savg.sperp[om][i];
			
			cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f  %+ .8f  %+ .8f %+ .8f  %+ .8f") %qx %qy %qz %sxx %sxy %sxz %syx %syy %syz %szx %szy %szz %sperp<<endl;
		 }				
	 }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void mc_pyrochlore_get_thermal_config_and_spin_wave(double spin,
//		  int L, int nstarts, int64_t nsamples, int64_t nburn, string start_config, 
//		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
//		   double J1, double J2, double J3, double J4, double Jnnn,
//		   double disorder_strength,
//		   double gxy, double gz) 
//{
//	/////////////////////////////////////////////////////////////////////////
//	// MC related quantities
//	// Units - Assume J in meV (millielectron volts)
//        //         Assume h in T   (tesla)
//        //         Convert h to meV units in total_h_energy and local_h_energy
//        // Convert temperature in Kelvin to temperature in meV
//	// Set 16xLxLxL pyrochlore lattice
//	int nsites=16*L*L*L;
//	//int nsites=4;
//       	double kB=1.38064852*1e-23;
//        double NA=6.02214179*1e23;
//        double JpermeV=1.60218*1e-22; 
//	double tempKelvin=temp;
//        temp=tempKelvin*0.08621738;
//
//	/////////////////////////////////////////////////////////////////////////
//	// J and g matrices
//	
//	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
//	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
//	RMatrix gmat0,gmat1,gmat2,gmat3; 
//        RMatrix bond_disorder_matrix(nsites,nsites);
//	// Given J1, J2, J3, J4 - make J mats 
//	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
//	// Given g's - make g mats 
//	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
//	std::vector<RMatrix> gmats;
//	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);
//
//	/////////////////////////////////////////////////////////////////////////
//	// MC related quantities
//	// Set nburn
//	nburn=nsamples;
//	cout<<"Nsamples = "<<nsamples<<endl;
//	cout<<"Nburn    = "<<nburn<<endl;
//	double e4avg,mx4avg,my4avg,mz4avg;
//	
//	/////////////////////////////////////////////////////////////////////////
//	// Set 16xLxLxL pyrochlore lattice
//	std::vector< std::vector<int> > neighbors(nsites);	
//	std::vector< std::vector<int> > nneighbors(nsites);	
//	std::vector< std::vector<double> > fullcoords;	
//	std::vector< std::vector<int> > ijkt;	
//       
//	//////////////////////////////////////
//	// Make pyrochlore lattice 
//	time_t start,end; 
//	time (&start);	
//        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,false);
//        //make_single_tet(L,fullcoords,ijkt,neighbors,nneighbors,false);
//	// Make disorder matrix to be added to the Heisenberg couplings 
//	make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
//	time (&end);	
//	double seconds=difftime(end,start);
//    	cout<<"Time to make pyrochlore = "<<seconds<<" seconds"<<endl;
//	time (&start);	
//	//////////////////////////////////////
//	
//	
//	/////////////////////////////////////
//	// Make phases and q values	
//	std::vector<std::vector<double> > qvals;
//	make_qs(L,qvals);
//	int numqs=qvals.size();
//	Matrix phases;
//	make_phases(fullcoords, qvals,phases);
//	time (&end);	
//	seconds=difftime(end,start);
//    	cout<<"Time to make q, phases = "<<seconds<<" seconds"<<endl;
//	/////////////////////////////////////
//
//	for (int nstart=0;nstart<nstarts;nstart++) // Number of times i.e. starts to run MC
//	{	
//		time(&start);
//		/////////////////////////////////////////////////////////////////////////
//		// Make random configuration of spins or selected type 
//		ntemps=1;
//		std::vector<QMC_Info> infos;
//		cout<<"Ntemps   = "<<ntemps<<endl;
//		double exponent=pow(100.0,1.0/double(ntemps)); // Tmax/Tmin=30
//		cout<<"exponent = "<<exponent<<endl;
//		for (int b=0;b<ntemps;b++)
//		{
//			QMC_Info qmc;
//			qmc.init(L,nsites,tempKelvin*pow(exponent,double(b)),false);
//			infos.push_back(qmc);
//			cout<<"Beta = "<<infos[b].beta<<endl;
//			if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
//			if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
//			if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
//		}
//		for (int b=0;b<ntemps;b++)
//		{
//			//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
//			double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
//									 neighbors,nneighbors,
//									 Jmat01, Jmat10, 
//									 Jmat02, Jmat20, 
//									 Jmat03, Jmat30, 
//									 Jmat12, Jmat21, 
//									 Jmat13, Jmat31, 
//									 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
//			double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
//								   hx,hy,hz, gmats , ijkt);
//			std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
//									      gmats, ijkt);
//			infos[b].energy=energyj+energyh;
//			infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
//			cout<<"Total J energy                  ="<<energyj<<endl;
//			cout<<"Total h energy                  ="<<energyh<<endl;
//			cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
//			cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
//		}
//
//		cout<<endl;	
//		cout<<endl;	
//		cout<<"========================================================"<<endl;
//		cout<<"Energy history "<<endl;	
//		/////////////////////////////////////////////////////////////////////////
//		// Accept reject Metropolis
//		double accept=0.0;
//		double reject=0.0;
//		double nswaps=0.0;
//		double nreplicatries=0.0;
//		double elowest=0.0;
//		for (int64_t n=0; n<(nsamples+nburn);n++) // Begin MC
//		{
//			///////////////////////////////////////////////////////////////////////
//			// Usual Moves of a serial Metropolis Monte Carlo
//			///////////////////////////////////////////////////////////////////////
//			int b=0;
//			// Very Small moves needed at low temperatures to increase acceptance rates
//			double move_size=0.3;  
//			// Choose random site
//			int site=uniform_rand_int(0,nsites);
//			int t=ijkt[site][3];
//			double sxnew,synew,sznew;
//			// Current sx,sy,sz on chosen site
//			double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
//			if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
//			// Choose a completely random direction within a cone
//			if (mcmove=="conical") conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
//			if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
//			if (mcmove=="largeD")
//			{
//				double tempnum=uniform_rnd();
//				if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
//				else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
//			}	
//			// Normalize new direction
//			double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
//			sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;
//
//			// Calculate local energy of old and new configs
//			double eff_field_x, eff_field_y, eff_field_z;
//			local_j_energy(spin, site,t,sx,sy,sz,
//				       infos[b].configx,infos[b].configy,infos[b].configz, 
//				       neighbors, nneighbors,
//				       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
//				       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
//				       eff_field_x, eff_field_y,eff_field_z);
//			double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
//			double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
//			double mxdiff=(sxnew-sx)*spin;	
//			double mydiff=(synew-sy)*spin;
//			double mzdiff=(sznew-sz)*spin;
//					// Jterms                        // hterms - no field for now
//			double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
//			double beta=infos[b].beta; 
//			if (n<nburn/2) beta=infos[b].beta*(2.0*double(n)/double(nburn));
//			double prob=exp(-beta*ediff);
//			double rand=uniform_rnd();
//			if (rand<prob) // Metropolis Accept-reject for a given temperature
//			{
//				if (b==0 and n>(nburn)) accept=accept+1.0;
//				// Reset configs to new configs, because accepted
//				infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
//				infos[b].energy+=ediff;
//				if (infos[b].energy<elowest) elowest=infos[b].energy;
//				infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
//			}
//			else
//			{
//				if (b==0 and n>(nburn)) reject=reject+1.0;
//				// Still move spin if rejected, but this conserves energy!!
//				double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
//				norm=sqrt(norm);
//				eff_field_x=eff_field_x/norm;
//				eff_field_y=eff_field_y/norm;
//				eff_field_z=eff_field_z/norm;
//				double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
//				double sintheta=sqrt(1.0-(costheta*costheta));	
//				double phi=uniform_rnd()*2.0*3.14159; // random in range 0 to 2 pi
//				double x1,y1,z1, x2,y2,z2;	
//				get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
//							     x1, y1, z1, 
//							     x2, y2, z2);
//				double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
//				double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
//				double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
//				double mxdiff=(sxnew-sx)*spin;	
//				double mydiff=(synew-sy)*spin;
//				double mzdiff=(sznew-sz)*spin;
//				infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
//				infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
//				// ediff=0 by construction!!
//			}
//			if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
//		} // End nsamples and i.e. MC run
//		cout<<"========================================================"<<endl;
//		cout<<endl;
//		spin_wave(spin, L, infos[0].configx,infos[0].configy,infos[0].configz,neighbors, 
//			  nneighbors, Jmat01, Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
//			  Jmat13, Jmat31,Jmat23, Jmat32, Jnnn, bond_disorder_matrix,fullcoords,ijkt);
//                //if (nstart>0) savg.update_totals(smunu);
//		//else	 { savg.init(int(smunu.qvals.size()),numomegas); savg.copy(smunu);}
//		time(&end);
//		double seconds=difftime(end,start);
//		cout<<"Done with run "<<nstart<<" in  "<<seconds<<" seconds"<<endl;
//	}  // End nstarts
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mc_pyrochlore_slowcool(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	double tempKelvin=temp;
        temp=tempKelvin*0.08621738;

	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,measure_corrs);
	// Make disorder matrix to be added to the Heisenberg couplings 
	make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	ntemps=1;
	std::vector<QMC_Info> infos;
	cout<<"Ntemps   = "<<ntemps<<endl;
	double exponent=pow(100.0,1.0/double(ntemps)); // Tmax/Tmin=30
	cout<<"exponent = "<<exponent<<endl;
	for (int b=0;b<ntemps;b++)
	{
		QMC_Info qmc;
		qmc.init(L,nsites,tempKelvin*pow(exponent,double(b)),measure_corrs);
		infos.push_back(qmc);
		cout<<"Beta = "<<infos[b].beta<<endl;
		if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
	}
	Matrix phases;
	make_phases(fullcoords,infos[0].qvals,phases);

	
	for (int b=0;b<ntemps;b++)
	{
		//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
		double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
								 neighbors,nneighbors,
							   	 Jmat01, Jmat10, 
							   	 Jmat02, Jmat20, 
							   	 Jmat03, Jmat30, 
							   	 Jmat12, Jmat21, 
							   	 Jmat13, Jmat31, 
							   	 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
		double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
							   hx,hy,hz, gmats , ijkt);
		std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
								      gmats, ijkt);
		fourier_transforms(spin,phases,infos[b].configx,infos[b].configy,infos[b].configz,infos[b].sxq,infos[b].syq,infos[b].szq);
		infos[b].energy=energyj+energyh;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J energy                  ="<<energyj<<endl;
		cout<<"Total h energy                  ="<<energyh<<endl;
		cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
		cout<<"Total magnetization (from FT is)="<<endl; cout<<infos[b].sxq[0]<<"  "<<infos[b].syq[0]<<"  "<<infos[b].szq[0]<<endl;
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	// Accept reject Metropolis
	double accept=0.0;
	double reject=0.0;
	double nswaps=0.0;
	double nreplicatries=0.0;
	double elowest=0.0;
	for (int64_t n=0; n<(nsamples+nburn);n++)
	{
		///////////////////////////////////////////////////////////////////////
		// Usual Moves of a serial Metropolis Monte Carlo
		///////////////////////////////////////////////////////////////////////
		int b=0;
		// Very Small moves needed at low temperatures to increase acceptance rates
		double move_size=0.3;  
		// Choose random site
		int site=uniform_rand_int(0,nsites);
		int t=ijkt[site][3];
		double sxnew,synew,sznew;
		// Current sx,sy,sz on chosen site
		double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
		if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
		// Choose a completely random direction within a cone
		if (mcmove=="conical") conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
		if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
		if (mcmove=="largeD")
		{
			double tempnum=uniform_rnd();
			if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
			else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
		}	
		// Normalize new direction
		double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
		sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

		// Calculate local energy of old and new configs
		double eff_field_x, eff_field_y, eff_field_z;
		local_j_energy(spin, site,t,sx,sy,sz,
			       infos[b].configx,infos[b].configy,infos[b].configz, 
			       neighbors, nneighbors,
			       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
			       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
			       eff_field_x, eff_field_y,eff_field_z);
		double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
		double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
		double mxdiff=(sxnew-sx)*spin;	
		double mydiff=(synew-sy)*spin;
		double mzdiff=(sznew-sz)*spin;
				// Jterms                        // hterms - no field for now
		double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
		double beta=infos[b].beta; 
		if (n<nburn/2) beta=infos[b].beta*(2.0*double(n)/double(nburn));
		double prob=exp(-beta*ediff);
		double rand=uniform_rnd();
		if (rand<prob) // Metropolis Accept-reject for a given temperature
		{
			if (b==0 and n>(nburn)) accept=accept+1.0;
			// Reset configs to new configs, because accepted
			infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
			infos[b].energy+=ediff;
			if (infos[b].energy<elowest) elowest=infos[b].energy;
			infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
			if (n>nburn and measure_corrs==true) // Updating needed only when measuring
			{
				#pragma omp parallel for
				for (int i=0;i<infos[b].numqs;i++)
				{
					infos[b].sxq[i]+=(mxdiff*phases(i,site));
					infos[b].syq[i]+=(mydiff*phases(i,site));
					infos[b].szq[i]+=(mzdiff*phases(i,site));
				}
			}
		}
		else
		{
			if (b==0 and n>(nburn)) reject=reject+1.0;
			// Still move spin if rejected, but this conserves energy!!
			double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
			norm=sqrt(norm);
			eff_field_x=eff_field_x/norm;
			eff_field_y=eff_field_y/norm;
			eff_field_z=eff_field_z/norm;
			double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
			double sintheta=sqrt(1.0-(costheta*costheta));	
			double phi=uniform_rnd()*2.0*3.14159; // random in range 0 to 2 pi
			double x1,y1,z1, x2,y2,z2;	
			get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
						     x1, y1, z1, 
						     x2, y2, z2);
			double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
			double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
			double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
			double mxdiff=(sxnew-sx)*spin;	
			double mydiff=(synew-sy)*spin;
			double mzdiff=(sznew-sz)*spin;
			infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
			infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
			if (n>nburn and measure_corrs==true) // Updating needed only when measuring
			{
				#pragma omp parallel for
				for (int i=0;i<infos[b].numqs;i++)
				{
					infos[b].sxq[i]+=(mxdiff*phases(i,site));
					infos[b].syq[i]+=(mydiff*phases(i,site));
					infos[b].szq[i]+=(mzdiff*phases(i,site));
				}
			}
			// ediff=0 by construction!!
		}
		// Structure factors not needed before measurement stage
		if (n==nburn and measure_corrs==true) fourier_transforms(spin,phases,infos[b].configx,infos[b].configy,infos[b].configz,infos[b].sxq,infos[b].syq,infos[b].szq);
		if (n>nburn and n%nsites==0) infos[0].update_totals(); // Update averages - after 1 sweep and after equilibration done 
		if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}
	cout<<"========================================================"<<endl;
	cout<<endl;
	accept=accept/(accept+reject);
	// Only the lowest temperature is relevant for averages we are interested in 
	infos[0].average();
	cout<<"elowest    	= "<<boost::format("%+ .15f") %elowest<<endl;
	cout<<"Nmeas    	= "<<boost::format("%+ .15f") %infos[0].nmeas<<endl;
	cout<<"accept    	= "<<boost::format("%+ .15f") %accept<<endl;
	cout<<"nswaps (0)    	= "<<boost::format("%+ .15f") %nswaps<<endl;
	cout<<"mxavg     	= "<<boost::format("%+ .15f") %infos[0].mxavg<<endl;
	cout<<"myavg     	= "<<boost::format("%+ .15f") %infos[0].myavg<<endl;
	cout<<"mzavg     	= "<<boost::format("%+ .15f") %infos[0].mzavg<<endl;
	cout<<"mx2avg    	= "<<boost::format("%+ .15f") %infos[0].mx2avg<<endl;
	cout<<"my2avg    	= "<<boost::format("%+ .15f") %infos[0].my2avg<<endl;
	cout<<"mz2avg    	= "<<boost::format("%+ .15f") %infos[0].mz2avg<<endl;
	cout<<"mx4avg    	= "<<boost::format("%+ .15f") %infos[0].mx4avg<<endl;
	cout<<"my4avg    	= "<<boost::format("%+ .15f") %infos[0].my4avg<<endl;
	cout<<"mz4avg    	= "<<boost::format("%+ .15f") %infos[0].mz4avg<<endl;
	cout<<"eavg      	= "<<boost::format("%+ .15f") %infos[0].eavg<<endl;
	cout<<"e2avg     	= "<<boost::format("%+ .15f") %infos[0].e2avg<<endl;
	cout<<"e4avg     	= "<<boost::format("%+ .15f") %infos[0].e4avg<<endl;
	cout<<"Cv        	= "<<boost::format("%+ .15f") %infos[0].spheat<<endl;
	cout<<"Cv/Ni(J/molK)	= "<<boost::format("%+ .15f") %infos[0].spheatpersite<<endl;
	cout<<"Cvps(no dim)	= "<<boost::format("%+ .15f") %infos[0].spheatpersitenodim<<endl;

	if (measure_corrs==true)
	{	
		cout<<"========================================================================================================================="<<endl;
		cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)         "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<infos[0].numqs;i++)
		{
		cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f  %+ .5f  %+ .5f %+ .5f") %infos[0].qvals[i][0] %infos[0].qvals[i][1] %infos[0].qvals[i][2] %infos[0].sxsxtot[i] %infos[0].sxsytot[i] %infos[0].sxsztot[i] %infos[0].sysxtot[i] %infos[0].sysytot[i] %infos[0].sysztot[i] %infos[0].szsxtot[i] %infos[0].szsytot[i] %infos[0].szsztot[i]<<endl;
			
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void mc_pyrochlore_ptall(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double tminKelvin, double tmaxKelvin, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,measure_corrs);
	// Make disorder matrix to be added to the Heisenberg couplings 
	// Till here the seed will give the same sequence of random numbers
        make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	std::vector<QMC_Info> infos;
	cout<<"Ntemps   = "<<ntemps<<endl;
	double tmax_over_tmin=tmaxKelvin/tminKelvin;

	double exponent=1.0; 
	if (ntemps>1) exponent=pow(tmax_over_tmin,1.0/double(ntemps-1));
	
	cout<<"exponent = "<<exponent<<endl;
	for (int b=0;b<ntemps;b++)
	{
		QMC_Info qmc;
		qmc.init(L,nsites,tminKelvin*pow(exponent,double(b)),measure_corrs);
		infos.push_back(qmc);
		cout<<"Beta = "<<infos[b].beta<<endl;
		if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
	}
	//Matrix phases;
	//make_phases(fullcoords,infos[0].qvals,phases);

	for (int b=0;b<ntemps;b++)
	{
		//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
		double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
								 neighbors,nneighbors,
							   	 Jmat01, Jmat10, 
							   	 Jmat02, Jmat20, 
							   	 Jmat03, Jmat30, 
							   	 Jmat12, Jmat21, 
							   	 Jmat13, Jmat31, 
							   	 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
		double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
							   hx,hy,hz, gmats , ijkt);
		std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
								      gmats, ijkt);
		infos[b].energy=energyj+energyh;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J energy                  ="<<energyj<<endl;
		cout<<"Total h energy                  ="<<energyh<<endl;
		cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	double nreplicatries=0.0;
	for (int64_t n=0; n<(nsamples+(nburn));n++)
	{
	   if (n%nsites!=0 or ntemps==1) // Do usual Metropolis MC
	   {
			///////////////////////////////////////////////////////////////////////
			// Usual Moves of a serial Metropolis Monte Carlo
			///////////////////////////////////////////////////////////////////////
			std::vector<double> rnd1,rnd2,rnd3,rnd4;
			std::vector<int> rndints;
		        // random numbers generated in advance
			for (int b=0;b<ntemps;b++)
			{
				bool cond=false;
				double r1,r2,d1,d2;
				while (cond==false) // Problem if fixed random numbers given !!!!!!!!!
				{
					r1=uniform_rnd();
					r2=uniform_rnd();
					d1=(2.0*r1 - 1);	
					d2=(2.0*r2 - 1);	
					if (d1*d1 + d2*d2 <=1.0) cond=true; 
				}
				rnd1.push_back(r1);
				rnd2.push_back(r2);
				// First 2 rnds are drawn in a circle for the conical move to work
				rnd3.push_back(uniform_rnd());
				rnd4.push_back(uniform_rnd());
				rndints.push_back(uniform_rand_int(0,nsites));
			}
			# pragma omp parallel for
			for (int b=0;b<ntemps;b++)
			{
					// Very Small moves needed at low temperatures to increase acceptance rates
					//double move_size=min(0.3,0.1*(infos[0].beta)/(infos[b].beta));  
					double move_size=0.2;  
					// Choose random site
					int site=rndints[b];
					int t=ijkt[site][3];
					double sxnew,synew,sznew;
					// Current sx,sy,sz on chosen site
					double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
					//if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
					// Choose a completely random direction within a cone
					if (mcmove=="conical") conical_move_continuous_spin_rnds_provided(move_size,rnd1[b],rnd2[b],sx,sy,sz,sxnew,synew,sznew);
					//if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					/*if (mcmove=="largeD")
					{
						double tempnum=rnd1[b];
						if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
						else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					}*/	
					// Normalize new direction
					double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
					sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

					// Calculate local energy of old and new configs
					double eff_field_x, eff_field_y, eff_field_z;
					local_j_energy(spin, site,t,sx,sy,sz,
						       infos[b].configx,infos[b].configy,infos[b].configz, 
						       neighbors, nneighbors,
						       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
						       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
						       eff_field_x, eff_field_y,eff_field_z);
					double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
					double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
					double mxdiff=(sxnew-sx)*spin;	
					double mydiff=(synew-sy)*spin;
					double mzdiff=(sznew-sz)*spin;
							// Jterms                        // hterms - no field for now
					double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
					double beta; 
					beta=infos[b].beta;
					double prob=exp(-beta*ediff);
					double rand=rnd3[b]; // random number previously generated
					if (rand<prob) // Metropolis Accept-reject for a given temperature
					{
						if (n>(nburn)) infos[b].accept=infos[b].accept+1.0;
						// Reset configs to new configs, because accepted
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].energy+=ediff;
						//if (infos[b].energy<elowest) elowest=infos[b].energy;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
					}
					else
					{
					       	if (n>(nburn)) infos[b].reject=infos[b].reject+1.0;
						// Still move spin if rejected, but this conserves energy!!
						double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
						norm=sqrt(norm);
						eff_field_x=eff_field_x/norm;
						eff_field_y=eff_field_y/norm;
						eff_field_z=eff_field_z/norm;
						double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
						double sintheta=sqrt(1.0-(costheta*costheta));	
						double phi=rnd4[b]*2.0*3.14159; // random in range 0 to 2 pi 
										// random number previously generated
						double x1,y1,z1, x2,y2,z2;	
						get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
		       							     x1, y1, z1, 
		       							     x2, y2, z2);
						double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
						double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
						double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
						double mxdiff=(sxnew-sx)*spin;	
						double mydiff=(synew-sy)*spin;
						double mzdiff=(sznew-sz)*spin;
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
						// ediff=0 by construction!!
					}
			}
	}	
	else if (ntemps>1)  //50 % chance of doing replica exchange
	{
		nreplicatries+=1.0;	
		///////////////////////////////////////////////////////////////////////
		////////// Attempted exchange moves of parallel tempering
		///////////////////////////////////////////////////////////////////////
		// Metropolis move done, try swapping every 2 sweeps
		for (int which1=0;which1<ntemps;which1++)
		{
			int which2=(which1+1)%(ntemps);
			// Slight bias, if i try to swap serially ?
			double rand=uniform_rnd();	
			double beta_i=infos[which1].beta;
			double beta_j=infos[which2].beta;
			double energy_i=infos[which1].energy;
			double energy_j=infos[which2].energy;
			double power=(beta_j-beta_i)*(energy_i - energy_j);
			double ratio=exp(-power);
			if (rand<ratio)
			{
				// SWAP quantities which are being saved
				std::vector<double> tempv;
				double temp_energy, temp_mx, temp_my, temp_mz;
				tempv=infos[which1].configx;	
				infos[which1].configx=infos[which2].configx;	
				infos[which2].configx=tempv;	
				
				tempv=infos[which1].configy;	
				infos[which1].configy=infos[which2].configy;	
				infos[which2].configy=tempv;	

				tempv=infos[which1].configz;	
				infos[which1].configz=infos[which2].configz;	
				infos[which2].configz=tempv;	
				
				temp_energy=infos[which1].energy;	
				infos[which1].energy=infos[which2].energy;	
				infos[which2].energy=temp_energy;	
				
				temp_mx=infos[which1].mx;	
				infos[which1].mx=infos[which2].mx;	
				infos[which2].mx=temp_mx;	
				
				temp_my=infos[which1].my;	
				infos[which1].my=infos[which2].my;	
				infos[which2].my=temp_my;	
				
				temp_mz=infos[which1].mz;	
				infos[which1].mz=infos[which2].mz;	
				infos[which2].mz=temp_mz;
				infos[which1].nswaps+=1.0;	
				infos[which2].nswaps+=1.0;	
				//nswaps+=1.0;	
			}
		}
	}
	if (n%nsites==0 and n>nburn) 
	{
		#pragma omp parallel for 
		for (int num=0;num<infos.size();num++)
		{
			infos[num].update_totals();
		}
	}
	if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}
	cout<<"========================================================"<<endl;
	cout<<endl;
	#pragma omp parallel for 
	for (int num=0;num<infos.size();num++)
	{
		infos[num].average();
	}

	for (int num=0;num<infos.size();num++)
	{	
		cout<<"================================================================================"<<endl;
		cout<<"Nmeas    	= "<<boost::format("%+ .15f") %infos[num].nmeas<<endl;
		cout<<"accept    	= "<<boost::format("%+ .15f") %(infos[num].accept/(infos[num].accept+infos[num].reject))<<endl;
		cout<<"nswaps     	= "<<boost::format("%+ .15f") %(infos[num].nswaps/nreplicatries)<<endl;
		cout<<"mxavg     	= "<<boost::format("%+ .15f") %infos[num].mxavg<<endl;
		cout<<"myavg     	= "<<boost::format("%+ .15f") %infos[num].myavg<<endl;
		cout<<"mzavg     	= "<<boost::format("%+ .15f") %infos[num].mzavg<<endl;
		cout<<"mx2avg    	= "<<boost::format("%+ .15f") %infos[num].mx2avg<<endl;
		cout<<"my2avg    	= "<<boost::format("%+ .15f") %infos[num].my2avg<<endl;
		cout<<"mz2avg    	= "<<boost::format("%+ .15f") %infos[num].mz2avg<<endl;
		cout<<"mx4avg    	= "<<boost::format("%+ .15f") %infos[num].mx4avg<<endl;
		cout<<"my4avg    	= "<<boost::format("%+ .15f") %infos[num].my4avg<<endl;
		cout<<"mz4avg    	= "<<boost::format("%+ .15f") %infos[num].mz4avg<<endl;
		cout<<"eavg      	= "<<boost::format("%+ .15f") %infos[num].eavg<<endl;
		cout<<"e2avg     	= "<<boost::format("%+ .15f") %infos[num].e2avg<<endl;
		cout<<"e4avg     	= "<<boost::format("%+ .15f") %infos[num].e4avg<<endl;
		cout<<"T (K)        	= "<<boost::format("%+ .15f") %infos[num].tempKelvin<<endl;
		cout<<"Cv        	= "<<boost::format("%+ .15f") %infos[num].spheat<<endl;
		cout<<"Cv/Ni(J/molK)	= "<<boost::format("%+ .15f") %infos[num].spheatpersite<<endl;
		cout<<"Cvps(no dim)	= "<<boost::format("%+ .15f") %infos[num].spheatpersitenodim<<endl;
		cout<<"================================================================================"<<endl;
	}
}


void mc_pyrochlore_pt(double spin, int L,int64_t nsamples, int64_t nburn, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	double tempKelvin=temp;
        temp=tempKelvin*0.08621738;

	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=nsamples;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,measure_corrs);
	// Make disorder matrix to be added to the Heisenberg couplings 
	// Till here the seed will give the same sequence of random numbers
        make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	std::vector<QMC_Info> infos;
	cout<<"Ntemps   = "<<ntemps<<endl;
	double tmax=2.0;
	double tmin=tempKelvin;
	double tmax_over_tmin=tmax/tmin;

	double exponent=pow(tmax_over_tmin,1.0/double(ntemps));
	cout<<"exponent = "<<exponent<<endl;
	for (int b=0;b<ntemps;b++)
	{
		QMC_Info qmc;
		qmc.init(L,nsites,tempKelvin*pow(exponent,double(b)),measure_corrs);
		infos.push_back(qmc);
		cout<<"Beta = "<<infos[b].beta<<endl;
		if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
	}
	//Matrix phases;
	//make_phases(fullcoords,infos[0].qvals,phases);

	for (int b=0;b<ntemps;b++)
	{
		//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
		double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
								 neighbors,nneighbors,
							   	 Jmat01, Jmat10, 
							   	 Jmat02, Jmat20, 
							   	 Jmat03, Jmat30, 
							   	 Jmat12, Jmat21, 
							   	 Jmat13, Jmat31, 
							   	 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
		double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
							   hx,hy,hz, gmats , ijkt);
		std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
								      gmats, ijkt);
		infos[b].energy=energyj+energyh;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J energy                  ="<<energyj<<endl;
		cout<<"Total h energy                  ="<<energyh<<endl;
		cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	// Accept reject Metropolis
	double accept=0.0;
	double reject=0.0;
	double nswaps=0.0;
	double nreplicatries=0.0;
	double elowest=0.0;
	
	for (int64_t n=0; n<(nsamples+(nburn));n++)
	{
	   if (n%nsites!=0 or ntemps==1) // Do usual Metropolis MC
	   {
			///////////////////////////////////////////////////////////////////////
			// Usual Moves of a serial Metropolis Monte Carlo
			///////////////////////////////////////////////////////////////////////
			std::vector<double> rnd1,rnd2,rnd3,rnd4;
			std::vector<int> rndints;
		        // random numbers generated in advance
			for (int b=0;b<ntemps;b++)
			{
				bool cond=false;
				double r1,r2,d1,d2;
				while (cond==false) // Problem if fixed random numbers given !!!!!!!!!
				{
					r1=uniform_rnd();
					r2=uniform_rnd();
					d1=(2.0*r1 - 1);	
					d2=(2.0*r2 - 1);	
					if (d1*d1 + d2*d2 <=1.0) cond=true; 
				}
				rnd1.push_back(r1);
				rnd2.push_back(r2);
				// First 2 rnds are drawn in a circle for the conical move to work
				rnd3.push_back(uniform_rnd());
				rnd4.push_back(uniform_rnd());
				rndints.push_back(uniform_rand_int(0,nsites));
			}
			# pragma omp parallel for
			for (int b=0;b<ntemps;b++)
			{
					// Very Small moves needed at low temperatures to increase acceptance rates
					double move_size=min(0.3,0.1*(infos[0].beta)/(infos[b].beta));  
					// Choose random site
					int site=rndints[b];
					int t=ijkt[site][3];
					double sxnew,synew,sznew;
					// Current sx,sy,sz on chosen site
					double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
					//if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
					// Choose a completely random direction within a cone
					if (mcmove=="conical") conical_move_continuous_spin_rnds_provided(move_size,rnd1[b],rnd2[b],sx,sy,sz,sxnew,synew,sznew);
					//if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					/*if (mcmove=="largeD")
					{
						double tempnum=rnd1[b];
						if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
						else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					}*/	
					// Normalize new direction
					double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
					sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

					// Calculate local energy of old and new configs
					double eff_field_x, eff_field_y, eff_field_z;
					local_j_energy(spin, site,t,sx,sy,sz,
						       infos[b].configx,infos[b].configy,infos[b].configz, 
						       neighbors, nneighbors,
						       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
						       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
						       eff_field_x, eff_field_y,eff_field_z);
					double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
					double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
					double mxdiff=(sxnew-sx)*spin;	
					double mydiff=(synew-sy)*spin;
					double mzdiff=(sznew-sz)*spin;
							// Jterms                        // hterms - no field for now
					double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
					double beta; 
					beta=infos[b].beta;
					double prob=exp(-beta*ediff);
					double rand=rnd3[b]; // random number previously generated
					if (rand<prob) // Metropolis Accept-reject for a given temperature
					{
						if (b==0 and n>(nburn)) accept=accept+1.0;
						// Reset configs to new configs, because accepted
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].energy+=ediff;
						//if (infos[b].energy<elowest) elowest=infos[b].energy;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
					}
					else
					{
					       	if (b==0 and n>(nburn)) reject=reject+1.0;
						// Still move spin if rejected, but this conserves energy!!
						double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
						norm=sqrt(norm);
						eff_field_x=eff_field_x/norm;
						eff_field_y=eff_field_y/norm;
						eff_field_z=eff_field_z/norm;
						double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
						double sintheta=sqrt(1.0-(costheta*costheta));	
						double phi=rnd4[b]*2.0*3.14159; // random in range 0 to 2 pi 
										// random number previously generated
						double x1,y1,z1, x2,y2,z2;	
						get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
		       							     x1, y1, z1, 
		       							     x2, y2, z2);
						double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
						double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
						double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
						double mxdiff=(sxnew-sx)*spin;	
						double mydiff=(synew-sy)*spin;
						double mzdiff=(sznew-sz)*spin;
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
						// ediff=0 by construction!!
					}
			}
	}	
	else if (ntemps>1)  //50 % chance of doing replica exchange
	{
		nreplicatries+=1.0;	
		///////////////////////////////////////////////////////////////////////
		////////// Attempted exchange moves of parallel tempering
		///////////////////////////////////////////////////////////////////////
		// Metropolis move done, try swapping every 2 sweeps
		for (int which1=0;which1<ntemps;which1++)
		{
			int which2=(which1+1)%(ntemps);
			// Slight bias, if i try to swap serially ?
			double rand=uniform_rnd();	
			double beta_i=infos[which1].beta;
			double beta_j=infos[which2].beta;
			double energy_i=infos[which1].energy;
			double energy_j=infos[which2].energy;
			double power=(beta_j-beta_i)*(energy_i - energy_j);
			double ratio=exp(-power);
			if (rand<ratio)
			{
				// SWAP quantities which are being saved
				std::vector<double> tempv;
				double temp_energy, temp_mx, temp_my, temp_mz;
				tempv=infos[which1].configx;	
				infos[which1].configx=infos[which2].configx;	
				infos[which2].configx=tempv;	
				
				tempv=infos[which1].configy;	
				infos[which1].configy=infos[which2].configy;	
				infos[which2].configy=tempv;	

				tempv=infos[which1].configz;	
				infos[which1].configz=infos[which2].configz;	
				infos[which2].configz=tempv;	
				
				temp_energy=infos[which1].energy;	
				infos[which1].energy=infos[which2].energy;	
				infos[which2].energy=temp_energy;	
				
				temp_mx=infos[which1].mx;	
				infos[which1].mx=infos[which2].mx;	
				infos[which2].mx=temp_mx;	
				
				temp_my=infos[which1].my;	
				infos[which1].my=infos[which2].my;	
				infos[which2].my=temp_my;	
				
				temp_mz=infos[which1].mz;	
				infos[which1].mz=infos[which2].mz;	
				infos[which2].mz=temp_mz;
				if (which1==0 or which2==0) nswaps+=1.0;	
				//nswaps+=1.0;	
			}
		}
	}
	if (n%nsites==0 and n>nburn) infos[0].update_totals();
	if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}
	cout<<"========================================================"<<endl;
	cout<<endl;
	accept=accept/(accept+reject);
	nswaps=nswaps/nreplicatries;
	infos[0].average();
	cout<<"elowest    	= "<<boost::format("%+ .15f") %elowest<<endl;
	cout<<"Nmeas    	= "<<boost::format("%+ .15f") %infos[0].nmeas<<endl;
	cout<<"accept    	= "<<boost::format("%+ .15f") %accept<<endl;
	cout<<"nswaps (0)    	= "<<boost::format("%+ .15f") %nswaps<<endl;
	cout<<"mxavg     	= "<<boost::format("%+ .15f") %infos[0].mxavg<<endl;
	cout<<"myavg     	= "<<boost::format("%+ .15f") %infos[0].myavg<<endl;
	cout<<"mzavg     	= "<<boost::format("%+ .15f") %infos[0].mzavg<<endl;
	cout<<"mx2avg    	= "<<boost::format("%+ .15f") %infos[0].mx2avg<<endl;
	cout<<"my2avg    	= "<<boost::format("%+ .15f") %infos[0].my2avg<<endl;
	cout<<"mz2avg    	= "<<boost::format("%+ .15f") %infos[0].mz2avg<<endl;
	cout<<"mx4avg    	= "<<boost::format("%+ .15f") %infos[0].mx4avg<<endl;
	cout<<"my4avg    	= "<<boost::format("%+ .15f") %infos[0].my4avg<<endl;
	cout<<"mz4avg    	= "<<boost::format("%+ .15f") %infos[0].mz4avg<<endl;
	cout<<"eavg      	= "<<boost::format("%+ .15f") %infos[0].eavg<<endl;
	cout<<"e2avg     	= "<<boost::format("%+ .15f") %infos[0].e2avg<<endl;
	cout<<"e4avg     	= "<<boost::format("%+ .15f") %infos[0].e4avg<<endl;
	cout<<"Cv        	= "<<boost::format("%+ .15f") %infos[0].spheat<<endl;
	cout<<"Cv/Ni(J/molK)	= "<<boost::format("%+ .15f") %infos[0].spheatpersite<<endl;
	cout<<"Cvps(no dim)	= "<<boost::format("%+ .15f") %infos[0].spheatpersitenodim<<endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mc_pyrochlore_groundstate(double spin, int L,int64_t nsamples, int64_t nburn, int wait, string start_config, 
		   string mcmove, double temp, int ntemps, double hx, double hy, double hz, 
		   double J1, double J2, double J3, double J4, double Jnnn,
		   double disorder_strength,
		   double gxy, double gz, 
		   double & eavg, 
		   double &mxavg, double &myavg, double &mzavg, double &e2avg, 
		   double &mx2avg, double &my2avg, double &mz2avg, bool &measure_corrs)
{
	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Units - Assume J in meV (millielectron volts)
        //         Assume h in T   (tesla)
        //         Convert h to meV units in total_h_energy and local_h_energy
        // Convert temperature in Kelvin to temperature in meV
	// Set 16xLxLxL pyrochlore lattice
	int nsites=16*L*L*L;
       	double kB=1.38064852*1e-23;
        double NA=6.02214179*1e23;
        double JpermeV=1.60218*1e-22; 
	double tempKelvin=temp;
        temp=tempKelvin*0.08621738;

	/////////////////////////////////////////////////////////////////////////
	// J and g matrices
	
	RMatrix Jmat01,Jmat02,Jmat03,Jmat12,Jmat13,Jmat23;
	RMatrix Jmat10,Jmat20,Jmat30,Jmat21,Jmat31,Jmat32;
	RMatrix gmat0,gmat1,gmat2,gmat3; 
        RMatrix bond_disorder_matrix(nsites,nsites);
	// Given J1, J2, J3, J4 - make J mats 
	make_J_mats(J1,J2,J3,J4,Jmat01,Jmat10, Jmat02, Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, Jmat13, Jmat31, Jmat23, Jmat32);
	// Given g's - make g mats 
	make_g_mats(gxy,gz,gmat0,gmat1,gmat2,gmat3);
	std::vector<RMatrix> gmats;
	gmats.push_back(gmat0);gmats.push_back(gmat1);gmats.push_back(gmat2);gmats.push_back(gmat3);

	/////////////////////////////////////////////////////////////////////////
	// MC related quantities
	// Set nburn
	nburn=1;
	cout<<"Nsamples = "<<nsamples<<endl;
	cout<<"Nburn    = "<<nburn<<endl;
	cout<<"Wait     = "<<wait<<endl;
	double e4avg,mx4avg,my4avg,mz4avg;
	
	/////////////////////////////////////////////////////////////////////////
	// Set 16xLxLxL pyrochlore lattice
	std::vector< std::vector<int> > neighbors(nsites);	
	std::vector< std::vector<int> > nneighbors(nsites);	
	std::vector< std::vector<double> > fullcoords;	
	std::vector< std::vector<int> > ijkt;	
        make_pyrochlore(L,fullcoords,ijkt,neighbors,nneighbors,measure_corrs);
	// Make disorder matrix to be added to the Heisenberg couplings 
	// Till here the seed will give the same sequence of random numbers
        make_bond_disorder_matrix(disorder_strength,neighbors,bond_disorder_matrix);
	
	// Now we EFFECTIVELY want a different seed. So we introduce a wait parameter
	// which skips wait random numbers
	//
	//
	for (int n=0;n<wait;n++) double r=uniform_rnd();
	/////////////////////////////////////////////////////////////////////////
	// Make random configuration of spins or selected type 
	std::vector<QMC_Info> infos;
	cout<<"Ntemps   = "<<ntemps<<endl;
	double exponent=pow(100.0,1.0/double(ntemps)); // Tmax/Tmin=30
	cout<<"exponent = "<<exponent<<endl;
	for (int b=0;b<ntemps;b++)
	{
		QMC_Info qmc;
		qmc.init(L,nsites,tempKelvin*pow(exponent,double(b)),measure_corrs);
		infos.push_back(qmc);
		cout<<"Beta = "<<infos[b].beta<<endl;
		if (start_config=="random") make_random_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="111")    make_111_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz); 
		if (start_config=="x")      make_x_config(nsites,infos[b].configx,infos[b].configy,infos[b].configz);
	}
	//Matrix phases;
	//make_phases(fullcoords,infos[0].qvals,phases);

	for (int b=0;b<ntemps;b++)
	{
		//if (measure_corrs) correlations(spin,configx,configy,configz,sxsxcorrs,sysycorrs,szszcorrs,sxsycorrs,sxszcorrs,syszcorrs);
		double 		    energyj=total_j_energy(spin, infos[b].configx,infos[b].configy,infos[b].configz,
								 neighbors,nneighbors,
							   	 Jmat01, Jmat10, 
							   	 Jmat02, Jmat20, 
							   	 Jmat03, Jmat30, 
							   	 Jmat12, Jmat21, 
							   	 Jmat13, Jmat31, 
							   	 Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt);
		double 		    energyh=total_h_energy(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
							   hx,hy,hz, gmats , ijkt);
		std::vector<double> magnetization=total_magnetization(spin, infos[b].configx, infos[b].configy, infos[b].configz, 
								      gmats, ijkt);
		infos[b].energy=energyj+energyh;
		infos[b].mx=magnetization[0];infos[b].my=magnetization[1];infos[b].mz=magnetization[2];
		cout<<"Total J energy                  ="<<energyj<<endl;
		cout<<"Total h energy                  ="<<energyh<<endl;
		cout<<"Total J energy + total h energy ="<<infos[b].energy<<endl;
		cout<<"Total magnetization is          ="<<endl; print_vec_acc(magnetization,true);
	}

	cout<<endl;	
	cout<<endl;	
	cout<<"========================================================"<<endl;
	cout<<"Energy history "<<endl;	
	/////////////////////////////////////////////////////////////////////////
	// Accept reject Metropolis
	double accept=0.0;
	double reject=0.0;
	double nswaps=0.0;
	double nreplicatries=0.0;
	double elowest=0.0;
	
	for (int64_t n=0; n<(nsamples+(nburn));n++)
	{
	   if (n%nsites!=0 or ntemps==1) // Do usual Metropolis MC
	   {
			///////////////////////////////////////////////////////////////////////
			// Usual Moves of a serial Metropolis Monte Carlo
			///////////////////////////////////////////////////////////////////////
			std::vector<double> rnd1,rnd2,rnd3,rnd4;
			std::vector<int> rndints;
		        // random numbers generated in advance
			for (int b=0;b<ntemps;b++)
			{
				bool cond=false;
				double r1,r2,d1,d2;
				while (cond==false) // Problem if fixed random numbers given !!!!!!!!!
				{
					r1=uniform_rnd();
					r2=uniform_rnd();
					d1=(2.0*r1 - 1);	
					d2=(2.0*r2 - 1);	
					if (d1*d1 + d2*d2 <=1.0) cond=true; 
				}
				rnd1.push_back(r1);
				rnd2.push_back(r2);
				// First 2 rnds are drawn in a circle for the conical move to work
				rnd3.push_back(uniform_rnd());
				rnd4.push_back(uniform_rnd());
				rndints.push_back(uniform_rand_int(0,nsites));
			}
			# pragma omp parallel for
			for (int b=0;b<ntemps;b++)
			{
					// Very Small moves needed at low temperatures to increase acceptance rates
					double move_size=min(0.3,0.1*(infos[0].beta)/(infos[b].beta));  
					// Choose random site
					int site=rndints[b];
					int t=ijkt[site][3];
					double sxnew,synew,sznew;
					// Current sx,sy,sz on chosen site
					double sx=infos[b].configx[site];double sy=infos[b].configy[site];double sz=infos[b].configz[site];						// Choose a completely random direction - This is INEFFICENT at low temps
					//if (mcmove=="random")  random_move_continuous_spin(sxnew,synew,sznew);
					// Choose a completely random direction within a cone
					if (mcmove=="conical") conical_move_continuous_spin_rnds_provided(move_size,rnd1[b],rnd2[b],sx,sy,sz,sxnew,synew,sznew);
					//if (mcmove=="infDspecial")    infD_move_special_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					/*if (mcmove=="largeD")
					{
						double tempnum=rnd1[b];
						if (tempnum>0.5) conical_move_continuous_spin(move_size,sx,sy,sz,sxnew,synew,sznew);
						else	         big_move_continuous_spin(sx,sy,sz,sxnew,synew,sznew);
					}*/	
					// Normalize new direction
					double norm=sqrt(sxnew*sxnew + synew*synew + sznew*sznew); 
					sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;

					// Calculate local energy of old and new configs
					double eff_field_x, eff_field_y, eff_field_z;
					local_j_energy(spin, site,t,sx,sy,sz,
						       infos[b].configx,infos[b].configy,infos[b].configz, 
						       neighbors, nneighbors,
						       Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
						       Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
						       eff_field_x, eff_field_y,eff_field_z);
					double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
					double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
					double mxdiff=(sxnew-sx)*spin;	
					double mydiff=(synew-sy)*spin;
					double mzdiff=(sznew-sz)*spin;
							// Jterms                        // hterms - no field for now
					double ediff=(local_energyj2-local_energyj1); //+ (local_energyh2-local_energyh1);
					double beta; 
					beta=infos[b].beta;
					double prob=exp(-beta*ediff);
					double rand=rnd3[b]; // random number previously generated
					if (rand<prob) // Metropolis Accept-reject for a given temperature
					{
						if (b==0 and n>(nburn)) accept=accept+1.0;
						// Reset configs to new configs, because accepted
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].energy+=ediff;
						//if (infos[b].energy<elowest) elowest=infos[b].energy;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
					}
					else
					{
					       	if (b==0 and n>(nburn)) reject=reject+1.0;
						// Still move spin if rejected, but this conserves energy!!
						double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
						norm=sqrt(norm);
						eff_field_x=eff_field_x/norm;
						eff_field_y=eff_field_y/norm;
						eff_field_z=eff_field_z/norm;
						double costheta=(sx*eff_field_x)+(sy*eff_field_y)+(sz*eff_field_z);
						double sintheta=sqrt(1.0-(costheta*costheta));	
						double phi=rnd4[b]*2.0*3.14159; // random in range 0 to 2 pi 
										// random number previously generated
						double x1,y1,z1, x2,y2,z2;	
						get_two_normalized_orth_dirs(eff_field_x, eff_field_y, eff_field_z, 
		       							     x1, y1, z1, 
		       							     x2, y2, z2);
						double sxnew=(costheta*eff_field_x)  + (sintheta*cos(phi)*x1)+ (sintheta*sin(phi)*x2);
						double synew=(costheta*eff_field_y)  + (sintheta*cos(phi)*y1)+ (sintheta*sin(phi)*y2);
						double sznew=(costheta*eff_field_z)  + (sintheta*cos(phi)*z1)+ (sintheta*sin(phi)*z2);
						double mxdiff=(sxnew-sx)*spin;	
						double mydiff=(synew-sy)*spin;
						double mzdiff=(sznew-sz)*spin;
						infos[b].configx[site]=sxnew;infos[b].configy[site]=synew;infos[b].configz[site]=sznew;
						infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
						/*for (int i=0;i<infos[b].numqs;i++)
						{
							infos[b].sxq[i]+=(mxdiff*phases(i,site));
							infos[b].syq[i]+=(mydiff*phases(i,site));
							infos[b].szq[i]+=(mzdiff*phases(i,site));
						}*/
						// ediff=0 by construction!!
					}
			}
	}	
	else if (ntemps>1)  //50 % chance of doing replica exchange
	{
		nreplicatries+=1.0;	
		///////////////////////////////////////////////////////////////////////
		////////// Attempted exchange moves of parallel tempering
		///////////////////////////////////////////////////////////////////////
		// Metropolis move done, try swapping every 2 sweeps
		for (int which1=0;which1<ntemps;which1++)
		{
			int which2=(which1+1)%(ntemps);
			// Slight bias, if i try to swap serially ?
			double rand=uniform_rnd();	
			double beta_i=infos[which1].beta;
			double beta_j=infos[which2].beta;
			double energy_i=infos[which1].energy;
			double energy_j=infos[which2].energy;
			double power=(beta_j-beta_i)*(energy_i - energy_j);
			double ratio=exp(-power);
			if (rand<ratio)
			{
				// SWAP quantities which are being saved
				std::vector<double> temp;
				double temp_energy, temp_mx, temp_my, temp_mz;
				temp=infos[which1].configx;	
				infos[which1].configx=infos[which2].configx;	
				infos[which2].configx=temp;	
				
				temp=infos[which1].configy;	
				infos[which1].configy=infos[which2].configy;	
				infos[which2].configy=temp;	

				temp=infos[which1].configz;	
				infos[which1].configz=infos[which2].configz;	
				infos[which2].configz=temp;	
				
				temp_energy=infos[which1].energy;	
				infos[which1].energy=infos[which2].energy;	
				infos[which2].energy=temp_energy;	
				
				temp_mx=infos[which1].mx;	
				infos[which1].mx=infos[which2].mx;	
				infos[which2].mx=temp_mx;	
				
				temp_my=infos[which1].my;	
				infos[which1].my=infos[which2].my;	
				infos[which2].my=temp_my;	
				
				temp_mz=infos[which1].mz;	
				infos[which1].mz=infos[which2].mz;	
				infos[which2].mz=temp_mz;
				if (which1==0 or which2==0) nswaps+=1.0;	
				//nswaps+=1.0;	
			}
		}
	}
	if (n%nsites==0 or n==nburn+nsamples-1) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}
	cout<<"========================================================"<<endl;
	cout<<endl;
	accept=accept/(accept+reject);
	nswaps=nswaps/nreplicatries;

	///////////////////////////////////////////////////////////////////////////////////////////
	cout<<"Started iterative minimization on lowest configurations in parallel"<<endl;
	for (int64_t n=0; n<(nsamples+(nburn));n++)
	{
		//cout<<"Started iterative minimzation for configuration number "<<b<<" at beta = "<<infos[b].beta<<endl;
		// random numbers generated in advance
		std::vector<double> rnd1;
		std::vector<int> rndints;
		for (int b=0;b<ntemps;b++)
		{
			rnd1.push_back(uniform_rnd());
			rndints.push_back(uniform_rand_int(0,nsites));
		}
		# pragma omp parallel for
		for (int b=0;b<min(100,ntemps);b++)
		{
			int site=rndints[b];
			double sx=infos[b].configx[site];
			double sy=infos[b].configy[site];
			double sz=infos[b].configz[site];
			int t=ijkt[site][3];
			// Calculate local energy of old and new configs
			double eff_field_x, eff_field_y, eff_field_z;
			local_j_energy(spin, site,t,sx,sy,sz,
			infos[b].configx,infos[b].configy,infos[b].configz, 
			neighbors, nneighbors,
			Jmat01,Jmat10, Jmat02,Jmat20, Jmat03, Jmat30, Jmat12, Jmat21, 
			Jmat13,Jmat31, Jmat23, Jmat32, Jnnn, bond_disorder_matrix, ijkt,
			eff_field_x, eff_field_y,eff_field_z);
			double norm=(eff_field_x*eff_field_x)+(eff_field_y*eff_field_y)+(eff_field_z*eff_field_z);
			norm=sqrt(norm);
			double r=0.5*rnd1[b];
			double sxnew=((-eff_field_x/norm))*r + sx ;
			double synew=((-eff_field_y/norm))*r + sy ;
			double sznew=((-eff_field_z/norm))*r + sz ;
			norm=(sxnew*sxnew) + (synew*synew) + (sznew*sznew);
			norm=sqrt(norm);
			sxnew=sxnew/norm; synew=synew/norm; sznew=sznew/norm;
			double local_energyj1=( (sx*eff_field_x) + (sy*eff_field_y)+ (sz*eff_field_z))*spin*spin; 
			double local_energyj2=( (sxnew*eff_field_x) + (synew*eff_field_y)+ (sznew*eff_field_z))*spin*spin; 
			infos[b].energy+=(local_energyj2-local_energyj1);
			infos[b].configx[site]=sxnew;
			infos[b].configy[site]=synew;
			infos[b].configz[site]=sznew;
			double mxdiff=(sxnew-sx)*spin;	
			double mydiff=(synew-sy)*spin;
			double mzdiff=(sznew-sz)*spin;
			infos[b].mx+=mxdiff;infos[b].my+=mydiff;infos[b].mz+=mzdiff;
		}
		if ((n%nsites==0 or n==nburn+nsamples-1) ) cout<<boost::format("%+ .5f") %infos[0].energy<<endl;
	}

	for (int b=0;b<min(100,ntemps);b++)
	{	
		cout<<"========================================================================================================================="<<endl;
		cout<<"For configuration number "<<b<<" at beta = "<<infos[b].beta<<endl;
		cout<<"========================================================================================================================="<<endl;
		///////////////////////////////////////////////////////////////////////////////////////////
		fourier_transforms_slow(spin,fullcoords,infos[b].qvals,infos[b].configx,infos[b].configy,infos[b].configz,infos[b].sxq,infos[b].syq,infos[b].szq);
		infos[b].update_totals(); // Update totals - only last snapshot*/
		// Only the lowest temperature is relevant for averages we are interested in 
		infos[b].average();

		std::vector<double> configx=infos[b].configx;
		std::vector<double> configy=infos[b].configy;
		std::vector<double> configz=infos[b].configz;

		double invsqrt12=1.0/sqrt(12.0);
		cout<<"========================================================================================================================="<<endl;
		cout<<"    x            y             z             sx             sy             sz                                            "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<nsites;i++) // Loop over all sites
		{
			cout<<boost::format("%+ .10f  %+ .10f  %+ .10f  %+ .10f   %+ .10f   %+ .10f") %fullcoords[i][0] %fullcoords[i][1] %fullcoords[i][2] %(configx[i]*spin) %(configy[i]*spin) %(configz[i]*spin) <<endl;

		}
		cout<<endl;
		cout<<endl;

		cout<<"========================================================================================================================="<<endl;
		cout<<"    s01            s02             s03             s12             s13             s23            f0             f1          f2"<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<nsites;i+=4) // Loop over all up tetrahedra
		{
			double s00=(configx[i]*configx[i])+(configy[i]*configy[i])+(configz[i]*configz[i]);
			double s01=(configx[i]*configx[i+1])+(configy[i]*configy[i+1])+(configz[i]*configz[i+1]);
			double s02=(configx[i]*configx[i+2])+(configy[i]*configy[i+2])+(configz[i]*configz[i+2]);
			double s03=(configx[i]*configx[i+3])+(configy[i]*configy[i+3])+(configz[i]*configz[i+3]);
			double s11=(configx[i+1]*configx[i+1])+(configy[i+1]*configy[i+1])+(configz[i+1]*configz[i+1]);
			double s12=(configx[i+1]*configx[i+2])+(configy[i+1]*configy[i+2])+(configz[i+1]*configz[i+2]);
			double s13=(configx[i+1]*configx[i+3])+(configy[i+1]*configy[i+3])+(configz[i+1]*configz[i+3]);
			double s22=(configx[i+2]*configx[i+2])+(configy[i+2]*configy[i+2])+(configz[i+2]*configz[i+2]);
			double s23=(configx[i+2]*configx[i+3])+(configy[i+2]*configy[i+3])+(configz[i+2]*configz[i+3]);
			double s33=(configx[i+3]*configx[i+3])+(configy[i+3]*configy[i+3])+(configz[i+3]*configz[i+3]);

			double f0=((s00+s11+s22+s33)+(2.0*(s01+s02+s03+s12+s13+s23)))*spin*spin;
			double f1=(s03+s12+s13+s02-(2.0*s01)-(2.0*s23))*invsqrt12*spin*spin;
			double f2=((s13+s02-s03-s12)/2.0)*spin*spin;
			
			cout<<boost::format("%+ .10f  %+ .10f  %+ .10f  %+ .10f   %+ .10f   %+ .10f   %+ .10f   %+ .10f   %+ .10f") %s01 %s02 %s03 %s12 %s13 %s23 %f0 %f1 %f2<<endl;

		}
		cout<<endl;
		cout<<endl;

		//cout<<"elowest    	= "<<boost::format("%+ .15f") %elowest<<endl;
		cout<<"Nmeas    	= "<<boost::format("%+ .15f") %infos[b].nmeas<<endl;
		cout<<"accept    	= "<<boost::format("%+ .15f") %accept<<endl;
		cout<<"nswaps (0)    	= "<<boost::format("%+ .15f") %nswaps<<endl;
		cout<<"mxavg     	= "<<boost::format("%+ .15f") %infos[b].mxavg<<endl;
		cout<<"myavg     	= "<<boost::format("%+ .15f") %infos[b].myavg<<endl;
		cout<<"mzavg     	= "<<boost::format("%+ .15f") %infos[b].mzavg<<endl;
		cout<<"mx2avg    	= "<<boost::format("%+ .15f") %infos[b].mx2avg<<endl;
		cout<<"my2avg    	= "<<boost::format("%+ .15f") %infos[b].my2avg<<endl;
		cout<<"mz2avg    	= "<<boost::format("%+ .15f") %infos[b].mz2avg<<endl;
		cout<<"mx4avg    	= "<<boost::format("%+ .15f") %infos[b].mx4avg<<endl;
		cout<<"my4avg    	= "<<boost::format("%+ .15f") %infos[b].my4avg<<endl;
		cout<<"mz4avg    	= "<<boost::format("%+ .15f") %infos[b].mz4avg<<endl;
		cout<<"eavg      	= "<<boost::format("%+ .15f") %infos[b].eavg<<endl;
		cout<<"e2avg     	= "<<boost::format("%+ .15f") %infos[b].e2avg<<endl;
		cout<<"e4avg     	= "<<boost::format("%+ .15f") %infos[b].e4avg<<endl;
		cout<<"Cv        	= "<<boost::format("%+ .15f") %infos[b].spheat<<endl;
		cout<<"Cv/Ni(J/molK)	= "<<boost::format("%+ .15f") %infos[b].spheatpersite<<endl;
		cout<<"Cvps(no dim)	= "<<boost::format("%+ .15f") %infos[b].spheatpersitenodim<<endl;
		
		cout<<"========================================================================================================================="<<endl;
		cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)         "<<endl;
		cout<<"========================================================================================================================="<<endl;
		for (int i=0;i<infos[b].numqs;i++)
		{
		cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f  %+ .5f  %+ .5f %+ .5f") %infos[b].qvals[i][0] %infos[b].qvals[i][1] %infos[b].qvals[i][2] %infos[b].sxsxtot[i] %infos[b].sxsytot[i] %infos[b].sxsztot[i] %infos[b].sysxtot[i] %infos[b].sysytot[i] %infos[b].sysztot[i] %infos[b].szsxtot[i] %infos[b].szsytot[i] %infos[b].szsztot[i]<<endl;
			
		}
		cout<<endl;
		cout<<endl;
	}

	
	// AVERAGED DATA	
	int numqs=infos[0].qvals.size();
	std::vector<complex<double> > allsxsxtot(numqs);
	std::vector<complex<double> > allsxsytot(numqs);
	std::vector<complex<double> > allsxsztot(numqs);
	
	std::vector<complex<double> > allsysxtot(numqs);
	std::vector<complex<double> > allsysytot(numqs);
	std::vector<complex<double> > allsysztot(numqs);
	
	std::vector<complex<double> > allszsxtot(numqs);
	std::vector<complex<double> > allszsytot(numqs);
	std::vector<complex<double> > allszsztot(numqs);
	
	for (int b=0;b<min(100,ntemps);b++)
	{
		#pragma omp parallel for
		for (int nq=0;nq<infos[b].qvals.size();nq++)
		{
			allsxsxtot[nq]+=infos[b].sxsxtot[nq];
			allsxsytot[nq]+=infos[b].sxsytot[nq];
			allsxsztot[nq]+=infos[b].sxsztot[nq];
		
			allsysxtot[nq]+=infos[b].sysxtot[nq];
			allsysytot[nq]+=infos[b].sysytot[nq];
			allsysztot[nq]+=infos[b].sysztot[nq];
			
			allszsxtot[nq]+=infos[b].szsxtot[nq];
			allszsytot[nq]+=infos[b].szsytot[nq];
			allszsztot[nq]+=infos[b].szsztot[nq];
		}
	}

	int nconfigs=min(100,ntemps);
	
	#pragma omp parallel for
	for (int nq=0;nq<numqs;nq++)
	{
		allsxsxtot[nq]=allsxsxtot[nq]/double(nconfigs);
		allsxsytot[nq]=allsxsytot[nq]/double(nconfigs);
		allsxsztot[nq]=allsxsztot[nq]/double(nconfigs);
	
		allsysxtot[nq]=allsysxtot[nq]/double(nconfigs);
		allsysytot[nq]=allsysytot[nq]/double(nconfigs);
		allsysztot[nq]=allsysztot[nq]/double(nconfigs);
		
		allszsxtot[nq]=allszsxtot[nq]/double(nconfigs);
		allszsytot[nq]=allszsytot[nq]/double(nconfigs);
		allszsztot[nq]=allszsztot[nq]/double(nconfigs);
	}
	
	cout<<"========================================================================================================================="<<endl;
	cout<<"                                  AVERAGED  DATA                                                                         "<<endl;
	cout<<"========================================================================================================================="<<endl;
	cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)         "<<endl;
	cout<<"========================================================================================================================="<<endl;
	for (int i=0;i<numqs;i++)
	{
		cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f   %+ .5f  %+ .5f  %+ .5f %+ .5f") %infos[0].qvals[i][0] %infos[0].qvals[i][1] %infos[0].qvals[i][2] %allsxsxtot[i] %allsxsytot[i] %allsxsztot[i] %allsysxtot[i] %allsysytot[i] %allsysztot[i] %allszsxtot[i] %allszsytot[i] %allszsztot[i]<<endl;
			
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


