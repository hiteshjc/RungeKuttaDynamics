#include"global.h"
#include"mtrand.h"

// Search
#include"search_for.h"

// Hamiltonian
#include"hamiltonian.h"
#include"oleg.h"
// MC
#include"mc.h"
#include"mc_finite_D.h"
#include"mc_pyrochlore.h"
#include"number_functions.h"
using namespace std;


int main(int argc, char *argv[])
{
    time_t start,end;
    double dif;
    bool found;
    int seed;
    MTRand irand;
    string str_ret,neigs_str,nkrylov_str,filename;
    if (argc <= 1)
    {
        cout << "Usage: " << argv[0] << " <Filename>" << endl;
        exit(1);
    }
 
    filename=argv[1];
    search_for(string("seed"),filename,str_ret,found);
    if (found) 
    {
		if (found) {seed=str_to_int(str_ret);}
		else {seed=1;}
    }	
    irand.seed(seed);
    double e,mx,my,mz,e2,mx2,my2,mz2;
    bool Lfound,nstatesfound,measurecorrsfound,disorderfound, nreplicasfound, waitfound; 
    bool spinfound, tempfound,Dfound, Dpfound, mcmovefound,configfound, nburnfound, nsamplesfound, latticefound, Jnnnfound; 
    bool hxfound,hyfound,hzfound, J1found, J2found, J3found, J4found, Jvalfound, gxyfound, gzfound, alphaLfound, kcutfound;
    int L,nstates=2;
    double spin, Jval, J1, J2, J3, J4, Jnnn, gxy, gz, hx,hy,hz,temp,D,Dp,alphaL;
    string mcmove,start_config,lattice;
    double nsamplesd, nburnd; 
    int nreplicas; 
    int kcut;
    int wait;
    bool measure_corrs=true;
    double disorder_strength=0.0;
   
    search_for(string("nstates"),filename,str_ret,nstatesfound);
    if (nstatesfound){nstates=str_to_int(str_ret);} else{nstates=2;}
    
    search_for(string("spin"),filename,str_ret,spinfound);
    if (spinfound){spin=str_to_d(str_ret);} else{spin=0.5;}
    
    search_for(string("L"),filename,str_ret,Lfound);
    if (Lfound){L=str_to_int(str_ret);} else{L=4;}
    
    search_for(string("nsamples"),filename,str_ret,nsamplesfound);
    if (nsamplesfound){nsamplesd=str_to_d(str_ret);} else{nsamplesd=1000000000;}
    
    search_for(string("nburn"),filename,str_ret,nburnfound);
    if (nburnfound){nburnd=str_to_d(str_ret);} else{nburnd=100000000;}
    int64_t nsamples=int64_t(nsamplesd);
    int64_t nburn=int64_t(nburnd);
    
    search_for(string("nreplicas"),filename,str_ret,nreplicasfound);
    if (nreplicasfound){nreplicas=str_to_int(str_ret);} else{nreplicas=1;}
    
    search_for(string("Jval"),filename,str_ret,Jvalfound);
    if (Jvalfound){Jval=str_to_d(str_ret);} else{Jval=0.0;}
    
    search_for(string("J1"),filename,str_ret,J1found);
    if (J1found){J1=str_to_d(str_ret);} else{J1=0;}
    
    search_for(string("J2"),filename,str_ret,J2found);
    if (J2found){J2=str_to_d(str_ret);} else{J2=0;}
    
    search_for(string("J3"),filename,str_ret,J3found);
    if (J3found){J3=str_to_d(str_ret);} else{J3=0;}
    
    search_for(string("J4"),filename,str_ret,J4found);
    if (J4found){J4=str_to_d(str_ret);} else{J4=0;}
    
    search_for(string("Jnnn"),filename,str_ret,Jnnnfound);
    if (Jnnnfound){Jnnn=str_to_d(str_ret);} else{Jnnn=0;}
    
    search_for(string("disorder"),filename,str_ret,disorderfound);
    if (disorderfound){disorder_strength=str_to_d(str_ret);} else{disorder_strength=0.0;}
    
    search_for(string("gxy"),filename,str_ret,gxyfound);
    if (gxyfound){gxy=str_to_d(str_ret);} else{gxy=0;}
    
    search_for(string("gz"),filename,str_ret,gzfound);
    if (gzfound){gz=str_to_d(str_ret);} else{gz=0;}
    
    search_for(string("hx"),filename,str_ret,hxfound);
    if (hxfound){hx=str_to_d(str_ret);} else{hx=0;}
    
    search_for(string("hy"),filename,str_ret,hyfound);
    if (hyfound){hy=str_to_d(str_ret);} else{hy=0;}
    
    search_for(string("hz"),filename,str_ret,hzfound);
    if (hzfound){hz=str_to_d(str_ret);} else{hz=0;}
    
    search_for(string("D"),filename,str_ret,Dfound);
    if (Dfound){D=str_to_d(str_ret);} else{D=0.0;}
    
    search_for(string("Dp"),filename,str_ret,Dpfound);
    if (Dpfound){Dp=str_to_d(str_ret);} else{Dp=0.0;}
    
    search_for(string("alphaL"),filename,str_ret,alphaLfound);
    if (alphaLfound){alphaL=str_to_d(str_ret);} else{alphaL=10.0;}
    
    search_for(string("nwait"),filename,str_ret,waitfound);
    if (waitfound){wait=str_to_int(str_ret);} else{wait=0;}
    
    int nstarts;
    bool nstartsfound;
    search_for(string("nstarts"),filename,str_ret,nstartsfound);
    if (nstartsfound){nstarts=str_to_int(str_ret);} else{nstarts=1;}
    
    search_for(string("kcut"),filename,str_ret,kcutfound);
    if (kcutfound){kcut=str_to_int(str_ret);} else{kcut=5;}
    
    search_for(string("temp"),filename,str_ret,tempfound);
    if (tempfound){temp=str_to_d(str_ret);} else{temp=1;}
   
    double tottime;
    bool   tottimefound; 
    search_for(string("tottime"),filename,str_ret,tottimefound);
    if (tottimefound){tottime=str_to_d(str_ret);} else{tottime=50;}
    
    bool   deltatfound;
    double deltat;
    search_for(string("deltat"),filename,str_ret,deltatfound);
    if (deltatfound){deltat=str_to_d(str_ret);} else{deltat=0.01;}
    
    bool   omegafound;
    double omega;
    search_for(string("omega"),filename,str_ret, omegafound);
    if (omegafound){omega=str_to_d(str_ret);} else{omega=0.00;}
   
   
    bool tminfound,tmaxfound; 
    double tmin, tmax;
 
    search_for(string("tmin"),filename,str_ret,tminfound);
    if (tminfound){tmin=str_to_d(str_ret);} else{tmin=1;}
    
    search_for(string("tmax"),filename,str_ret,tmaxfound);
    if (tmaxfound){tmax=str_to_d(str_ret);} else{tmax=1;}
    
    search_for(string("mcmove"),filename,str_ret,mcmovefound);
    if (mcmovefound){mcmove=str_ret;} else{mcmove="random";}
    
    search_for(string("start_config"),filename,str_ret,configfound);
    if (configfound){start_config=str_ret;} else{start_config="random";}
    
    search_for(string("lattice"),filename,str_ret,latticefound);
    if (latticefound){lattice=str_ret;} else{lattice="cubic";}
    
    string method="finitetemp";
    bool methodfound;
    search_for(string("method"),filename,str_ret,methodfound);
    if (methodfound){method=str_ret;}
    
    search_for(string("measure_corrs"),filename,str_ret,measurecorrsfound);
    if (measurecorrsfound){measure_corrs=str_to_bool(str_ret);} else{measure_corrs=true;}

    cout<<"lattice   		= "<<lattice<<endl;    
    cout<<"spin      		= "<<spin<<endl;    
    cout<<"L         		= "<<L<<endl;    
    cout<<"temp (K)  		= "<<temp<<endl;    
    cout<<"tmin (K)  		= "<<tmin<<endl;    
    cout<<"tmax (K)  		= "<<tmax<<endl;    
    cout<<"gxy       		= "<<gxy<<endl;    
    cout<<"gz        		= "<<gz<<endl;    
    cout<<"J1   (meV)		= "<<J1<<endl;    
    cout<<"J2   (meV)		= "<<J2<<endl;    
    cout<<"J3   (meV)		= "<<J3<<endl;    
    cout<<"J4   (meV)		= "<<J4<<endl;    
    cout<<"Jnnn (meV)		= "<<Jnnn<<endl;    
    cout<<"disorder (meV)  	= "<<disorder_strength<<endl;    
    cout<<"hx   (T)  		= "<<hx<<endl;    
    cout<<"hy   (T)  		= "<<hy<<endl;    
    cout<<"hz   (T)  		= "<<hz<<endl;    
    cout<<"D/J       		= "<<D<<endl;    
    cout<<"nstarts   		= "<<nstarts<<endl;    
    cout<<"nreplicas   		= "<<nreplicas<<endl;    
    cout<<"mcmove    		= "<<mcmove<<endl;    
    cout<<"st_config 		= "<<start_config<<endl;    
    cout<<"measure   		= "<<measure_corrs<<endl; 
    cout<<"deltat  (meV)^-1	= "<<deltat<<endl; 
    cout<<"omega (meV)   	= "<<omega<<endl; 
    cout<<"Tot-time (meV)^-1   	= "<<tottime<<endl; 
    int numprocs=nstarts;
    STensor savg,smunu;
    for (int i=0;i<numprocs;i++)
    {
	time(&start);
	mc_pyrochlore_get_thermal_config_and_time_evolve(spin, deltat, omega, tottime, L, nsamples,nburn,start_config,mcmove,temp,nreplicas, hx,hy,hz,J1,J2,J3,J4,Jnnn,disorder_strength,gxy,gz,smunu);
	if (i>0) savg.update_totals(smunu);
	else	 { savg.init(int(smunu.qvals.size())); savg.copy(smunu);}
	time(&end);
	double seconds=difftime(end,start);
    	cout<<"Done with run "<<i<<" in  "<<seconds<<" seconds"<<endl;
    }
    savg.average();
    
    int numqs=int(savg.qvals.size());
    cout<<"======================================================================================================================================="<<endl;
    cout<<"                                              AVERAGED DATA                                                                            "<<endl;
    cout<<"======================================================================================================================================="<<endl;
    cout<<" h      k      l     SXX(Q)    SXY(Q)    SXZ(Q)    SYX(Q)     SYY(Q)     SYZ(Q)     SZX(Q)     SZY(Q)     SZZ(Q)          Sperp(Q)     "<<endl;
    cout<<"======================================================================================================================================="<<endl;
    for (int i=0;i<numqs;i++)
    {
	double qx=savg.qvals[i][0];double qy=savg.qvals[i][1];double qz=savg.qvals[i][2];
	complex<double> sxx=savg.sxx[i];complex<double> syy=savg.syy[i];complex<double> szz=savg.szz[i];
	
	complex<double> sxy=savg.sxy[i];complex<double> syx=savg.syx[i];
	complex<double> sxz=savg.sxz[i];complex<double> szx=savg.szx[i];
	complex<double> syz=savg.syz[i];complex<double> szy=savg.szy[i];
	complex<double> sperp=savg.sperp[i];
	
	cout<<boost::format("%+ .5f  %+ .5f  %+ .5f  %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f   %+ .8f  %+ .8f  %+ .8f %+ .8f  %+ .8f") %qx %qy %qz %sxx %sxy %sxz %syx %syy %syz %szx %szy %szz %sperp<<endl;
   }				

}

/////////////////////////////////////////////////////////////////////////////
