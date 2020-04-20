#include "math_utilities.h"
#include "mtrand.h"
#include"number_functions.h"
#include"search_for.h"

using namespace std;

double uniform_rnd()
{   MTRand drand;
    return drand();
}

int uniform_rand_int(int lowest,int highest)
{
	double rand=uniform_rnd();
	return lowest+int((highest-lowest)*rand);
}
////////////////////////////////////////////////////////////////////////////
double mean(double sum_x, int nsamples)
{   return (sum_x)/double(nsamples);}

double variance(double sum_x, double sum_x_2, int nsamples)
{   double var;
    var=sum_x_2-((sum_x*sum_x)/double(nsamples));
    if (var<0) {cout<<"Some bug!"<<endl;cout<<"var="<<var<<endl;}
    return var/double(nsamples);
}
////////////////////////////////////////////////////////////////////////////

double sigma(double sum_x, double sum_x_2, int nsamples)
{return sqrt(abs(variance(sum_x,sum_x_2,nsamples)));}

////////////////////////////////////////////////////////////////////////////

double sigma_over_root_n(double sum_x, double sum_x_2, int nsamples)
{return sqrt(variance(sum_x,sum_x_2,nsamples))/sqrt(double(nsamples));}

////////////////////////////////////////////////////////////////////////////

double error_in_mean(double sum_x, double sum_x_2, int nsamples)
// error in mean is sigma over root n
{return sigma_over_root_n(sum_x, sum_x_2, nsamples);}

////////////////////////////////////////////////////////////////////////////

int closest_int(double nearint)
{
   if (nearint>0){return int(nearint+0.5);}
   else{return int(nearint-0.5);}	
}

////////////////////////////////////////////////////////////////////////////

void convert_to_g6(int nsites,
                   std::vector< std::vector<int> > &ordered_pairs, 
                   std:: string &tmp_string)
{
    std::vector<int> all_bits,bit_rep,new_bit_rep;
    std::vector<int> word(6);
    int i,j,p,excess_needed,nwords,num=0;
    char c='a';
    stringstream ss;
    string s;

    tmp_string=string("");
    c=nsites+63;

    ss<<c;ss>>s;tmp_string+=s;

    all_bits.clear();
    for (p=1;p<nsites;p++)
    {
        bit_rep.push_back(0);new_bit_rep=bit_rep;
        for (i=0;i<ordered_pairs.size();i++)
        {
            if (ordered_pairs[i][1]==p)
            {new_bit_rep[ordered_pairs[i][0]]=1;}
        }
        all_bits.insert(all_bits.end(),new_bit_rep.begin(),new_bit_rep.end());
    }
        
    if (nsites*(nsites-1)/2 % 6 ==0){nwords=(nsites*(nsites-1)/12);}
    else
    {
        nwords=(nsites*(nsites-1)/12) + 1;
        excess_needed=6-(nsites*(nsites-1)/2 %6);
        for (j=0;j<excess_needed;j++) {all_bits.push_back(0);}
    }
    
    for (i=0;i<nwords;i++)
    {
        stringstream ss2;
        word.clear();
        for (j=0;j<6;j++){word.push_back(all_bits[(6*i)+j]);}
        
        num=convert_vec_to_num(word,2)+63; c=num;
        ss2<<c;ss2>>s;tmp_string+=s;
    }

    return;
}
////////////////////////////////////////////////////////////////////////////

void convert_all_to_g6(std::string read_file)
{
    bool found;
    int i,j,ctr,nsites,max;
    char *filechar;
    std::string str,str_ret,tmp_string,new_file;
    std::vector< std::vector<int> > pairs;

    ctr=0;found=true;
    new_file=read_file+string(".g6");
    filechar=new char [new_file.size()+1];
    strcpy (filechar,new_file.c_str());
    ofstream clus_file(filechar);

    while (found)
    {
        str="pairs_cluster_";
        str+=to_string(ctr);
        search_for(str,read_file,str_ret,found);
        ctr+=1;

        if (found)
        {
              if (str_ret.substr(0,1)==string("[")){pairs=convert_string_to_vec_of_vec(str_ret);}
              
              max=0;
              for (i=0;i<pairs.size();i++)
              { for (j=0;j<2;j++){if (pairs[i][j]>max) {max=pairs[i][j];} } }
              nsites=max+1; 
              
	      convert_to_g6(nsites,pairs,tmp_string);

              if (clus_file.is_open())
              {clus_file<<tmp_string;clus_file<<string("\n");}  
              else
              {cout<<"I could not dump g6 graphs! "<<endl;}
        }
    }
    delete[] filechar;
}

//////////////////////////////////////////////////////////////////////////////



