#include"printing_functions.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// Aesthetic/ Printing functions - programmed by HJC on Nov 7 2009
//////////////////////////////////////////////////////////////////////////////
// Print a vector with integer entries

template<class T>
void print_vec(T &vec,bool vert)
{
    if (not vert)
    {		
    	for (int i=0;i<vec.size();i++){cout<<vec[i]<<"  ";}
    	cout<<endl;
    }
    else
    {   for (int i=0;i<vec.size();i++){cout<<vec[i]<<endl;} }
}

void print_vec_acc(std::vector<double> &vec,bool vert,int size)
{
    if (size==0) {size=vec.size();}
    
    size=min(size,int(vec.size()));

    if (not vert)
    {		
    	for (int i=0;i<size;i++){cout<<boost::format("%+.10f") %vec[i]<<"  ";}
    	cout<<endl;
    }
    else
    {   for (int i=0;i<size;i++){cout<<boost::format("%+.10f") %vec[i]<<endl;} }
}

template void print_vec< std::vector<double> > (std::vector<double> &, bool);
template void print_vec< std::vector<int> > (std::vector<int> &, bool);
template void print_vec< std::vector<complex<double> >  > (std::vector<complex<double> > &, bool);

// Print a vector of vector with integer entries
void print_mat_int(  std::vector< std::vector<int> > const &mat)
{ print_mat(mat);}

// Print a vector of vector with double entries
void print_mat_double(std::vector< std::vector<double> > const &mat)
{print_mat(mat);}

template<class T>
void print_mat( T &mat)
{
    for (int i=0;i<mat.size();i++)
    { 
        for (int j=0;j<mat[i].size();j++)
        {cout<<mat[i][j]<<"  ";}
        cout<<endl;
     }
}


void print_mathematica_pairs(std::vector< std::vector<int> > const &pairs)
{
	cout<<"{";
	for (int i=0;i<pairs.size()-1;i++)
	{cout<<pairs[i][0]<<"->"<<pairs[i][1]<<",";}
	cout<<pairs[pairs.size()-1][0]<<"->"<<pairs[pairs.size()-1][1]<<"}"<<endl;
}

void print_mathematica_vector(std::vector<double> const &vec)
{
	cout<<"{";
	for (int i=0;i<vec.size()-1;i++){cout<<vec[i]<<",";}
	cout<<vec[vec.size()-1]<<"}"<<endl;
}

void print_real_mat(RMatrix &mat)
{
    int i,j;
    for (i=0;i<mat.NRows();i++)
    {
        for (j=0;j<mat.NCols();j++) cout<<boost::format("%+.3f") % mat(i,j)<<" ";
        cout<<endl;
    }
    cout<<endl;
}



