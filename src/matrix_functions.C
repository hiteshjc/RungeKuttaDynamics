#include"bethe_lapack_interface.h"
#include"matrix_functions.h"
#include"printing_functions.h"
#include"math_utilities.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////////
void symmetric_diagonalize(Matrix &a, 
                 std::vector<double> &eigs, 
                 Matrix &eigenvecs)
{
    int w,info=0;
    int n=int(a.NRows());
    std::vector< complex<double> > work(4*n);
    std::vector<double> rwork(4*n);

    //cout<<"Making copy of matrix"<<endl;
    eigenvecs=a;
    work.assign(4*n,(0.,0.));
     
    //print_real_mat(a);

    //cout<<"Actual diagonalization"<<endl;
    zheev('V','U',n,&*eigenvecs.begin(),n,
           &*eigs.begin(),&*work.begin(),int(work.size()),&*rwork.begin(),info);
    
    //cout<<"Finished diagonalizing (now returning)"<<endl;

    // Sorted the eigenvalues and corresponding eigenvectors
    // Eigenvectors are columns of the matrix named eigenvecs

}

void real_symmetric_diagonalize(RMatrix &a, 
                 std::vector<double> &eigs, 
                 RMatrix &eigenvecs)
{
    int w,info=0;
    int n=int(a.NRows());
    std::vector<double> work(4*n);

    //cout<<"Making copy of matrix"<<endl;
    eigenvecs=a;
    work.assign(4*n,0.);

    //print_real_mat(a);

    //cout<<"Actual diagonalization"<<endl;
    dsyev('V','U',n,&*eigenvecs.begin(),n,
           &*eigs.begin(),&*work.begin(),int(work.size()),info);
    
    //cout<<"Finished diagonalizing (now returning)"<<endl;

    // Sorted the eigenvalues and corresponding eigenvectors
    // Eigenvectors are columns of the matrix named eigenvecs

}

/////////////////////////////////////////////////////////////
void real_matrix_multiply( RMatrix &a, RMatrix &b, RMatrix &c)
{
    double alpha=1.0;
    double beta=0.0;
    

    if (a.NCols()!=b.NRows())
    {
	cout<<"Problem in matrix multiply"<<endl;
	return;
    }
    c.resize(a.NRows(),b.NCols());

    /*dgemm('N', 'N' , int(a.NRows()) ,int(b.NCols()) ,int(a.NCols()) ,alpha ,
          &*a1.begin(), int(a.NRows()), &*b1.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));  */
    
    dgemm('N', 'N' , int(a.NRows()) ,int(b.NCols()) ,int(a.NCols()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));
}


/////////////////////////////////////////////////////////////
void real_matrix_multiply_atb( RMatrix &a, RMatrix &b, RMatrix &c)
{
    double alpha=1.0;
    double beta=0.0;
    
    if (a.NRows()!=b.NRows())
    {
	cout<<"Problem in matrix multiply atb"<<endl;
	return;
    }
    c.resize(a.NCols(),b.NCols());

    dgemm('T', 'N' , int(a.NCols()) ,int(b.NCols()) ,int(a.NRows()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));
}

/////////////////////////////////////////////////////////////
void real_matrix_multiply_abt( RMatrix &a, RMatrix &b, RMatrix &c)
{
    double alpha=1.0;
    double beta=0.0;

    if (a.NCols()!=b.NCols())
    {
	cout<<"Problem in matrix multiply abt"<<endl;
	return;
    }
    c.resize(a.NRows(),b.NRows());
    
    dgemm('N', 'T' , int(a.NRows()) ,int(b.NRows()) ,int(a.NCols()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), int(b.NRows()),
          beta,&*c.begin(),int(c.NRows()));
}

//////////////////////////////////////////////////////////////////////////////

void real_matrix_times_vector(RMatrix &a,
                              std::vector<double> &b,
                              std::vector<double> &c)
{
    double alpha=1.0;
    double beta=0.0;
    int brows=b.size();
    
    //Matrix a1;
    //std::vector<double> b1;
    //a1=a;
    //b1=b;

    c.resize(brows);

    dgemm('N', 'N' , int(a.NRows()) ,1 ,int(a.NCols()) ,alpha ,
          &*a.begin(), int(a.NRows()), &*b.begin(), brows,
          beta,&*c.begin(),brows);   
}

//////////////////////////////////////////////////////////////////////////////

void real_vector_times_vector(std::vector<double> &a,           
                              std::vector<double> &b,
                              RMatrix &c)
{
    double alpha=1.0;
    double beta=0.0;
    int brows=b.size();

    c.resize(brows,brows);

    dgemm('N', 'T' , a.size(),b.size(),1,alpha ,
          &*a.begin(), a.size(), &*b.begin(), b.size(),
          beta,&*c.begin(),b.size());   
}

//////////////////////////////////////////////////////////////////////////////

void svd( RMatrix &a, 
          std::vector<double> &eigs, 
          RMatrix &u, RMatrix &vt)
{
    int info=0;
    int m=int(a.NRows());
    int n=int(a.NRows());
    std::vector<double> work(4*m);
    RMatrix a_copy;

    a_copy=a;
    work.assign(4*m,0.);

    u.resize(m,m);
    vt.resize(n,n);
    eigs.resize(min(m,n));

    cout<<"Performing SVD"<<endl;
    dgesvd('A','A',m,n,&*a_copy.begin(),m,
           &*eigs.begin(),&*u.begin(),m,
	   &*vt.begin(),n,
	   &*work.begin(),int(work.size()),info);
    
    cout<<"Finished SVD"<<endl;

}

double yt_A_x(RMatrix &y, RMatrix &A, RMatrix &x)
{
	RMatrix temp,temp2;
	real_matrix_multiply_atb( y, A, temp);
	real_matrix_multiply( temp, x, temp2);
	return temp2(0,0);
}
