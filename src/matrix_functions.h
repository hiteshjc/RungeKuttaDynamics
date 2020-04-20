#ifndef MATRIX_FUNCTIONS_H
#define MATRIX_FUNCTIONS_H

#include"global.h"
#include"matrix.h"

using namespace std;
double yt_A_x(RMatrix &y, RMatrix &A, RMatrix &x);

// Diagonalize A and report eigenvalues and orthogonal eigenvecs
void real_symmetric_diagonalize(RMatrix &a, 
                 std::vector<double> &eigs, 
                 RMatrix &eigenvecs);

void general_real_diagonalize
                (Matrix &a, 
                 std::vector< complex<double> > &eigs, 
                 Matrix &eigenvecs);

// Diagonalize symmetric matrix A and report eigenvalues and orthogonal eigenvecs

void symmetric_diagonalize
                (Matrix &a, 
                 std::vector<double> &eigs, 
                 Matrix &eigenvecs);

// Orthogonalize eigenvectors

void gram_schmidt(Matrix &eigenvecs, std::vector< complex<double> > &eigs);

void eta_orthogonalize( std::vector< std::vector<double> > &eigenvecs,
			std::vector<double> const &eta);

// C=A*B

void real_matrix_multiply( RMatrix &a, 
                           RMatrix &b, 
                           RMatrix &c);

void real_matrix_multiply_atb( RMatrix &a, RMatrix &b, RMatrix &c);
void real_matrix_multiply_abt( RMatrix &a, RMatrix &b, RMatrix &c);

void perform_lu(Matrix &mat, Matrix &lu, std::vector<int> &ipiv);
void invert_real_matrix(Matrix &mat, Matrix &inverse);

// Multiply matrix times a vector

void real_matrix_times_vector(RMatrix &a,
                              std::vector<double> &b,
                              std::vector<double> &c);

// Multiply vector times vector
void real_vector_times_vector(std::vector<double> &a,
                              std::vector<double> &b,
                              RMatrix &c);

// Singular Value decomposition

void svd( RMatrix &a, 
          std::vector<double> &eigs, 
          RMatrix &u, RMatrix &vt);

void matrix_lanczos(int 			      		iterations,
		    int 			      		how_many_eigenvecs, 
                    Matrix 					&a,
		    std::vector<double> 			&eigs,
		    Matrix 					&eigenvecs,
	            bool 					ipr);
#endif
