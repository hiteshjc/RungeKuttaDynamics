#ifndef MATRIX_H
#define MATRIX_H

// Evgenii B. Rudnyi, http:// MatrixProgramming.com
#include<complex>
#include<vector>
#include<cstring>
using namespace std;

class Matrix : public std::vector< complex<double> >
{
  size_t m;
  size_t n;

public:
   Matrix() : m(0), n(0) {}
   Matrix(size_t m_, size_t n_) : std::vector< complex<double> >(m_*n_), m(m_), n(n_) {}
   Matrix(size_t m_, size_t n_, double val) : std::vector< complex<double> >(m_*n_, val), m(m_), n(n_) {}

  void resize(size_t m_, size_t n_)
    {m = m_; n = n_; std::vector< complex<double> >::resize(m_*n_);}
  void reserve(size_t m_, size_t n_)
    {std::vector< complex<double> >::reserve(m_*n_);}
  void clear()
    {m = n = 0; std::vector< complex<double> >::clear();}

  size_t NRows() const {return m;}
  size_t NCols() const {return n;}

  complex<double>& operator()(size_t i, size_t j)
  {
    return operator[](i + j*m);
  }
  const complex<double>& operator()(size_t i, size_t j) const
  {
    return operator[](i + j*m);
  }
	void swap( Matrix &y)
	{
		std::vector< complex<double> >::swap(y);
		std::swap(n, y.n);
		std::swap(m, y.m);
	}
	void clearMemory()
	{
		 Matrix empty;
	  swap(empty);
	}
};

class RMatrix : public std::vector<double>
{
  size_t m;
  size_t n;

public:
   RMatrix() : m(0), n(0) {}
   RMatrix(size_t m_, size_t n_) : std::vector<double>(m_*n_), m(m_), n(n_) {}
   RMatrix(size_t m_, size_t n_, double val) : std::vector<double>(m_*n_, val), m(m_), n(n_) {}

  void resize(size_t m_, size_t n_)
    {m = m_; n = n_; std::vector<double>::resize(m_*n_);}
  void reserve(size_t m_, size_t n_)
    {std::vector<double>::reserve(m_*n_);}
  void clear()
    {m = n = 0; std::vector<double>::clear();}

  size_t NRows() const {return m;}
  size_t NCols() const {return n;}

  double& operator()(size_t i, size_t j)
  {
    return operator[](i + j*m);
  }
  const double& operator()(size_t i, size_t j) const
  {
    return operator[](i + j*m);
  }
	void swap( RMatrix &y)
	{
		std::vector<double>::swap(y);
		std::swap(n, y.n);
		std::swap(m, y.m);
	}
	void clearMemory()
	{
		RMatrix empty;
	  swap(empty);
	}
};


#endif
