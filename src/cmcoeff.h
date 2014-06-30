#ifndef _CMCOEFF_H_
#define _CMCOEFF_H_

#define N_POLES 20

#include <complex>
using namespace std;

class CMCoeff
{
public:
  CMCoeff(REAL charges[], const CSpPnt pos[], int num,
	  int p, REAL kappa) : m_p(p);
  CMCoeff() : m_p(N_POLES) { reset(); }
  CMCoeff(int p) m_p(p) {};
  CMCoeff(REAL theta, REAL phi, int p) : m_p(p);
  CMCoeff & operator=(const CMCoeff & M)
  {
    m_p = M.m_p;
    for (int n = 0; n < m_p; n++)
      for (int m = 0; m < n; m++)
	m_M[n][m] = M.m_M[n][m];
  }
  CMCoeff(const CMCoeff & M)
  { *this = M; }
    
  void reset()
  {
    for (int n = 0; n < p; n++)
      for (int m = 0; m < n; m++)
	m_M[n][m] = Complex(0.0, 0.0);
  }
  
  const Complex & operator()(int n, int m) const 
  {
    assert(abs(m) <= n && n < p);
    return (m >= 0 ? m_M[n][m] : recip(conj(m_M[n][m]), m); }
  }

  Complex & operator()(int n, int m)
  {
    assert(m >= 0 && m <= n && n < p);
    return m_M[n][m];
  }

  void operator*=(REAL C[N_POLES])
  {
    for (int n = 0; n < m_p; n++)
      for (int m = 0; m <= n; m++)
	m_M[n][m] *= C[n];
  }
  
  friend REAL inprod(const CMCoeff & M1,  const CMCoeff & M2);
  friend CMCoeff operator+(const CMCoeff & M1,  const CMCoeff & M2);

  static void legendre(REAL P[N_POLES][N_POLES], int n, REAL xval) const;
  static void besselk(REAL K[], int n, REAL val) const;
  static void besseli(REAL I[], int n, REAL val) const;
  
private:
  Complex recip(const Complex & c, int m)
  { return (m % 2 == 0 ? c : -c); }

  Complex m_M[N_POLES][N_POLES];
  int m_p;
};











#endif
