#ifndef _GRADCOEFF_H_
#define _GRADCOEFF_H_

#define dRHO 0
#define dTHETA 1
#define dPHI 2

#define dX 0
#define dY 1
#define dZ 2

#include "tricoeff.h"

using namespace std;

//!  The GradCoeff expansion class
/*!
		The class that contains all details for gradient coefficients.  
*/
class CGradCoeff : public CTriCoeff
{
 public: 
  CGradCoeff(int p = 0, int res = N_POLES) : CTriCoeff(p, res) {} 
  CGradCoeff(const CTriCoeff & G) : CTriCoeff(G) {} 
  CGradCoeff(const CMCoeff & M, const CSpPnt & c); 
  
  CGradCoeff sphToCart(const CPnt * R);
  
  friend ostream & operator<<(ostream & out, const CGradCoeff & G);

  const Complex dr(int n, int m) const
  { return m_M[dRHO](n,m); }
  const Complex dt(int n, int m) const
  { return m_M[dTHETA](n,m); }
  const Complex dp(int n, int m) const
  { return m_M[dPHI](n,m); }

  Complex & dr(int n, int m)
  { return m_M[dRHO](n,m); }
  Complex & dt(int n, int m)
  { return m_M[dTHETA](n,m); }
  Complex & dp(int n, int m)
  { return m_M[dPHI](n,m); }
};

inline CGradCoeff 
CGradCoeff::sphToCart(const CPnt * R)
{
  CGradCoeff G(getOrder());

  for (int i = 0; i < 3; i++)
    {
      G.m_M[0] += R[i].x()*m_M[i];
      G.m_M[1] += R[i].y()*m_M[i];
      G.m_M[2] += R[i].z()*m_M[i];
    }

  return G;
}

#endif
