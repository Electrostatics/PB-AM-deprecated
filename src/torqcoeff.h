#ifndef _TORQCOEFF_H_
#define _TORQCOEFF_H_

#define X 0
#define Y 1
#define Z 2

#include "tricoeff.h"

class CTorqCoeff : public CTriCoeff
{
 public:
  CTorqCoeff(int p = 0, int res = 0) : CTriCoeff(p, res) {} 
  CTorqCoeff(const CTriCoeff & G) : CTriCoeff(G) {} 
  CTorqCoeff(const vector<REAL> & charges, const vector<CPnt>&  pos, 
	     int p, REAL rad);

  /*  
const Complex dx(int n, int m) const
  { return m_M[0](n,m); }
  const Complex dy(int n, int m) const
  { return m_M[1](n,m); }
  const Complex dz(int n, int m) const
  { return m_M[2](n,m); }

  Complex & dx(int n, int m)
  { return m_M[0](n,m); }
  Complex & dy(int n, int m)
  { return m_M[1](n,m); }
  Complex & dz(int n, int m)
  { return m_M[2](n,m); }
  */
};




#endif
