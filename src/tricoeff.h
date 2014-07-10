#ifndef _TRICOEFF_H_
#define _TRICOEFF_H_

#include "mcoeff.h"

using namespace std;

	//!  The CTriCoeff class 
/*!
		The tri coefficient class. Contains a trio of matrix coefficient objects
*/
class CTriCoeff
{
 public:
  CTriCoeff(int p = 0, int res = N_POLES); 
  CTriCoeff(const CTriCoeff & G) { *this = G; }
  
  const CTriCoeff & operator=(const CTriCoeff & G)
    { m_M[0] = G.m_M[0]; m_M[1] = G.m_M[1]; m_M[2] = G.m_M[2]; }
  void copy(const CTriCoeff & G, int p);
  void copy_p(const CTriCoeff & G, int p);
  
  void rotate(const CQuat & Q, int p);
  void incRotate(const CQuat & Q);
  void reset(int p = 0)
    { m_M[0].reset(p); m_M[1].reset(p); m_M[2].reset(p); }
  void recip()
    { m_M[0].recip(); m_M[1].recip(); m_M[2].recip(); }

  CTriCoeff & operator+=(const CTriCoeff & E)
    { m_M[0]+=E.m_M[0]; m_M[1]+=E.m_M[1]; m_M[2]+=E.m_M[2]; return *this; } 
  CTriCoeff & operator-=(const CTriCoeff & E)
    { m_M[0]-=E.m_M[0]; m_M[1]-=E.m_M[1]; m_M[2]-=E.m_M[2]; return *this; } 
  CTriCoeff & operator*=(const REAL * C)
    { m_M[0] *= C; m_M[1] *= C; m_M[2] *= C; return *this; } 
  
  friend CTriCoeff operator*(const CTriCoeff & G, const REAL C[N_POLES])
    { CTriCoeff N = G; N *= C; return N; }
  friend CTriCoeff operator+(const CTriCoeff & G1, const CTriCoeff & G2);
  friend CTriCoeff operator-(const CTriCoeff & G1, const CTriCoeff & G2);
  friend CPnt cross(const CTriCoeff & G1, const CTriCoeff & G2);
  friend CPnt inprod(const CMCoeff & M1, const CTriCoeff & M2);
  friend ostream & operator<<(ostream & out, const CTriCoeff & exp);
  friend CTriCoeff conj(const CTriCoeff & G);
  
  int getOrder() const
    { return m_M[0].getOrder(); }
  void setOrder(int p)
    { m_M[0].setOrder(p); m_M[1].setOrder(p); m_M[2].setOrder(p); }
  const CMCoeff & operator[](int c) const
    { assert(c >= 0 && c < 3); return m_M[c]; }
  CMCoeff & operator[](int c)
    { assert(c >= 0 && c < 3); return m_M[c]; }
  
  virtual void output() const;
  static REAL computeDev(const CTriCoeff & G1, const CTriCoeff & G2);

  void saveUndo();
  void undo();
  
protected:
  CMCoeff m_M[3];
};

///////////////////////////////////////////
/////// Inline functions

//!  The MPE setOrder function
/*!
		Function that sets the order ??  
*/
inline void
CTriCoeff::saveUndo()
{
  m_M[0].saveUndo();
  m_M[1].saveUndo();
  m_M[2].saveUndo();
}

inline void
CTriCoeff::undo()
{
  m_M[0].undo();
  m_M[1].undo();
  m_M[2].undo();
}

inline void 
CTriCoeff::copy(const CTriCoeff & G, int p)
{ 
  m_M[0].copy(G.m_M[0],p); 
  m_M[1].copy(G.m_M[1],p); 
  m_M[2].copy(G.m_M[2],p);
}

inline void 
CTriCoeff::copy_p(const CTriCoeff & G, int p)
{ 
  m_M[0].copy_p(G.m_M[0],p); 
  m_M[1].copy_p(G.m_M[1],p); 
  m_M[2].copy_p(G.m_M[2],p);
}

inline CTriCoeff
conj(const CTriCoeff & G)
{
  CTriCoeff H(G.getOrder());
  H.m_M[0] = conj(G.m_M[0]);
  H.m_M[1] = conj(G.m_M[1]);
  H.m_M[2] = conj(G.m_M[2]);

  return H;
}

inline CTriCoeff
operator+(const CTriCoeff & G1, const CTriCoeff & G2)
{
  CTriCoeff G(G1.getOrder());

  G.m_M[0] = G1.m_M[0] + G2.m_M[0];
  G.m_M[1] = G1.m_M[1] + G2.m_M[1];
  G.m_M[2] = G1.m_M[2] + G2.m_M[2];

  return G;
}

inline CTriCoeff
operator-(const CTriCoeff & G1, const CTriCoeff & G2)
{
  CTriCoeff G(G1.getOrder());

  G.m_M[0] = G1.m_M[0] - G2.m_M[0];
  G.m_M[1] = G1.m_M[1] - G2.m_M[1];
  G.m_M[2] = G1.m_M[2] - G2.m_M[2];

  return G;
}

inline CPnt
inprod(const CMCoeff & M1, const CTriCoeff & M2)
{
  assert(M1.getOrder() == M2.getOrder());

  return CPnt(inprod(M1, M2.m_M[0]), inprod(M1, M2.m_M[1]), 
	      inprod(M1, M2.m_M[2]));
}

inline CPnt 
cross(const CTriCoeff & G1, const CTriCoeff & G2)
{
  assert(G1.getOrder() == G2.getOrder());
  CPnt res(inprod(G1.m_M[1], G2.m_M[2]) - inprod(G1.m_M[2],G2.m_M[1]),
	   inprod(G1.m_M[2], G2.m_M[0]) - inprod(G1.m_M[0],G2.m_M[2]),
	   inprod(G1.m_M[0], G2.m_M[1]) - inprod(G1.m_M[1],G2.m_M[0]));

  return res;
}

inline REAL
CTriCoeff::computeDev(const CTriCoeff & G1, const CTriCoeff & G2)
{
  return (CMCoeff::computeDev(G1.m_M[0], G2.m_M[0]) +
	  CMCoeff::computeDev(G1.m_M[1], G2.m_M[1]) +
	  CMCoeff::computeDev(G1.m_M[2], G2.m_M[2]));
}

#endif
