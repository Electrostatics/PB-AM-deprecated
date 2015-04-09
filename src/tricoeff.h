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
 
//!  The CTriCoeff class constructor
/*!
		The tri coefficient class constructor.
		\param p is an integer number of poles
		\param res is an integer number of maximum poles */
  CTriCoeff(int p = 0, int res = N_POLES); 

//!  The CTriCoeff class constructor
/*!
		The tri coefficient class constructor from another TriCoeff object.
		\param G is an object of the TriCoeff class */
  CTriCoeff(const CTriCoeff & G) { *this = G; }
  
// Rotation functions	
  void rotate(const CQuat & Q, int p);
  void incRotate(const CQuat & Q);
	
// Printing functions
  friend ostream & operator<<(ostream & out, const CTriCoeff & exp);	
  virtual void output() const;

//!  The CTriCoeff reset function
/*! Function that resets the number of poles in each matrix coeff. */	
  void reset(int p = 0)
    { m_M[0].reset(p); m_M[1].reset(p); m_M[2].reset(p); }
		
  void recip()
    { m_M[0].recip(); m_M[1].recip(); m_M[2].recip(); }

//!  The CTriCoeff getOrder function
/*! Function that returns the order (npoles) of the TriCoefficients. */	
  int getOrder() const
    { return m_M[0].getOrder(); }
//!  The CTriCoeff setOrder function
/*! Function that sets the order (npoles) of the TriCoefficients. */			
  void setOrder(int p)
    { m_M[0].setOrder(p); m_M[1].setOrder(p); m_M[2].setOrder(p); }

// Operators
//!  The CTriCoeff = operator
/*! Function that sets on TriCoeff equal to another. */
  void operator=(const CTriCoeff & G)
    { m_M[0] = G.m_M[0]; m_M[1] = G.m_M[1]; m_M[2] = G.m_M[2]; }
//!  The CTriCoeff [] operator
/*! Function that returns the indexed matrix in the object. */
  const CMCoeff & operator[](int c) const
    { assert(c >= 0 && c < 3); return m_M[c]; }
//!  The CTriCoeff [] operator
/*! Function that returns the indexed matrix in the object. */
  CMCoeff & operator[](int c)
    { assert(c >= 0 && c < 3); return m_M[c]; }
//!  The CTriCoeff += operator
/*! Function that adds two TriCoeff objects together. */
  CTriCoeff & operator+=(const CTriCoeff & E)
    { m_M[0]+=E.m_M[0]; m_M[1]+=E.m_M[1]; m_M[2]+=E.m_M[2]; return *this; } 
//!  The CTriCoeff -= operator
/*! Function that subtracts two TriCoeff objects. */
  CTriCoeff & operator-=(const CTriCoeff & E)
    { m_M[0]-=E.m_M[0]; m_M[1]-=E.m_M[1]; m_M[2]-=E.m_M[2]; return *this; } 
//!  The CTriCoeff *= operator
/*! Function that multiplies a TriCoeff by a scalar C. */
  CTriCoeff & operator*=(const REAL * C)
    { m_M[0] *= C; m_M[1] *= C; m_M[2] *= C; return *this; }
//!  The CTriCoeff * operator
/*! Function that multiplies a TriCoeff by A vector of length N poles. */
  friend CTriCoeff operator*(const CTriCoeff & G, const REAL C[N_POLES])
    { CTriCoeff N = G; N *= C; return N; }

// Inline functions described below
  friend CTriCoeff conj(const CTriCoeff & G);
  friend CPnt cross(const CTriCoeff & G1, const CTriCoeff & G2);
  friend CPnt inprod(const CMCoeff & M1, const CTriCoeff & M2);
  friend CTriCoeff operator+(const CTriCoeff & G1, const CTriCoeff & G2);
  friend CTriCoeff operator-(const CTriCoeff & G1, const CTriCoeff & G2);	
  static REAL computeDev(const CTriCoeff & G1, const CTriCoeff & G2);

// Saving/Undoing the object
  void saveUndo();
  void undo();
	
// Copying the objects
  void copy(const CTriCoeff & G, int p);
  void copy_p(const CTriCoeff & G, int p);
  
protected:
  CMCoeff m_M[3];		//!< The three matrix coefficients that comprise the object
};

///////////////////////////////////////////
/////// Inline functions

//!  The CTriCoeff saveUndo function
/*! Function that saves the parameters of each of the three matrix objects
	within the TriCoeff object. */
inline void
CTriCoeff::saveUndo()
{
  m_M[0].saveUndo();
  m_M[1].saveUndo();
  m_M[2].saveUndo();
}

//!  The CTriCoeff undo function
/*! Function that reverts to saved states from each of the three 
	matrix objects within the TriCoeff object. */
inline void
CTriCoeff::undo()
{
  m_M[0].undo();
  m_M[1].undo();
  m_M[2].undo();
}

//!  The CTriCoeff copy function
/*! Function that copies an input triCoeff object into each of 
	the three matrix objects within the TriCoeff object. 
	\param G a triCoeff object to copy 
	\param p an integer number of poles in the TriCoeff object */
inline void 
CTriCoeff::copy(const CTriCoeff & G, int p)
{ 
  m_M[0].copy(G.m_M[0],p); 
  m_M[1].copy(G.m_M[1],p); 
  m_M[2].copy(G.m_M[2],p);
}

//!  The CTriCoeff undo function
/*! Function that copies an input triCoeff object into each of 
	the three matrix objects within the TriCoeff object. 
	\param G a triCoeff object to copy 
	\param p an integer number of poles in the TriCoeff object */
inline void 
CTriCoeff::copy_p(const CTriCoeff & G, int p)
{ 
  m_M[0].copy_p(G.m_M[0],p); 
  m_M[1].copy_p(G.m_M[1],p); 
  m_M[2].copy_p(G.m_M[2],p);
}

//!  The CTriCoeff conj function
/*! Function that returns the conjugate of an input G TriCoeff object
	\param G a TriCoeff object 
	\return H a Tricoeff object that is the conjugate of G  */
inline CTriCoeff
conj(const CTriCoeff & G)
{
  CTriCoeff H(G.getOrder());
  H.m_M[0] = conj(G.m_M[0]);
  H.m_M[1] = conj(G.m_M[1]);
  H.m_M[2] = conj(G.m_M[2]);

  return H;
}

//!  The CTriCoeff + operator
/*! Function that adds each matrix of the TriCoeff
		object to the second TriCoeff object.
		\param G1 first TriCoeff object to add
		\param G2 second TriCoeff to add
		\return G a TriCoeff object that is the sum of G1+G2 */
inline CTriCoeff
operator+(const CTriCoeff & G1, const CTriCoeff & G2)
{
  CTriCoeff G(G1.getOrder());

  G.m_M[0] = G1.m_M[0] + G2.m_M[0];
  G.m_M[1] = G1.m_M[1] + G2.m_M[1];
  G.m_M[2] = G1.m_M[2] + G2.m_M[2];

  return G;
}

//!  The CTriCoeff - operator
/*! Function that subtracts each matrix of the TriCoeff
		object to the second TriCoeff object.
		\param G1 first TriCoeff object to subtract from
		\param G2 second TriCoeff to subtract
		\return G a TriCoeff object that is the subtraction of G1-G2 */
inline CTriCoeff
operator-(const CTriCoeff & G1, const CTriCoeff & G2)
{
  CTriCoeff G(G1.getOrder());

  G.m_M[0] = G1.m_M[0] - G2.m_M[0];
  G.m_M[1] = G1.m_M[1] - G2.m_M[1];
  G.m_M[2] = G1.m_M[2] - G2.m_M[2];

  return G;
}

//!  The CTriCoeff inprod function
/*! Function that computes the inner product between a matrix coeff
			object and a TriCoeff object: [ M1 dot M2[0], M1 dot M2[1], M1 dot M2[2] ]
			\param M1 a MCoeff object
			\param M2 a TriCoeff object  
			\return A cartesian vector object of inner products, depicted above */
inline CPnt
inprod(const CMCoeff & M1, const CTriCoeff & M2)
{
  assert(M1.getOrder() == M2.getOrder());

  return CPnt(inprod(M1, M2.m_M[0]), inprod(M1, M2.m_M[1]), 
	      inprod(M1, M2.m_M[2]));
}

//!  The CTriCoeff cross function
/*! Function that computes the cross product between two TriCoeff matrices.
[ G1[1]*G2[2] - G1[2]*G2[1], G1[2]*G2[0] - G1[0]*G2[2], G1[0]*G2[1] - G1[1]*G2[0] ] 
		\param G1 the first TriCoeff object
		\param G2 the second TriCoeff object
		\return A CPnt object of the cross product between the two TriCoeffs*/
inline CPnt 
cross(const CTriCoeff & G1, const CTriCoeff & G2)
{
  assert(G1.getOrder() == G2.getOrder());
  CPnt res(inprod(G1.m_M[1], G2.m_M[2]) - inprod(G1.m_M[2],G2.m_M[1]),
	   inprod(G1.m_M[2], G2.m_M[0]) - inprod(G1.m_M[0],G2.m_M[2]),
	   inprod(G1.m_M[0], G2.m_M[1]) - inprod(G1.m_M[1],G2.m_M[0]));

  return res;
}

//!  The CTriCoeff computeDev function
/*! Function that computes the deviation between each mat object
		in two TriCoeff objects.
		\param G1 the first TriCoeff object
		\param G2 the second TriCoeff object 
		\return a floating point deviation, the sum of all three deviations */
inline REAL
CTriCoeff::computeDev(const CTriCoeff & G1, const CTriCoeff & G2)
{
  return (CMCoeff::computeDev(G1.m_M[0], G2.m_M[0]) +
	  CMCoeff::computeDev(G1.m_M[1], G2.m_M[1]) +
	  CMCoeff::computeDev(G1.m_M[2], G2.m_M[2]));
}

#endif
