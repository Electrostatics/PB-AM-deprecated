#ifndef _MCOEFF_H_
#define _MCOEFF_H_

#define N_POLES 30
#define PREC_LIMIT (1e-30)

#include <complex>
#include <vector>
#include "util.h"

using namespace std;

class CTriCoeff;

//!  The mcoeff expansion class
/*!
		The class that contains all details for a Matrix expansion
		and its coefficients 
*/

class CMCoeff
{
public:

//! Enum of MCoeff object type.
  enum TYPE {
							MPOLE, 
							MPOLE_K, 
							LOCAL, 
							LOCAL_K
	};

// Functions in mcoeff.cpp
	
	//!  The MCoeff class constructor
	/*!
		Initialize a MCoeff object.
		\param charges a vector of floating point charges
		\param pos a vector of cartesian coordinates, one per charge
		\param p an integer that generally represents the number of poles
		\param rad a floating point number that represents the radius of the CG molecule. Default is zero.
		\return an object of the MCoeff class */
  CMCoeff(const vector<REAL> & charges, const vector<CPnt>&  pos, 
					int p, REAL rad = 0.0);
	//!  The MCoeff class constructor
	/*!
		Initialize a MCoeff object.
		\param charges a floating point number representin a partial charge
		\param pos a vector of cartesian coordinates for the XYZ coord of the charge
		\param p an integer that generally represents the number of poles
		\param type an Enum indicating the type of expansion class to be made
		\return an object of the MCoeff class */
  CMCoeff(REAL ch, const CPnt & pos,int p, TYPE type);

	// Functions for computing derivatives! 
  CMCoeff dMdr(REAL r) const;
  CMCoeff dMdt(REAL theta, REAL phi) const;
  CMCoeff dMdp() const;
	
	// Printing functions
  void output(int p);	  
  friend ostream & operator<<(ostream & out, const CMCoeff & p);

	
// Inline functions
	//!  The mcoeff class constructor
	/*!
		Initialize a MCoeff object.
		\param p is an int that is a minimal number of poles. Default is zero.
		\param res is an in that describes the total number of residuals (poles) for the expansion. Default is N_POLES
		\param type is an enum of TYPE that gives type of expansion coefficients.  Default is MPOLE
		\return an object of the MCoeff class */
  CMCoeff(int p = 0, int res = N_POLES, TYPE type = MPOLE) : 
	m_p(0), m_type(type)
  { reserve(res); reset(p); }
	
	//!  The mcoeff class constructor
	/*!
		Initialize a MCoeff object.
		\param M is an input MCoeff object
		\return an object of the MCoeff class */
  CMCoeff(const CMCoeff & M) : m_type(M.m_type)
  { *this = M; }
 
	//!  The mcoeff class = operator
	/*! Set the object to the left of = sign to the MCoeff object on the right of it */	 
  void operator=(const CMCoeff & M)
  { m_M.assign(M.m_M.begin(), M.m_M.begin()+IDX[M.m_p]); m_p = M.m_p; }
	
	//!  The mcoeff class recip fuction
	/*! Multiplies each value of the Matrix object by -1 */
  void recip()
  { for (int l = 0; l < IDX[m_p]; l++) m_M[l] = -m_M[l]; }

	//!  The mcoeff class reserve fuction
	/*! Reallocates the size of the coefficient matrix
			\param p an int of the position in the index to indicate
							what size allocation is desired for coeff matrix */	
	void reserve(int p)
	{if (p > 0) m_M.reserve(IDX[p]); }
  
	//!  The mcoeff class () operator
	/*! Returns the n,mth index of the coefficient matrix
			\param n an index of the matrix, runs from 0 to npoles
			\param m an index of the matrix, runs from -n to +n */	
  const Complex operator()(int n, int m) const 
	{
		assert(abs(m) <= n && n < m_p);
		return (m >= 0 ? m_M[IDX[n]+m] : conj(m_M[IDX[n]-m])); 
	}
	//!  The mcoeff class () operator
	/*! Returns the n,mth index of the coefficient matrix
			\param n an index of the matrix, runs from 0 to npoles
			\param m an index of the matrix, runs from -n to +n */	
  Complex & operator()(int n, int m)
	{ assert(m >= 0 && m <= n && n < m_p); return m_M[IDX[n]+m]; }
  
  friend CMCoeff operator*(const CMCoeff & M, const REAL C[N_POLES])
  { CMCoeff N = M; N *= C; return N; }
 
	//!  The mcoeff class getOrder fuction
	/*! Returns the number of poles in the MCoeff object */	 
  int getOrder() const
  { return m_p; }
	//!  The mcoeff class setOrder fuction
	/*! Sets the order of poles and resizes the matrix
		\param p an integer of number of desired poles  */
  void setOrder(int p)
  { 
    m_M.resize(IDX[p]);
    m_p = p; 
  }


// Inline functions, written after class below
	// Save and undo
  void saveUndo();
  void undo();

// Reset Matrix coefficients
	void reset(int p = 0);
	
	// Copy functions
  void copy(const CMCoeff & M, int p);
  void copy_p(const CMCoeff & M, int p);

	// Operators
  void operator*=(const REAL C[N_POLES]);
  void operator+=(const CMCoeff & E);
  void operator-=(const CMCoeff & E);
  friend CMCoeff operator*(REAL s, const CMCoeff & M);
  friend CMCoeff operator+(const CMCoeff & M1, const CMCoeff & M2);
  friend CMCoeff operator-(const CMCoeff & M1, const CMCoeff & M2);
  friend CPnt inprod(const CMCoeff & M1, const CTriCoeff & M2);
  friend REAL inprod(const CMCoeff & M1, const CMCoeff & M2);
	friend CMCoeff conj(const CMCoeff & M);
	
	// For computing change in iteration rounds
  static REAL computeDev(const CMCoeff & M1, const CMCoeff & M2);
  REAL change() const;
	
// Variables
	
  static REAL RS;							//!< Scaling factor for system (~ protein radius)
	static REAL IRS;						//!< Inverse scaling factor (=1/RS)
	static REAL KAPPA;					//!< Inverse Debye length
	
protected:
  static int IDX[2*N_POLES+1];//<! An index vector for 2*Pol + 1 values
  TYPE m_type;								//<! Enum indicating the type of MCoeff given.
  vector<Complex> m_M;				//!< A vector of complex numbers of matrix coefficients
	vector<Complex> m_MU;				//!< A vector of complex numbers of matrix coefficients, saved incase of undo
  int m_p;										//!< Number of poles for the expansion.
	int m_pU;										//!< Number of poles for the expansion, saved incase of undo
}; // end CMCoeff

//////////////////////////////////////////
// Inline functions
//////////////////////////////////////////

//!  The mcoeff saveUndo function
/*! A function that saves the M matrix in the matrix MU
for storage in the case of an undo.  It also saves the current number 
of poles */
inline void
CMCoeff::saveUndo()
{
  m_MU.assign(m_M.begin(), m_M.begin()+IDX[m_p]); 
  m_pU = m_p; 
}

//!  The mcoeff undo function
/*! A function that reverts the M matrix to the matrix MU.  
		It also reverts the number of poles */
inline void
CMCoeff::undo()
{
  m_M.assign(m_MU.begin(), m_MU.begin()+IDX[m_pU]); 
  m_p = m_pU; 
}

//!  The mcoeff reset function
/*! A function that resets the matrix coefficients at input integer
		to all zero, and sets the number of poles to given input.
		\param p an integer of desired number of poles. */
inline void 
CMCoeff::reset(int p)
{
  m_M.assign(IDX[p], Complex());
  m_p = p;
}

//!  The mcoeff copy function
/*! A function that copies the number of poles and matrix
	coefficients from an input MCoeff object to current MCoeff
	object.
		\param M an input of a coefficient matrix to copy
		\param p an integer of desired number of poles. */
inline void 
CMCoeff::copy(const CMCoeff & M, int p)
{
  assert(p <= M.m_p);
  m_M.assign(M.m_M.begin(), M.m_M.begin()+IDX[p]);
  m_p = p;
  m_type = M.m_type;
}

//!  The mcoeff copy_p function
/*! A function that copies one level of poles to the end 
		of the MCoeff being operated on. 
		\param M a MCoeff object to copy from
		\param p an int of the pole number to copy */
inline void 
CMCoeff::copy_p(const CMCoeff & M, int p)
{
  assert(p <= M.m_p && (p == m_p || p == m_p+1));
  m_M.resize(IDX[p-1]);
  m_M.insert(m_M.end(), M.m_M.begin()+IDX[p-1], M.m_M.begin()+IDX[p]);
  m_p = p;
}

//!  The MCoeff expansion class operator *=
/*! Function that multiplies the MCoeff by the Gamma operator */
inline void 
CMCoeff::operator*=(const REAL C[N_POLES])
{
  for (int n = 0; n < m_p; n++)
    for (int m = 0; m <= n; m++)
      m_M[IDX[n]+m] *= C[n];
}

//!  The MCoeff expansion class operator +=
/*! Function that adds one MCoeff object to another */
inline  void 
CMCoeff::operator+=(const CMCoeff & E)
{
  if (m_p < E.m_p)
	{
		m_M.resize(IDX[E.m_p]);
		m_p = E.m_p;
	}
	
  for (int l = 0; l < IDX[E.m_p]; l++)
    m_M[l] += E.m_M[l];
}

//!  The MCoeff expansion class operator -=
/*! Function that subtracts one MCoeff object from another */
inline  void 
CMCoeff::operator-=(const CMCoeff & E)
{
	if (m_p < E.m_p)
	{
		m_M.resize(IDX[E.m_p]);
		m_p = E.m_p;
	}
	
  for (int l = 0; l < IDX[E.m_p]; l++)
    m_M[l] -= E.m_M[l];
}

//!  The MCoeff expansion class operator *
/*! Function that multiplies the MCoeff by a floating point */
inline CMCoeff 
operator*(REAL s, const CMCoeff & M)
{
  CMCoeff N(M.m_p);
	
  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = s * M.m_M[l];
	
  return N;
}

//!  The MCoeff expansion class operator +
/*! Function that adds one MCoeff object to another */
inline CMCoeff
operator+(const CMCoeff & M1, const CMCoeff & M2)
{
  assert(M1.m_type == M2.m_type);
	
  if (M1.m_p >= M2.m_p)
	{
		CMCoeff S(M1);
		S += M2;
		return S;
	}
  else
	{
		CMCoeff S(M2);
		S += M1;
		return S;
	}
}

//!  The MCoeff expansion class operator -
/*! Function that subtracts one MCoeff object from another */
inline CMCoeff
operator-(const CMCoeff & M1, const CMCoeff & M2)
{
  assert(M1.m_type == M2.m_type);
	
  if (M1.m_p >= M2.m_p)
	{
		CMCoeff S(M1);
		S -= M2;
		return S;
	}
  else
	{
		CMCoeff S(M2.m_p);
		S -= M2;
		S += M1;
		return S;
	}
}

//!  The mcoeff expansion conj function
/*! Function that computes and returns the complex
conjugate of a matrix. */
inline CMCoeff
conj(const CMCoeff & M)
{
  CMCoeff N(M.m_p, M.m_type);
  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = conj(M.m_M[l]);
	
  return N;
}

//!  The mcoeff inprod function
/*! Function that computes the inner product of two MCoeff objects
	\param M1 one MCoeff input
	\param M2 the other MCoeff input
*/
inline REAL 
inprod(const CMCoeff & M1,  const CMCoeff & M2)
{

  int p = (M1.m_p > M2.m_p ? M2.m_p : M1.m_p);
  REAL sum1(0.0), sum2(0.0);
  for (int n = 0; n < p; n++)
	{
		for (int m = 1; m <= n; m++)
			sum1 += M1(n,m).real()*M2(n,m).real() +
			M1(n,m).imag()*M2(n,m).imag();
		
		sum2 += M1(n,0).real()*M2(n,0).real();
	}
	
  return 2*sum1 + sum2;
}

//!  The mcoeff expansion computeDev function
/*!
		A function that computes the deviation between 2 Matrix coeff objects.
		Represents EQ 52 in paper, somewhat.
		\param M1 a pointer to an MCoeff object
		\param M2 a pointer to a second MCoeff object for comparison to first
		\return a floating point number of the deviation between the two objects'
								 matrix coefficients
*/
inline REAL
CMCoeff::computeDev(const CMCoeff & M1, const CMCoeff & M2)
{
  REAL sum = 0;
  
  assert(M1.m_p == M2.m_p);
	
  for (int n = 0; n < M2.m_p; n++)
    for (int m = 0; m <= n; m++)
		{
			if (M1(n,m) == Complex(0.0,0.0) && M2(n,m) == Complex(0.0,0.0))
				continue;
			
			Complex a;
			if (fabs(M1(n,m).real()) < PREC_LIMIT && 
					fabs(M1(n,m).imag()) <PREC_LIMIT)
				a = M2(n,m);
			else if (fabs(M2(n,m).real()) < PREC_LIMIT && 
							 fabs(M2(n,m).imag()) <PREC_LIMIT)
				a = M1(n,m);
			else
				a = 0.5*(M1(n,m) - M2(n,m))/(M1(n,m) + M2(n,m));
			
			sum += a.real()*a.real() + a.imag()*a.imag();
		}
  
  return sum;
}

//!  The mcoeff change function
/*! Computes the difference between the saved Matrix and the current one,
			This and the function computeDev account for EQ 52 in 
			the Lotan 2006 paper */
inline REAL
CMCoeff::change() const
{
  CMCoeff tmp;
  tmp.setOrder(m_p);
  tmp.m_M = m_MU;
	
  REAL dev = CMCoeff::computeDev(*this, tmp);
  return sqrt(dev/(m_p*(m_p+1)*0.5));
}


#endif
