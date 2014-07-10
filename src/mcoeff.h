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
  CMCoeff(REAL ch, const CPnt & pos,int p, TYPE type);
	
//!  The mcoeff class constructor
/*!
		Initialize a MCoeff object.
		\param p is an int that is a minimal number of poles. Default is zero.
		\param res is an in that describes the total number of residuals (poles) for the expansion. Default is N_POLES
		\param type is an enum of TYPE that gives type of expansion coefficients.  Default is MPOLE
		\return an object of the MCoeff class
*/
  CMCoeff(int p = 0, int res = N_POLES, TYPE type = MPOLE) : 
	m_p(0), m_type(type)
  { reserve(res); reset(p); }
  
  CMCoeff & operator=(const CMCoeff & M)
  { m_M.assign(M.m_M.begin(), M.m_M.begin()+IDX[M.m_p]); m_p = M.m_p; }
  void copy(const CMCoeff & M, int p);
  void copy_p(const CMCoeff & M, int p);
	
  CMCoeff(const CMCoeff & M) : m_type(M.m_type)
  { *this = M; }
  
  CMCoeff dMdr(REAL r) const;
  CMCoeff dMdt(REAL theta, REAL phi) const;
  CMCoeff dMdp() const;
  void recip()
  { for (int l = 0; l < IDX[m_p]; l++) m_M[l] = -m_M[l]; }
  void reset(int p = 0);
  void reserve(int p)
	{if (p > 0) m_M.reserve(IDX[p]); }
  
  const Complex operator()(int n, int m) const 
	{
		assert(abs(m) <= n && n < m_p);
		return (m >= 0 ? m_M[IDX[n]+m] : conj(m_M[IDX[n]-m])); 
	}
  Complex & operator()(int n, int m)
	{ assert(m >= 0 && m <= n && n < m_p); return m_M[IDX[n]+m]; }
  
  void operator*=(const REAL C[N_POLES]);
  void operator+=(const CMCoeff & E);
  void operator-=(const CMCoeff & E);
  
  friend ostream & operator<<(ostream & out, const CMCoeff & p);
  friend CPnt inprod(const CMCoeff & M1, const CTriCoeff & M2);
  friend REAL inprod(const CMCoeff & M1, const CMCoeff & M2);
  friend CMCoeff operator+(const CMCoeff & M1, const CMCoeff & M2);
  friend CMCoeff operator-(const CMCoeff & M1, const CMCoeff & M2);
  friend CMCoeff operator*(const CMCoeff & M, const REAL C[N_POLES])
  { CMCoeff N = M; N *= C; return N; }
  friend CMCoeff operator*(REAL s, const CMCoeff & M);
  friend CMCoeff conj(const CMCoeff & M);
  
  int getOrder() const
  { return m_p; }
  void setOrder(int p)
  { 
    m_M.resize(IDX[p]);
    m_p = p; 
  }
	
  void saveUndo();
  void undo();
  void output(int p);
	
  REAL change() const;
  static REAL computeDev(const CMCoeff & M1, const CMCoeff & M2);
	
  static REAL RS;							//!< Scaling factor for system (~ protein radius)
	static REAL IRS;						//!< Inverse scaling factor (=1/RS)
	static REAL KAPPA;					//!< Inverse Debye length
	
protected:
  static int IDX[2*N_POLES+1];
  TYPE m_type;								//<! Enum indicating the type of MCoeff given.
  vector<Complex> m_M;
	vector<Complex> m_MU;
  int m_p;										//!< Number of poles for the expansion.
	int m_pU;
}; // end CMCoeff

//////////////////////////////////////////
// Inline functions
//////////////////////////////////////////

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline void
CMCoeff::saveUndo()
{
  m_MU.assign(m_M.begin(), m_M.begin()+IDX[m_p]); 
  m_pU = m_p; 
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline void
CMCoeff::undo()
{
  m_M.assign(m_MU.begin(), m_MU.begin()+IDX[m_pU]); 
  m_p = m_pU; 
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline void 
CMCoeff::reset(int p)
{
  m_M.assign(IDX[p], Complex());
  m_p = p;
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline void 
CMCoeff::copy(const CMCoeff & M, int p)
{
  assert(p <= M.m_p);
  m_M.assign(M.m_M.begin(), M.m_M.begin()+IDX[p]);
  m_p = p;
  m_type = M.m_type;
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline void 
CMCoeff::copy_p(const CMCoeff & M, int p)
{
  assert(p <= M.m_p && (p == m_p || p == m_p+1));
  m_M.resize(IDX[p-1]);
  m_M.insert(m_M.end(), M.m_M.begin()+IDX[p-1], M.m_M.begin()+IDX[p]);
  m_p = p;
}

//!  The MCoeff expansion class operator *
/*!
		Function that multiplies the MCoeff by the Gamma operator  
*/
inline void 
CMCoeff::operator*=(const REAL C[N_POLES])
{
  for (int n = 0; n < m_p; n++)
    for (int m = 0; m <= n; m++)
      m_M[IDX[n]+m] *= C[n];
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
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

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
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

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline CMCoeff 
operator*(REAL s, const CMCoeff & M)
{
  CMCoeff N(M.m_p);
	
  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = s * M.m_M[l];
	
  return N;
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
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

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
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

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
inline CMCoeff
conj(const CMCoeff & M)
{
  CMCoeff N(M.m_p, M.m_type);
  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = conj(M.m_M[l]);
	
  return N;
}

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
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
		A function that computes the deviation between 2 Matrix coeff objects
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

//!  The mcoeff expansion class
/*!
		The class that contains all details for a mcoeff ??  
*/
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
