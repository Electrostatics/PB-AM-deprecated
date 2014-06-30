#ifndef _MCOEFF_H_
#define _MCOEFF_H_

#define N_POLES 30
#define PREC_LIMIT (1e-30)

#include <complex>
#include <vector>
#include "util.h"

using namespace std;

class CTriCoeff;

class CMCoeff
{
public:
  enum TYPE {MPOLE, MPOLE_K, LOCAL, LOCAL_K};

  CMCoeff(const vector<REAL> & charges, const vector<CPnt>&  pos, 
	  int p, REAL rad = 0.0);
  CMCoeff(REAL ch, const CPnt & pos,int p, TYPE type);
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
 
  static REAL RS, IRS, KAPPA;
 
  

protected:
  static int IDX[2*N_POLES+1];
  TYPE m_type;
  vector<Complex> m_M, m_MU;
  int m_p, m_pU;
};

inline void
CMCoeff::saveUndo()
{
  m_MU.assign(m_M.begin(), m_M.begin()+IDX[m_p]); 
  m_pU = m_p; 
}

inline void
CMCoeff::undo()
{
  m_M.assign(m_MU.begin(), m_MU.begin()+IDX[m_pU]); 
  m_p = m_pU; 
}

inline void 
CMCoeff::reset(int p)
{
  m_M.assign(IDX[p], Complex());
  m_p = p;
}

inline void 
CMCoeff::copy(const CMCoeff & M, int p)
{
  assert(p <= M.m_p);
  m_M.assign(M.m_M.begin(), M.m_M.begin()+IDX[p]);
  m_p = p;
  m_type = M.m_type;
}

inline void 
CMCoeff::copy_p(const CMCoeff & M, int p)
{
  assert(p <= M.m_p && (p == m_p || p == m_p+1));
  m_M.resize(IDX[p-1]);
  m_M.insert(m_M.end(), M.m_M.begin()+IDX[p-1], M.m_M.begin()+IDX[p]);
  m_p = p;
}

inline void 
CMCoeff::operator*=(const REAL C[N_POLES])
{
  for (int n = 0; n < m_p; n++)
    for (int m = 0; m <= n; m++)
      m_M[IDX[n]+m] *= C[n];
}

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

inline CMCoeff 
operator*(REAL s, const CMCoeff & M)
{
  CMCoeff N(M.m_p);

  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = s * M.m_M[l];

  return N;
}

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


inline CMCoeff
conj(const CMCoeff & M)
{
  CMCoeff N(M.m_p, M.m_type);
  for (int l = 0; l < CMCoeff::IDX[M.m_p]; l++)
    N.m_M[l] = conj(M.m_M[l]);

  return N;
}

inline REAL 
inprod(const CMCoeff & M1,  const CMCoeff & M2)
{
  //assert(M1.m_type != M2.m_type);

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
