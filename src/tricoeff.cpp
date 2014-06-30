#include <cmath>
#include "tricoeff.h"

CTriCoeff::CTriCoeff(int p, int res)
{
  assert(res >= 0 && p >= 0);
  if (res > 0)
    {
      m_M[0].reserve(res); 
      m_M[1].reserve(res); 
      m_M[2].reserve(res); 
    }

  m_M[0].reset(p);
  m_M[1].reset(p);
  m_M[2].reset(p);
}

void
CTriCoeff::rotate(const CQuat & Q, int p)
{
  assert(p <= getOrder());
  CTriCoeff G(p);
  
  for (int n = 0; n < p; n++)
    for (int m = 0; m <= n; m++)
      {
	CPnt ar(m_M[0](n,m).real(), m_M[1](n,m).real(), m_M[2](n,m).real());
	CPnt ai(m_M[0](n,m).imag(), m_M[1](n,m).imag(), m_M[2](n,m).imag());
	CPnt br = Q*ar;
	CPnt bi = Q*ai;

	G.m_M[0](n,m) = Complex(br.x(), bi.x());
	G.m_M[1](n,m) = Complex(br.y(), bi.y());
	G.m_M[2](n,m) = Complex(br.z(), bi.z());
      }

  this->copy(G, p);
}

void
CTriCoeff::incRotate(const CQuat & Q)
{
  assert(getOrder() > 0);
  CTriCoeff G(getOrder());

  int n = getOrder()-1;
  for (int m = 0; m <= n; m++)
    {
      CPnt ar(m_M[0](n,m).real(), m_M[1](n,m).real(), m_M[2](n,m).real());
      CPnt ai(m_M[0](n,m).imag(), m_M[1](n,m).imag(), m_M[2](n,m).imag());
      CPnt br = Q*ar;
      CPnt bi = Q*ai;
      
      G.m_M[0](n,m) = Complex(br.x(), bi.x());
      G.m_M[1](n,m) = Complex(br.y(), bi.y());
      G.m_M[2](n,m) = Complex(br.z(), bi.z());
    }
  
  this->copy_p(G, getOrder());
}

ostream & 
operator<<(ostream & out, const CTriCoeff & T)
{
  cout << "\t---000---" << endl;
  cout << T.m_M[0] << endl;

  cout << "\t---111---" << endl;
  cout << T.m_M[1] << endl;

  cout << "\t---222---" << endl;
  cout << T.m_M[2] << endl;
  
  return out;
}

void
CTriCoeff::output() const
{
  cout << *this; 
}
