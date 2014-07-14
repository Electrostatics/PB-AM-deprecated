#include <cmath>

#include "mcoeff.h"
#include "shcoeff.h"

REAL CMCoeff::RS = 1.0;
REAL CMCoeff::IRS = 1.0;
REAL CMCoeff::KAPPA = 0.0;
int CMCoeff::IDX[2*N_POLES+1];

/******************************************************************/
/******************************************************************//**
* Construct matrix expansion coefficient object
******************************************************************/

CMCoeff::CMCoeff(const vector<REAL> & charges, const vector<CPnt> & pos, 
								 int p, REAL rad) : m_p(p)
{
  assert(p > 0);
  m_M.reserve(IDX[N_POLES]);
  m_M.resize(IDX[p]);
	
  if (rad == 0.0)
    m_type= MPOLE_K;
  else
    m_type = MPOLE;
	
  REAL I[p];
  CSHCoeff SH(p);
  for (int i = 0; i < pos.size(); i++)
	{
		if (charges[i] == 0.0)
			continue;
		
		CSpPnt spos = CartToSph(pos[i]);
		SH.reset(spos.theta(), spos.phi(), p);
		
		if (m_type == MPOLE_K)
			CSHCoeff::besseli(I, p, KAPPA*spos.rho());
		else
			CSHCoeff::besseli(I, p, 0.0);
		
		REAL k = charges[i];   
		REAL r = spos.rho()*IRS;
		for (int n = 0; n < p; n++)
		{
			I[n] *= k;
			k *= r;
		}
		
		(*this) += (SH * I);
	}
}

/******************************************************************/
/******************************************************************//**
*  Construct matrix expansion coefficient object
******************************************************************/

CMCoeff::CMCoeff(REAL ch, const CPnt & pos, int p, TYPE type) 
  : m_p(p), m_type(type)
{
  assert(p > 0);
  m_M.resize(IDX[p]);
  REAL k = ch; 
  CSpPnt spos = CartToSph(pos);
  CSHCoeff SH(p);
  SH.reset(spos.theta(), spos.phi(), p);
  REAL C[p];
	
  if (m_type == MPOLE || m_type == MPOLE_K)
	{
		if (m_type == MPOLE_K)
			CSHCoeff::besseli(C, p, KAPPA*spos.rho());
		else
			CSHCoeff::besseli(C, p, 0.0);  
		
		REAL r = spos.rho()*IRS;
		for (int n = 0; n < p; n++)
		{
			C[n] *= k;
			k *= r;
		}
	}
  else
	{
		REAL ir = 1.0/spos.rho();
		if (m_type == LOCAL_K)
		{
			CSHCoeff::besselk(C, p, KAPPA*spos.rho());
			k *=  exp(-KAPPA*spos.rho())*ir;
		}
		else
		{
			CSHCoeff::besselk(C, p, 0.0);
			k *= ir;
		}
		
		ir *= RS;  
		for (int n = 0; n < p; n++)
		{
			C[n] *= k;;
			k *= ir;
		}
	}
	
  (*this) += (SH * C);
}

/******************************************************************/
/******************************************************************//**
* Take partial derivative of matrix with respect to radius r
\return a MCoeff object of derivatives WRT r
******************************************************************/
CMCoeff
CMCoeff::dMdr(REAL r) const
{
  CMCoeff M(*this);
  REAL k[m_p];
  REAL ir = 1/r;
	
  if (m_type == MPOLE)
	{
		for (int n = 0; n < m_p; n++)
			k[n] = ir*n;
	}
  else if (m_type == LOCAL)
	{
		for (int n = 0; n < m_p; n++)
			k[n] = -ir*(n+1);
	}
  else if (m_type == MPOLE_K)
	{
		REAL I[m_p+1];
		REAL c = KAPPA*KAPPA*r*r;
		CSHCoeff::besseli(I, m_p+1, KAPPA*r);
		for (int n = 0; n < m_p; n++)
			k[n] = ir*(n + c*I[n+1]/(I[n]*(2*n+3)));
	}
  else
	{
		REAL K[m_p+1];
		CSHCoeff::besselk(K, m_p+1, KAPPA*r);
		
		for (int n = 0; n < m_p; n++)
			k[n] = ir*(n - (2*n+1)*K[n+1]/K[n]);
	}
	
  M *= k;
  return M;
}

/******************************************************************/
/******************************************************************//**
*  Take partial derivative of matrix with respect to theta angle
\return a MCoeff object of derivatives WRT theta
******************************************************************/

CMCoeff
CMCoeff::dMdt(REAL theta, REAL phi) const
{
  CMCoeff M(m_p);
  REAL cot_t = 1.0/tan(theta);
  Complex cis(cos(phi), -sin(phi));
	
  M(0,0) = 0.0;
  for (int n = 1; n < m_p; n++)
	{
		M(n,0) = -sqrt((REAL)n*(n+1))*(cis*m_M[IDX[n]+1]);
		for (int m = 1; m < n; m++)
			M(n,m) = m*cot_t*m_M[IDX[n]+m] -
			sqrt((REAL)(n-m)*(n+m+1))*(cis*m_M[IDX[n]+m+1]);
		
		M(n,n) = n*cot_t*m_M[IDX[n]+n];
	}
	
  return M; 
}

/******************************************************************/
/******************************************************************//**
*  Take partial derivative of matrix with respect to phi angle
\return a MCoeff object of derivatives WRT phi
******************************************************************/
CMCoeff
CMCoeff::dMdp() const
{
  CMCoeff M(m_p);
	
  for (int n = 0; n < m_p; n++)
    for (int m = 0; m <= n; m++)
      M(n,m) = Complex(-m_M[IDX[n]+m].imag()*m, m_M[IDX[n]+m].real()*m);
  
  return M;
}

/******************************************************************/
/******************************************************************//**
* Printing out the Matrix coefficient object
******************************************************************/

ostream & 
operator<<(ostream & out, const CMCoeff & M)
{
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  for (int n = 0; n < M.getOrder(); n++)
	{
		for (int m = 0; m <= n; m++)
		{
			REAL r = fabs(M(n,m).real())>1e-15 ? M(n,m).real() : 0;
			REAL im = fabs(M(n,m).imag())>1e-15 ? M(n,m).imag() : 0;
			cout << "(" << r << "," << im << ") | ";
		}
		cout << endl;
	}
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  return out;
}

/******************************************************************/
/******************************************************************//**
*  Printing out the Matrix coefficient object
******************************************************************/

void
CMCoeff::output(int p)
{
  int temp = m_p;
  m_p = p;
  cout << *this; 
  m_p = temp;
	
}

