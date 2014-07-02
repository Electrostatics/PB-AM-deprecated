#include "shcoeff.h"
#include "transcoeff.h"

REAL CTransCoeff::m_alpha[N_POLES*2][N_POLES];
REAL CTransCoeff::m_beta[N_POLES*2][N_POLES];
REAL CTransCoeff::m_gamma[N_POLES*2][2*N_POLES-1];
REAL CTransCoeff::m_delta[N_POLES*2][2*N_POLES-1];
bool CTransCoeff::m_even[4*N_POLES];

/******************************************************************/
/******************************************************************//**
* Initialize the rotation coefficients (includes coeff equs from 
* 2006 paper
******************************************************************/


void
CTransCoeff::initConstants()
{
  REAL kapsqr = CMCoeff::KAPPA*CMCoeff::KAPPA*CMCoeff::RS*CMCoeff::RS;  //!< K^2 * RS^2

  for (int n = 0; n < 2*N_POLES; n++)
    for (int m = 0; m < N_POLES; m++)
		{
			ALPHA(n,m) = sqrt((REAL)(n+m+1)*(n-m+1));						// Given as alpha(n,m) in Lotan, 2006 (eq 1.9)
			BETA(n,m)  = kapsqr*ALPHA(n,m)/((2*n+1)*(2*n+3));		// Given as beta(n,m) in Lotan, 2006 (eq 1.9)

			GAMMA(n,m) = sqrt((REAL)(n-m-1)*(n-m));							// Given as eta(n,m) in Lotan, 2006 (eq 1.9)
			DELTA(n,m) = kapsqr*GAMMA(n,m)/((2*n+1)*(2*n-1));		// Given as mu(n,m) in Lotan, 2006 (eq 1.9)
			if (m != 0)
			{
				GAMMA(n,-m) = -sqrt((REAL)(n+m-1)*(n+m));					// Given as eta(n,m) in Lotan, 2006 (eq 1.9)
				DELTA(n,-m) = kapsqr*GAMMA(n,-m)/((2*n+1)*(2*n-1));// Given as mu(n,m) in Lotan, 2006 (eq 1.9)
			}
		}

	
  EVEN(0) = true;
  for (int n = 1; n < 4*N_POLES; n++)
    EVEN(n) = !EVEN(n-1);								// determining whether each n is even or not
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/


CTransCoeff::CTransCoeff(bool bGrad) : m_bGrad(bGrad), m_p(1)
{
  m_T.resize(1);
  m_T[0] = new REAL*[1];
  m_T[0][0] = new REAL[1];

  if (m_bGrad)
    {
      m_dT.resize(1);
      m_dT[0] = new REAL*[1];
      m_dT[0][0] = new REAL[1];
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::allocate()
{
  assert(m_T.size() < 2*N_POLES-1);
  int p = m_T.size();
  
  m_T.resize(p+2);
  m_T[p] = new REAL*[p+1];
  m_T[p+1] = new REAL*[p+2];
  for (int j = 0; j < p+1; j++)
    m_T[p][j] = new REAL[p+1];
  for (int j = 0; j < p+2; j++)
    m_T[p+1][j] = new REAL[p+2];
  
  if (m_bGrad)
    {
      m_dT.resize(p+2);
      m_dT[p] = new REAL*[p+1];
      m_dT[p+1] = new REAL*[p+2];
      for (int j = 0; j < p+1; j++)
	m_dT[p][j] = new REAL[p+1];
      for (int j = 0; j < p+2; j++)
	m_dT[p+1][j] = new REAL[p+2];
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void 
CTransCoeff::deallocate()
{
  assert(m_T.size() > 1);
  int p = m_T.size();

  for (int j = 0; j < p; j++)
    delete[] m_T[p-1][j];
  for (int j = 0; j < p-1; j++)
    delete[] m_T[p-2][j];
  delete[] m_T[p-1];
  delete[] m_T[p-2];
  m_T.resize(p-2);

  if (m_bGrad)
    {
      for (int j = 0; j < p; j++)
	delete[] m_dT[p-1][j];
      for (int j = 0; j < p-1; j++)
	delete[] m_dT[p-2][j];
      delete[] m_dT[p-1];
      delete[] m_dT[p-2];
      m_dT.resize(p-2);
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::reallocate(int p)
{
  assert(p <= N_POLES && p >= 1);
  if (p - m_p > 0)
    for (int i = m_p; i < p; i++)
      allocate();
  else
    for (int i = m_p; i > p; i--)
      deallocate();
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::initParams(REAL d)
{
  m_d = d; 
  m_kd = CMCoeff::KAPPA * m_d;
  
  if (d > 0.0)
    {
      m_ir = 1.0/m_d;
      m_exkid = exp(-m_kd)*m_ir;  
    }

  CSHCoeff::besselk(m_K, 2*m_p, m_kd);
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void 
CTransCoeff::reset(REAL rho, int p)
{
  assert(p > 0 && p <= N_POLES);

  reallocate(p);

  m_p = p;
  initParams(rho);  
  computeCoeff();
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

// Apply the translation operator to the MP coeffs
void
CTransCoeff::translate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose)
{
  assert(p <= N_POLES && Min.getOrder() >= p);
  Mout.reset(p);
    
  for (int n = 0; n < p; n++)
    {
      for (int m = 0; m <= n; m++)
	{
	  int l;
	  if (tpose)
	    {
	      for (l = m; l <= n; l++)
		Mout(n,m) += TRANS(n,l,m)*Min(l,m);
	       	      
	      for (; l < p; l++)
		{
		  REAL T = (EVEN(n+l) ? TRANS(l,n,m) : -TRANS(l,n,m));
		  Mout(n,m) += T*Min(l,m);
		}
	    }
	  else
	    {
	      for (l = m; l <= n; l++)
		{
		  REAL T = (EVEN(n+l) ? TRANS(n,l,m) : -TRANS(n,l,m));
		  Mout(n,m) += T*Min(l,m);
		}
	      
	      for (; l < p; l++)
		Mout(n,m) += TRANS(l,n,m)*Min(l,m);
	    }
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

// Apply the derivative of the translation operator 
void
CTransCoeff::dTranslate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose)
{
  assert(p <= N_POLES && Min.getOrder() >= p);
  Mout.reset(p);
  
  for (int n = 0; n < p; n++)
    {
      for (int m = 0; m <= n; m++)
	{
	  int l;
	  if (tpose)
	    {
	      for (l = m; l <= n; l++)
		Mout(n,m) += dTRANS(n,l,m)*Min(l,m);
	      
	      for (; l < p; l++)
		{
		  REAL dT = (EVEN(n+l) ? dTRANS(l,n,m) : -dTRANS(l,n,m));
		  
		  Mout(n,m) += dT*Min(l,m);
		}
	    }
	  else
	    {
	      for (l = m; l <= n; l++)
		{
		  REAL dT = (EVEN(n+l) ? dTRANS(n,l,m) : -dTRANS(n,l,m));
		  
		  Mout(n,m) += dT*Min(l,m);
		}
	      
	      for (; l < p; l++)
		Mout(n,m) += dTRANS(l,n,m)*Min(l,m);
	    }
	}
    }
}


/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
// Apply the translation operator to transform the top level of MP coeeffs 
// Assuming translation based upon p-1 levels already done.
void
CTransCoeff::incTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose)
{
  assert(Min.getOrder() >= m_p && Mout.getOrder() == m_p-1);
  Mout.setOrder(m_p);

  for (int n = 0; n < m_p-1; n++)
    for (int m = 0; m <= n; m++)
      {
	int l = m_p-1;
	if (tpose)
	  {
	      REAL T = (EVEN(n+l) ? TRANS(l,n,m) : -TRANS(l,n,m));
	      Mout(n,m) += T*Min(l,m);
	  }
	else
	  Mout(n,m) += TRANS(l,n,m)*Min(l,m);
      }

  int n = m_p-1;
  for (int m = 0; m <= n; m++)
    {
      int l;
      Mout(n,m) = 0.0;
      if (tpose)
	{
	  for (l = m; l <= n; l++)
	    Mout(n,m) += TRANS(n,l,m)*Min(l,m);
	  
	  for (; l < m_p; l++)
	    {
	      REAL T = (EVEN(n+l) ? TRANS(l,n,m) : -TRANS(l,n,m));
	      Mout(n,m) += T*Min(l,m);
	    }
	}
      else
	{
	  for (l = m; l <= n; l++)
	    {
	      REAL T = (EVEN(n+l) ? TRANS(n,l,m) : -TRANS(n,l,m));
	      Mout(n,m) += T*Min(l,m);
	    }
	  
	  for (; l < m_p; l++)
	    Mout(n,m) += TRANS(l,n,m)*Min(l,m);
	}
    }  
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::computeCoeff()
{
  TRANS(0,0,0) = m_exkid*m_K[0];
  m_exkid *= m_ir;
  if (m_bGrad)
    dTRANS(0,0,0) = -m_exkid*m_K[1];
  m_exkid *= CMCoeff::RS;

  for (int l = 1; l < 2*m_p-1; l++)
    {
      TRANS(l,0,0) = m_exkid*m_K[l];
      m_exkid *= m_ir;
      if (m_bGrad)
	dTRANS(l,0,0) = m_exkid*(l*m_K[l] - (2*l+1)*m_K[l+1]);
      m_exkid *= CMCoeff::RS;
    }

  computeCoeff_(m_T);

  if (m_bGrad)
    computeCoeff_(m_dT);
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::computeCoeff_(vector<REAL**> & U)
{
  for (int l = 1; l < 2*m_p-2; l++)
    U[l][1][0] = -(BETA(l-1,0)*U[l-1][0][0] +
		   ALPHA(l,0)*U[l+1][0][0])/ALPHA(0,0);
  
  for (int n = 1; n < m_p-1; n++)
    for (int l = n+1; l < 2*m_p-2-n; l++)
      U[l][n+1][0] = -(BETA(l-1,0)*U[l-1][n][0] +
		       BETA(n-1,0)*U[l][n-1][0] +
		       ALPHA(l,0)*U[l+1][n][0])/ALPHA(n,0);  
  for (int m = 1; m < m_p; m++)
    {
      for (int l = m; l < 2*m_p-1-m; l++)
	U[l][m][m] = -(DELTA(l,-m)*U[l-1][m-1][m-1] +
		       GAMMA(l+1,m-1)*U[l+1][m-1][m-1])/GAMMA(m,-m);
      
      for (int n = m; n < m_p-1; n++)
	for (int l = n+1; l < 2*m_p-2-n; l++)
	  U[l][n+1][m] = -(BETA(l-1,m)*U[l-1][n][m] +
			   BETA(n-1,m)*U[l][n-1][m] +
			   ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::computeIncCoeff()
{
  for (int l = 2*m_p-3; l < 2*m_p-1; l++)
    {
      TRANS(l,0,0) = m_exkid*m_K[l];
      m_exkid *= m_ir;
      if (m_bGrad)
	dTRANS(l,0,0) = m_exkid*(l*m_K[l] - (2*l+1)*m_K[l+1]);
      m_exkid *=  CMCoeff::RS;
    }

  computeIncCoeff_(m_T);

  if (m_bGrad)
    computeIncCoeff_(m_dT);
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::computeIncCoeff_(vector<REAL**> & U)
{
  if (m_p == 2)
    U[1][1][0] = -(BETA(0,0)*U[0][0][0] +
		   ALPHA(1,0)*U[2][0][0])/ALPHA(0,0);
  else
    for (int l = 2*m_p-4; (l < 2*m_p-2 && l > 0); l++)
      U[l][1][0] = -(BETA(l-1,0)*U[l-1][0][0] +
		     ALPHA(l,0)*U[l+1][0][0])/ALPHA(0,0);
  
  for (int n = 1; n < m_p-2; n++)
    for (int l = 2*m_p-4-n; l < 2*m_p-2-n; l++)
      U[l][n+1][0] = -(BETA(l-1,0)*U[l-1][n][0] +
		       BETA(n-1,0)*U[l][n-1][0] +
		       ALPHA(l,0)*U[l+1][n][0])/ALPHA(n,0);

  int n = m_p - 2;
  for (int l = n+1; (l < 2*m_p-2-n && n > 0); l++)
    U[l][n+1][0] = -(BETA(l-1,0)*U[l-1][n][0] +
		     BETA(n-1,0)*U[l][n-1][0] +
		     ALPHA(l,0)*U[l+1][n][0])/ALPHA(n,0);
  
  for (int m = 1; m < m_p-1; m++)
    {
      for (int l = 2*m_p-3-m; l < 2*m_p-1-m; l++)
	U[l][m][m] = -(DELTA(l,-m)*U[l-1][m-1][m-1] +
		       GAMMA(l+1,m-1)*U[l+1][m-1][m-1])/GAMMA(m,-m);
      
      for (int n = m; n < m_p-2; n++)
	for (int l = 2*m_p-4-n; l < 2*m_p-2-n; l++)
	  U[l][n+1][m] = -(BETA(l-1,m)*U[l-1][n][m] +
			   BETA(n-1,m)*U[l][n-1][m] +
			   ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
      int n = m_p-2;
      for (int l = n+1; l < 2*m_p-2-n; l++)
	  U[l][n+1][m] = -(BETA(l-1,m)*U[l-1][n][m] +
			   BETA(n-1,m)*U[l][n-1][m] +
			   ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
    }

  int m = m_p-1;
  for (int l = m; l < 2*m_p-1-m; l++)
    U[l][m][m] = -(DELTA(l,-m)*U[l-1][m-1][m-1] +
		   GAMMA(l+1,m-1)*U[l+1][m-1][m-1])/GAMMA(m,-m);
  
  for (int n = m; n < m_p-1; n++)
    for (int l = n+1; l < 2*m_p-2-n; l++)
      U[l][n+1][m] = -(BETA(l-1,m)*U[l-1][n][m] +
		       BETA(n-1,m)*U[l][n-1][m] +
		       ALPHA(l,m)*U[l+1][n][m])/ALPHA(n,m);
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

REAL
CTransCoeff::computeError(const CMCoeff & tM1, const CMCoeff & tM2, int p)
{
  assert(p <= m_p);
  assert(tM1.getOrder() >= p && tM2.getOrder() >= p);

  REAL err = 0.0;

  int l = p-1;
  for (int n = 0; n < p; n++)
    {
      err += TRANS(l,n,0)*tM1(n,0).real()*tM2(l,0).real();	
      for (int m = 1; m <= n; m++)
	err += 2*TRANS(l,n,m)*(tM1(n,m).real()*tM2(l,m).real() +
				 tM1(n,m).imag()*tM2(l,m).imag());
    }

  int n = p-1;
  for (l = 0; l < p-1; l++)
    {
      REAL T = (EVEN(n+l) ? TRANS(n,l,0) : -TRANS(n,l,0));
      err += T*tM1(n,0).real()*tM2(l,0).real();
    }

  for (int m = 1; m <= n; m++)
    for (l = m; l < p-1; l++)
      {
	REAL T = (EVEN(n+l) ? TRANS(n,l,m) : -TRANS(n,l,m));
	err += 2*T*(tM1(n,m).real()*tM2(l,m).real() +
		      tM1(n,m).imag()*tM2(l,m).imag());
      }
  
  return err;
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::outputTrans(int p) const
{
  cout << "****TRANS COEFFICIENTS****" << endl;
  for (int m = 0; m < p; m++)
    {
      cout << "\t---m = " << m << "---" << endl;
      for (int l = m; l < p; l++)
	{
	  for (int n = m; n <= l; n++)
	    {
	      REAL r = fabs(TRANS(l,n,m))>1e-15 ? TRANS(l,n,m) : 0;
	      cout << r << " | ";
	    }
	  cout << endl;
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::outputdTrans(int p) const
{
  cout << "****dTRANS COEFFICIENTS****" << endl;
  for (int m = 0; m < p; m++)
    {
      cout << "\t---m = " << m << "---" << endl;
      for (int l = m; l < p; l++)
	{
	  for (int n = m; n <= l; n++)
	    {
	      REAL r = fabs(dTRANS(l,n,m))>1e-15 ? dTRANS(l,n,m) : 0;
	      cout << r << " | ";
	    }
	  cout << endl;
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CTransCoeff::exportMat(double ** T)
{
  int ind[m_p];
  ind[0] = 1;
  for (int j = 1; j < m_p; j++)
    ind[j] = ind[j-1] + 2*j-1;

  for (int n = 0; n < m_p; n++)
    for (int m = -n; m <= n; m++) 
      {
	int j = ind[n] + (m+n);
	for (int l = abs(m); l < m_p; l++)
	  {
	    int k = ind[l] + (l+m);
	    if (l > n)
	      T[j][k] = TRANS(l,n,abs(m));
	    else
	      T[j][k] = (EVEN(n+l) ? TRANS(n,l,abs(m)) : -TRANS(n,l,abs(m)));
	  }
      }
}
