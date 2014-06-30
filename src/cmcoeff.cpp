#include "cmcoeff.h"

CMCoeff
operator+(const CMCoeff & M1,  const CMCoeff & M2)
{
  assert(M1.m_p == M2.m_p);
  CMCoeff S(M1.m_p);

  for (int n = 0; n < m_p; n++)
    for (int m = 0; m < n; m++)
      S(n,m) = M1(n,m) + M2(n,m);
}

REAL inprod(const CMCoeff & M1,  const CMCoeff & M2)
{
  assert(M1.m_p == M2.m_p);
  
  REAL sum1(0.0), sum2(0.0);
  for (int n = 0; n < p; n++)
    {
      for (int m = 1; m <= n; m++)
	sum1 += M1(n,m).real()*M2(n,m).real() -
	  M1(n,m).imag()*M2(n,m).imag();
	
      sum2 += M1(n,0).real()*M2(n,0).real();
    }

  return 2*sum1 + sum2;
}

CMCoeff::CMCoeff(REAL charges[], const CSpPnt pos[], int num,
		 int p, REAL kappa)
{
  REAL I[p];
  for (int i = 0; i < num; i++)
    {
      Complex SH[N_POLES][N_POLES];
      CMPE::sphericalHarmonics(SH, pos[i].theta(), pos[i].phi(), p);
      
      CMPE::besseli(I, p, kappa*pos[i].rho());
      REAL r = 1.0;
      for (int n = 0; n < p; n++)
	{
	  REAL k = I[n]*r*charges[i];
	  for (int m = 0; m <= n; m++)
	    m_M[n][m] += k*SH[n][m];
	   
	  r *= pos[i].rho();
	}
    }  
}

CMCoeff::CMCoeff(REAL theta, REAL phi, int p)
{
  REAL P[N_POLES][N_POLES];
  CMCoeff::legendre(P, p, cos(theta));
  
  Complex cis[p];
  for (int m = 0; m <= p; m++)
    cis[m] = Complex(cos(m*phi), sin(m*phi));
		     
  for (int n = 0; n < p; n++)
    for (int m = 0; m <= n; m++)
      m_M[n][m] = (P[n][m]*CONST3[n][m])*cis[m];
}

static void 
CMCoeff::legendre(REAL P[N_POLES][N_POLES], int p, REAL xval) const
{                                       
  int m, l;                          
  REAL negterm, oddfact, nextodd, sqroot, sqrtterm;                
  
  negterm = 1.0;
  oddfact = 1.0;                                                   
  nextodd = 1.0;                                                   
  sqroot = sqrt(1.0 - xval*xval);                                  
  sqrtterm = 1.0;                                                  
  for(int m = 0; m < p ; m++)
    {                                           
      P[m][m] = negterm*oddfact*sqrtterm;                        
      negterm *= -1.0;                                             
      oddfact *= nextodd;                                          
      nextodd += 2.0;                                              
      sqrtterm *= sqroot;                                          
      if(m < p-1)
	{                                                
	  P[m+1][m] = xval * (REAL)(2*m+1) * P[m][m];           
	  for(n = m+2; n < p; n++)
	    P[n][m] = xval*CONST1*P[n-1][m] - CONST2*P[n-2][m];
	}                                                        
    }                                                                
}

static void 
CMCoeff::besselk(REAL K[], int p, REAL val) const
{
  if (p > 0)
    K[0] = 1.0;

  if (p > 1)
    K[1] = 1.0 + val; 

  REAL valsqr = val*val;
  for (int n = 2; n < p; n++)
    K[n] = K[n-1] + valsqr*K[n-2]*CONST4;
}

/*
static void 
CMCoeff::besseli(REAL I[], int p, REAL val) const
{
  for (int n = 0; n < p; n++)
     I[n] = 1.0;
  
  if (val != 0.0)
    {
      REAL z = 0.5*val*val;
      for (int n = 0; n < p; n++)
	{
	  REAL t = z/(2*n+3);
	  for (int j = 1; j <= 20; j++)
	    {
	      I[n] += t;
	      t *= (z/((j+1)*(2*(n+j)+3)));
	    }
        }
    }
}
*/

static void 
CMCoeff::besseli(REAL I[], int p, REAL val) const
{ 
  if (val != 0.0)
    {
      long double ival = 1.0/val;
      long double ivalsqr = ival*ival;
      long double J[p];
      
      if (p > 0)
        J[0] = sinh((long double)val)*ival;
      
      if (p > 1)
        J[1] = 3*ivalsqr*(cosh((long double)val) - J[0]);
  
      for (int j = 1; j < p-1; j++)
	J[j+1] = CONST5[j]*(J[j-1]-J[j])*ivalsqr;

      for (int n = 0; n < p; n++)
	I[n] = (REAL)J[n];
    }
  else
    for (int n = 0; n < p; n++)
      I[n] = 1.0;
}

static void
CMCoeff::init()
{
  REAL temp[2*N_POLES];
  temp[0] = 1.0;
  for (int i = 1; i < 2*N_POLES; i++)
    temp[i] = temp[i-1]*sqrt((REAL)i);
  
  for (int n = 0; n < N_POLES; n++)
    {
      for (int m = 0; m <= n; m++)
	{
	  CONST1[n][m] = (2*n-1)/(REAL)(n-m);
	  CONST2[n][m] = (n+m-1)/(REAL)(n-m);
	  CONST3[n][m] = temp[n-m]/temp[n+m];
      }

      CONST4[n] = 1.0/((2*n-1)*(2*n-3));
      CONST5[n] = (REAL)((2*n+1)*(2*n+3)); 
    }
