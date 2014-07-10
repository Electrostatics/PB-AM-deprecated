#include "shcoeff.h"

REAL CSHCoeff::CONST1[2*N_POLES][2*N_POLES];
REAL CSHCoeff::CONST2[2*N_POLES][2*N_POLES];
REAL CSHCoeff::CONST3[2*N_POLES][2*N_POLES];
REAL CSHCoeff::CONST4[2*N_POLES];  
//REAL CSHCoeff::CONST5[2*N_POLES];
REAL CSHCoeff::CONST6[2*N_POLES];

/******************************************************************/
/******************************************************************//**
* Initialize the spherical harmonics coefficients
* Inputs:   rs = scaling factor ~ avg mol radius in system
*						kappa = inverse debye length
******************************************************************/
void
CSHCoeff::init(REAL rs, REAL kappa)
{
  REAL temp[4*N_POLES];									// temp vector used to calculate sqrt((n-m)!/(n+m)!)
  temp[0] = 1.0;
  for (int i = 1; i < 4*N_POLES; i++)
    temp[i] = temp[i-1]*sqrt((REAL)i);
  
  for (int n = 0; n < 2*N_POLES; n++)
	{
		for (int m = 0; m <= n; m++)
		{
			CONST1[n][m] = (2*n-1)/(REAL)(n-m); // (2l-1)/(l-m) for use in legendre computation
			CONST2[n][m] = (n+m-1)/(REAL)(n-m); // (l+m-1)/(l-m) for use in legendre computation
			CONST3[n][m] = temp[n-m]/temp[n+m];	// sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
		}
		
		CONST4[n] = 1.0/((2*n-1)*(2*n-3));		// Used for besselK recursion EQ3 in Lotan 2006
//		CONST5[n] = (REAL)((2*n+1)*(2*n+3)); 
	}
	
  CONST6[0] = 1.0; CONST6[1] = 1.0;			// (2l-1)!! double factorial, for use in legendre recursion
  for (int n = 2; n < 2*N_POLES; n++)
    CONST6[n] = CONST6[n-1]*(2*n - 1);
  
  IDX[0] = 0;
  for (int n = 1; n <= 2*N_POLES; n++)
    IDX[n] = IDX[n-1]+n;
	
  RS = rs;
  IRS = 1.0/RS;
  KAPPA = kappa;
}

/******************************************************************/
/******************************************************************//**
* Initialize SHCoeff class
******************************************************************/
CSHCoeff::CSHCoeff(int p, int res) : CMCoeff(p,res)
{ 
  assert(p <= res);
  m_P.reserve(IDX[res]);
  m_cis.reserve(res);

  m_P.resize(IDX[p]); 
  m_cis.resize(p); 
}

/******************************************************************/
/******************************************************************//**
* Reset the spherical harmonics for a given theta and phi
	\param theta a given theta angle for computing SH
	\param phi a given phi angle for computing SH
	\param p an int with the number of poles for computing SH
******************************************************************/
void
CSHCoeff::reset(REAL theta, REAL phi, int p) 
{
  m_P.resize(IDX[p]); 
  m_cis.resize(p);
  CMCoeff::reset(p);
	
  m_theta = theta; 
  m_phi = phi;
  m_xval = cos(theta);
  CSHCoeff::legendre();
  
  for (int m = 0; m < m_p; m++)
    m_cis[m] = Complex(cos(m*m_phi), sin(m*m_phi));
	
  for (int n = 0; n < m_p; n++)
    for (int m = 0, s = 1; m <= n; m++, s = -s)
      m_M[IDX[n]+m] = (s*m_P[IDX[n]+m]*CONST3[n][m])*m_cis[m];
}

/******************************************************************/
/******************************************************************//**
* Function to increase the number of poles and create a legendre 
		polynomial for m_p+1 poles.
******************************************************************/
void
CSHCoeff::inc()
{
  assert(m_p+1 <= 2*N_POLES);
	
  setOrder(m_p+1);
  m_cis.resize(m_p);
  m_P.resize(IDX[m_p]); 
	
  incLegendre();
	
	
  m_cis[m_p-1] = Complex(cos((m_p-1)*m_phi), sin((m_p-1)*m_phi)); // exp( i*m*phi )
	
  for (int m = 0, s = 1; m < m_p; m++, s = -s)
    m_M[IDX[m_p-1]+m] = (s*m_P[IDX[m_p-1]+m]*CONST3[m_p-1][m])*m_cis[m]; // Exactly EQ(1) of Lotan 2006
										// (-1)^m * Pn,m(cos(theta) * sqrt((n-m)!/(n+m)!)*
}

/******************************************************************/
/******************************************************************//**
* Function to compute the Associated Legendre Polynomial of 
	m_xval = cos(theta).  This uses the recursion functions for the 
	polynomials, which are as follows:
	
	Pl,l (x) = (-1)^l * (2l-1)!! * (1-x^2)^(l/2)
	Pl,l+1 (x) = x * (2l+1) * Pl,l(x)
	Pl,m (x) = x * (2l-1)/(l-m) * Pl-1,m(x) - (l+m-1)/(l-m) * Pl-2,m(x) 
******************************************************************/
void 
CSHCoeff::legendre() 
{                         
  m_negterm = 1.0;
  m_sqroot = sqrt(1.0 - m_xval*m_xval);                                  
  m_sqrtterm = 1.0;                                                  
  for(int m = 0; m < m_p ; m++)
	{                                           
		m_P[IDX[m]+m] = m_negterm*CONST6[m]*m_sqrtterm;                        
		m_negterm = -m_negterm;                                             
		m_sqrtterm *= m_sqroot;                                          
		if(m < m_p-1)
		{                                                
			m_P[IDX[m+1]+m] = m_xval * (REAL)(2*m+1) * m_P[IDX[m]+m];           
			for(int n = m+2; n < m_p; n++)
				m_P[IDX[n]+m] = m_xval*CONST1[n][m]*m_P[IDX[n-1]+m] - 
													CONST2[n][m]*m_P[IDX[n-2]+m];
		}                                                        
	}                                                                
}

/******************************************************************/
/******************************************************************//**
* Creating a special spherical harmonic, from a legendre polynomial,
using the recursive relationship, for m=1:

	Pl,1 (val) = val * (2l-1)/(l-1) * Pl-1,1(val) - l/(l-1) * Pl-2,1(val)
	
	and finally: SH[l] = -1*sqrt((l-1)!/(l+1)!)*Pl,1(val)

******************************************************************/
void 
CSHCoeff::specialSH(REAL SH[], int n, REAL val)
{
  SH[0] = 0;
  SH[1] = -1;
  SH[2] = -val * 3;    

  for (int l = 3 ; l < n; l++)                                
    SH[l] = (val*(2*l-1)*SH[l-1] - l*SH[l-2])/(l-1);  

  // This part converts the legendre polynomial to a spherical harmonic.
  for (int l = 1; l < n; l++)
    SH[l] *= (-CONST3[l][1]);
}


/******************************************************************/
/******************************************************************//**
* increase Legendre.  The function to increase the legendre polynomial.
 This uses the recursion functions for the 
	polynomials, which are as follows:
	
	Pl,l (x) = (-1)^l * (2l-1)!! * (1-x^2)^(l/2)
	Pl,l+1 (x) = x * (2l+1) * Pl,l(x)
	Pl,m (x) = x * (2l-1)/(l-m) * Pl-1,m(x) - (l+m-1)/(l-m) * Pl-2,m(x) 
******************************************************************/
void 
CSHCoeff::incLegendre() 
{                         
  int n = m_p-1, m;                          
 
  m_P[IDX[n]+n] = m_negterm*CONST6[n]*m_sqrtterm;                        
  m_negterm = -m_negterm;                                             
  m_sqrtterm *= m_sqroot;                                          

  m_P[IDX[n]+n-1] = m_xval * (REAL)(2*n-1) * m_P[IDX[n-1]+n-1];           
  for (int m = 0; m < m_p-2; m++)
    m_P[IDX[n]+m] = m_xval*CONST1[n][m]*m_P[IDX[n-1]+m] -    
      CONST2[n][m]*m_P[IDX[n-2]+m];
}

/******************************************************************/
/******************************************************************//**
* Function to create the k modified spherical bessel function.
******************************************************************/
void 
CSHCoeff::besselk(REAL K[], int p, REAL val) 
{
  if (p > 0)
    K[0] = 1.0;
	
  if (p > 1)
    K[1] = 1.0 + val; 
	
  REAL valsqr = val*val;
  for (int n = 2; n < p; n++)
    K[n] = K[n-1] + valsqr*K[n-2]*CONST4[n];
}

/******************************************************************/
/******************************************************************//**
* Function to increase the stored k modified spherical bessel function 
******************************************************************/
void 
CSHCoeff::incBesselk(REAL K[], int p, REAL val) 
{
  if (p == 1)
    K[0] = 1.0;
  else if (p == 2)
    K[1] = 1.0 + val; 
  else
	{
		REAL valsqr = val*val;
		K[p-1] = K[p-2] + valsqr*K[p-3]*CONST4[p-1];
	}
}

/******************************************************************/
/******************************************************************//**
* Function to create the i modified spherical bessel function. 
******************************************************************/
void 
CSHCoeff::besseli(REAL * I, int p, REAL val) 
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
	      t *= (z/((j+1)*(2*(n+j)+3)));		// EQ 1.15
	      if (t < 1e-20)
					break;
	    }
		}
	}
}
