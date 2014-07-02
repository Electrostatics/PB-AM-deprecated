#include "rotcoeff.h"

#define EPS_SIN_THETA (1e-12)

REAL CRotCoeff::m_zeta[N_POLES*2][2*N_POLES-1];
REAL CRotCoeff::m_eta[N_POLES*2][4*N_POLES-1];
REAL CRotCoeff::m_Q[2*N_POLES-1][N_POLES][2];

/******************************************************************/
/******************************************************************//**
* Initialize the rotation coefficients (includes coeff equs from 
* 2006 paper
******************************************************************/

void
CRotCoeff::initConstants()
{
  for (int n = 0; n < 2*N_POLES; n++)
    for (int m = 0; m <= n; m++)
      ZETA(n,m) = sqrt((REAL)(n+m+1)*(n-m+1)/((2*n+1)*(2*n+3))); //!< Given as a(n,m) in Lotan, 2006 (eq 1.2)

  for (int n = 0; n < 2*N_POLES; n++)
    for (int m = 0; m <= n; m++)
      {
				ETA(n,m) = sqrt((REAL)(n-m-1)*(n-m)/((2*n+1)*(2*n-1))); //!< Given as b(n,m) in Lotan, 2006 (eq 1.2)
				if (m != 0)
					ETA(n,-m) = -sqrt((REAL)(n+m-1)*(n+m)/((2*n+1)*(2*n-1))); //!< Given as b(n,m) in Lotan, 2006 (eq 1.2)
      }

  for (int n = 0; n < 2*N_POLES; n++)
    computeQCoeff();
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

CRotCoeff::CRotCoeff(bool bGrad) : m_bSing(false),
				   m_SH(0, 1), m_bGrad(bGrad)
{
  m_R.resize(1);
  m_R[0] = new Complex*[1];
  m_R[0][0] = new Complex[1];

 if (m_bGrad)
   {
     m_dR.resize(1);
     m_dR[0] = new Complex*[1];
     m_dR[0][0] = new Complex[1];
   }
 
 m_p = 1;
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CRotCoeff::allocate()
{
  assert(m_R.size() < 2*N_POLES-1);

  int p = m_R.size();
  m_R.resize(p+2);

  m_R[p] = new Complex*[2*p+1];
  m_R[p+1] = new Complex*[2*p+3];
  for (int j = 0; j < 2*p+1; j++)
    m_R[p][j] = new Complex[p+1];
  for (int j = 0; j < 2*p+3; j++)
    m_R[p+1][j] = new Complex[p+2];
   
  if (m_bGrad)
    {
      m_dR.resize(p+2);
      
      m_dR[p] = new Complex*[2*p+1];
      m_dR[p+1] = new Complex*[2*p+3];
      for (int j = 0; j < 2*p+1; j++)
	m_dR[p][j] = new Complex[1+p];
      for (int j = 0; j < 2*p+3; j++)
	m_dR[p+1][j] = new Complex[p+2];
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CRotCoeff::deallocate()
{
  assert(m_R.size() > 1);
  int p = m_R.size();

  for (int j = 0; j < 2*p - 1; j++)
    delete[] m_R[p-1][j];
  for (int j = 0; j < 2*p - 3; j++)
    delete[] m_R[p-2][j];
  delete[] m_R[p-1];
  delete[] m_R[p-2];
  m_R.resize(p-2);

  if (m_bGrad)
    {
      for (int j = 0; j < 2*p-1; j++)
	delete[] m_dR[p-1][j];
      for (int j = 0; j < 2*p-3; j++)
	delete[] m_dR[p-2][j];
      delete[] m_dR[p-1];
      delete[] m_dR[p-2];
      m_dR.resize(p-2);
    }
}

/******************************************************************/
/******************************************************************//**
*  a function to reallocate enough space for the coefficient for 
each pole
******************************************************************/

void
CRotCoeff::reallocate(int p)
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
*  reset rotation coefficient given an angle for theta, phi and xi
		and a number of poles
		\param theta, a real number of angle theta
		\param phi, a real number for phi angle in rad
		\xi a 
******************************************************************/

void 
CRotCoeff::reset(REAL theta, REAL phi, REAL xi, int p)
{
  assert(p > 0 && p <= N_POLES);

  initParams(theta, phi, xi);
  reallocate(p);
  m_SH.reset(m_theta, m_phi, 2*p-1);

  m_p = p;  
  if (!isSingular())
    {
      computeCoeff();
      if (m_bGrad)
	computeGradCoeff();
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CRotCoeff::reset(const CPnt & axis, REAL ang, int p)
{
  assert(p > 0 && p <= N_POLES);
	
  CPnt ax = axis.normalize();
  REAL c = cos(ang);
  REAL s = sin(ang);
  REAL v = 1.0 - c;
  
  REAL kxkz = ax.x()*ax.z()*v;
  REAL kykz = ax.y()*ax.z()*v;
	
  REAL theta = acos(ax.z()*ax.z()*v + c);
  REAL phi, xi;
  if (fabs(theta) > 1e-12) 
	{ 
		phi = atan2(kykz + ax.x()*s, kxkz - ax.y()*s);
		xi = atan2(kykz - ax.x()*s, kxkz + ax.y()*s);
	}
  else
	{ 
		theta = 0.0;
		if (fabs(ax.x()) > 1e-12 || fabs(ax.y()) > 1e-12)
		{
			phi = atan2(ax.x(), -ax.y());
			xi = atan2(-ax.x(), ax.y());
		}
		else
		{
			phi = M_PI/2.0;
			xi = -M_PI/2.0;
		}
	}
	
  reset(theta, phi, xi, p);
  m_Quat = CQuat(ax, ang); 
}

/******************************************************************/
/******************************************************************//**
*  reset: Resets the rotation coefficients for a given number of poles
	\param Q a quaternion for rotation
	\param p an int describing the number of poles
******************************************************************/

void
CRotCoeff::reset(const CQuat & Q, int p)
{
  assert(p > 0 && p <= N_POLES);
	
  REAL theta = acos(1.0 - 2*Q.x()*Q.x() - 2*Q.y()*Q.y());
  REAL phi, xi;
  if (fabs(theta) > 1e-12) 
	{ 
		phi = atan2(Q.y()*Q.z() + Q.w()*Q.x(), Q.x()*Q.z() - Q.w()*Q.y());
		xi = atan2(Q.y()*Q.z() - Q.w()*Q.x(), Q.x()*Q.z() + Q.w()*Q.y());
	}
  else
	{
		theta = 0.0;
		phi = M_PI/2.0;
		xi = -M_PI/2.0;
	}
	
  reset(theta, phi, xi, p);
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CRotCoeff::initParams(REAL theta, REAL phi, REAL xi)
{
  m_theta = theta; m_phi = phi; m_xi = xi;
	
  m_sint = sin(m_theta);
  if (m_sint < EPS_SIN_THETA)
	{
		m_bSing = true;
		if (m_theta < M_PI/2)
		{
			m_theta = 0.0;
			m_phi = 0.0;
			m_sint = 0.0;
			m_cost = 1.0;
			m_exphi = Complex(1.0,0.0);
		}
		else
		{
			m_theta = M_PI;
			m_phi = 0.0;
			m_sint = 0.0;
			m_cost = -1.0;
			m_exphi = Complex(1.0,0.0);
		}
	}
  else
	{
		m_bSing = false;
		m_cost = cos(m_theta);
		m_cott = m_cost/m_sint;
		m_exphi = Complex(cos(phi),sin(phi));
	}
	
  m_exiphi = conj(m_exphi);
  m_exxi = Complex(cos(xi),sin(xi));	  
}

/******************************************************************/
/******************************************************************//**
*  Apply the rotation operator to the MP coeffs
******************************************************************/
void
CRotCoeff::rotate(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2, 
		  bool bFor)
{
  assert(Min.getOrder() >= p2 && p2 <= N_POLES && p2 <= m_p);
  assert(p1 == 1 || Mout.getOrder() == p1 - 1);
  assert(p1 != 0 && p1 <= p2);

  if (isSingular())
    {
      if (m_cost == -1.0)
	{
	  Mout.setOrder(p2);
	  for (int n = p1 - 1; n < p2; n++)
	    for (int m = 0; m <= n; m++) 
	      Mout(n,m) = (n % 2 == 0 ? Min(n,-m) : -Min(n,-m));
	}
      else
	{
	  if (p1 == 1)
	    Mout.setOrder(0);

	  for (int p = p1; p <= p2; p++)
	    Mout.copy_p(Min, p);
	}
    }
  else
    {
      Mout.setOrder(p2);
      if (bFor)
	{
	  for (int n = p1 - 1; n < p2; n++)
	    for (int m = 0; m <= n; m++)
	      {
		Mout(n,m) = ROT(n,0,m)*Min(n,0);
		
		for (int s = 1; s <= n; s++)
		  Mout(n,m) += (ROT(n,s,m)*Min(n,s) + ROT(n,-s,m)*Min(n,-s));
	      }
	}
      else
	{
	  for (int n = p1 - 1; n < p2; n++)
	    for (int s = 0; s <= n; s++)
	      {
		Mout(n,s) = ROT(n,-s,0)*Min(n,0);
		
		for (int m = 1; m <= n; m++)
		  Mout(n,s) += conj(ROT(n,s,m))*Min(n,m) + 
		    ROT(n,-s,m)*Min(n,-m);
	      }
	}
    }
}

/******************************************************************/
/******************************************************************//**
* Apply the derivative of the rotation operator with respect to THETA
    to the MP coeffs
******************************************************************/

void
CRotCoeff::dRotateT(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2,
		    bool bFor)
{
  if (!m_bGrad)
    {
      cout << "dRdt: This operator has no derivatives!!!" << endl;
      return;
    }

  assert(Min.getOrder() >= p2 && p2 <= N_POLES && p2 <= m_p);
  assert(p1 == 1 || Mout.getOrder() == p1 - 1);
  assert(p1 != 0 && p1 <= p2);

  Mout.setOrder(p2);
  
  if (isSingular())
    {
      dRotateTSing(Min, Mout, p1, p2);
      if (!bFor)
	Mout.recip();
    }
  else
    {
      if (bFor)
	{
	  for (int n = p1 - 1; n < p2; n++)
	    for (int m = 0; m <= n; m++)
	      {
		Mout(n,m) = dROT(n,0,m)*Min(n,0);
		
		for (int s = 1; s <= n; s++)
		  Mout(n,m) += (dROT(n,s,m)*Min(n,s) + dROT(n,-s,m)*Min(n,-s));
	      }
	}
      else
	{
	  for (int n =  p1 - 1; n < p2; n++)
	    for (int s = 0; s <= n; s++)
	      {
		Mout(n,s) = dROT(n,-s,0)*Min(n,0);
		
		for (int m = 1; m <= n; m++)
		  Mout(n,s) += conj(dROT(n,s,m))*Min(n,m) +
		    dROT(n,-s,m)*Min(n,-m);
	      }
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs
******************************************************************/
void
CRotCoeff::dRotateP(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2,
		    bool bFor)
{
  if (!m_bGrad)
    {
      cout << "dRdp: This operator has no derivatives!!!" << endl;
      return;
    }

  assert(Min.getOrder() >= p2 && p2 <= N_POLES && p2 <= m_p);
  assert(p1 == 1 || Mout.getOrder() == p1 - 1);
  assert(p1 != 0 && p1 <= p2);

  Mout.setOrder(p2);

  if (isSingular())
    {
      dRotatePSing(Min, Mout, p1, p2);
      if (!bFor && m_cost == 1.0)
	Mout.recip();
    }
  else
    {
      if (bFor)
	{
	  for (int n = p1 - 1; n < p2; n++)
	    for (int m = 0; m <= n; m++)
	      {
		Mout(n,m) = Complex(0.0,0.0);

		for (int s = 1; s <= n; s++)
		  Mout(n,m) -=  (ROT(n,s,m)*derv(Min(n,s),s) + 
				 ROT(n,-s,m)*derv(Min(n,-s),-s)); 
	      }
	}
      else
	{
	  for (int n = p1 - 1 ; n < p2; n++)
	    for (int s = 0; s <= n; s++)
	      {
		Mout(n,s) = ROT(n,-s,0)*derv(Min(n,0),s);
		
		for (int m = 1; m <= n; m++)
		  Mout(n,s) += conj(ROT(n,s,m))*derv(Min(n,m),s) + 
		    ROT(n,-s,m)*derv(Min(n,-m),s);
	      }
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs in the singular case sin(theta) = 0.0
******************************************************************/
void 
CRotCoeff::dRotateTSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2)
{
  assert(Min.getOrder() >= p2 && p2 <= N_POLES && p2 <= m_p);
  assert(p1 == 1 || Mout.getOrder() == p1 - 1);  
  assert(p1 != 0 && p1 <= p2);

  if (p1 == 1)
    {
      Mout(0,0) = Complex(0.0,0.0);
      p1++;
    }

  if (m_cost == 1.0)
    {
      for (int n = p1 - 1; n < p2; n++)
	{
	  Mout(n,0) = Complex(2.0*m_Q[n][0][1]*Min(n,1).real(),0.0);
	  for (int m = 1; m < n; m++)
	    Mout(n,m) = m_Q[n][m][0]*Min(n,m-1) + m_Q[n][m][1]*Min(n,m+1); 

	  Mout(n,n) = m_Q[n][n][0]*Min(n,n-1); 
	}
    }
  else
    {
      REAL s = -1.0;
      for (int n = p1 - 1; n < p2; n++, s = -s)
	{
	  Mout(n,0) = Complex(2.0*s*m_Q[n][0][1]*Min(n,1).real(),0.0);
	  for (int m = 1; m < n; m++)
	    Mout(n,m) = s*(m_Q[n][m][0]*Min(n,-m+1) + 
			   m_Q[n][m][1]*Min(n,-m-1)); 

	  Mout(n,n) = s*m_Q[n][n][0]*Min(n,-n+1);
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  Apply the derivative of the rotation operator with respect to PHI 
			to the MP coeffs in the singular case sin(theta) = 0.0
******************************************************************/
void 
CRotCoeff::dRotatePSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2)
{
  assert(Min.getOrder() >= p2 && p2 <= N_POLES && p2 <= m_p);
  assert(p1 == 1 || Mout.getOrder() == p1 - 1); 
  assert(p1 != 0 && p1 <= p2);

  if (p1 == 1)
    {
      Mout(0,0) = Complex(0.0,0.0);
      p1++;
    }
  
  if (m_cost == 1.0)
    {
      for (int n = p1 - 1; n < p2; n++)
	{
	  Mout(n,0) = Complex(2.0*m_Q[n][0][1]*Min(n,1).imag(), 0.0); 
	  for (int m = 1; m < n; m++)
	    Mout(n,m) = Complex(0.0,m_Q[n][m][0])*Min(n,m-1) - 
	      Complex(0.0,m_Q[n][m][1])*Min(n,m+1); 

	  Mout(n,n) = Complex(0.0, m_Q[n][n][0])*Min(n,n-1); 
	}
    }
  else
    {
      REAL s = 1.0;
      for (int n = p1 - 1; n < p2; n++, s = -s)
	{
	  Mout(n,0) = Complex(s*2*m_Q[n][0][1]*Min(n,1).imag(), 0.0); 
	  for (int m = 1; m < n; m++)
	    Mout(n,m) =  s*(Complex(0.0,m_Q[n][m][1])*Min(n,-m-1) -
			    Complex(0.0,m_Q[n][m][0])*Min(n,-m+1));
	  
	  Mout(n,n) = Complex(0.0,-s*m_Q[n][n][0])*Min(n,-n+1);
	}
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
void
CRotCoeff::computeCoeff()
{
  const CSHCoeff & SH = m_SH;
  Complex dum1, dum2, dum3;

  ROT(0,0,0) = SH(0,0);
  for (int n = 1; n < 2*m_p-1; n++)
    {
      ROT(n,0,0) = SH(n,0);
      for (int s = 1; s <  n; s++)
	{	 
	  ROT(n,s,0) = SH(n,-s);
	  ROT(n,-s,0) = SH(n,s);
	}

      ROT(n,n,0) = SH(n,-n);
      ROT(n,-n,0) = SH(n,n);
    }

  m_r1 = -0.5*(1 + m_cost);
  m_r2 = 0.5*(1 - m_cost);
  m_r3 = -m_sint;

  for (int m = 0; m < m_p-1; m++)
    for (int n = m+2; n < 2*m_p-m-1; n++)
      for (int s = -n+1; s < n; s++)
	{
	  dum1 = (m_r1*ETA(n,s-1))*ROT(n,s-1,m);
	  dum2 = (m_r2*ETA(n,-s-1))*ROT(n,s+1,m);
	  dum3 = (m_r3*ZETA(n-1,abs(s)))*ROT(n,s,m);
	  ROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi+dum2*m_exphi+dum3)/ETA(n,m);
	}
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
void
CRotCoeff::computeGradCoeff()
{
  const CSHCoeff & SH = m_SH;
  Complex dum1, dum2, dum3;

  dROT(0,0,0) = Complex(0.0,0.0);
  for (int n = 1; n < 2*m_p-1; n++)
    {
      dROT(n,0,0) = -sqrt((REAL)n*(n+1))*m_exiphi*SH(n,1);
      for (int s = 1; s <  n; s++)
	{	 
	  dROT(n,-s,0) = (s*m_cott)*SH(n,s) -
	    sqrt((REAL)(n-s)*(n+s+1))*m_exiphi*SH(n,s+1);
	  dROT(n,s,0) = conj(dROT(n,-s,0));
	}
      
      dROT(n,-n,0) = (n*m_cott)*SH(n,n);
      dROT(n,n,0) = conj(dROT(n,-n,0));
    }
  
  m_dr1 = 0.5*m_sint;
  m_dr2 = m_dr1;
  m_dr3 = -m_cost;
  for (int m = 0; m < m_p-1; m++)
    for (int n = m+2; n < 2*m_p-m-1; n++)
      for (int s = -n+1; s < n; s++)
	{
	  dum1 = ETA(n,s-1)*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
	  dum2 = ETA(n,-s-1)*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
	  dum3 = ZETA(n-1,abs(s))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
	  dROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi+dum2*m_exphi+dum3)/ETA(n,m);
	}
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/

void
CRotCoeff::computeIncCoeff()
{
  const CSHCoeff & SH = m_SH;
  Complex dum1, dum2, dum3;

  for (int n = 2*m_p-3; n < 2*m_p-1; n++)
    {
      ROT(n,0,0) = SH(n,0);
      for (int s = 1; s < n; s++)
	{	 
	  ROT(n,s,0) = SH(n,-s);
	  ROT(n,-s,0) = SH(n,s);
	}

      ROT(n,n,0) = SH(n,-n);
      ROT(n,-n,0) = SH(n,n);
    }

  for (int m = 0; m < m_p-2; m++)
    for (int n = 2*m_p-m-3; n < 2*m_p-m-1; n++)
      for (int s = -n+1; s < n; s++)
	{
	  dum1 = (m_r1*ETA(n,s-1))*ROT(n,s-1,m);
	  dum2 = (m_r2*ETA(n,-s-1))*ROT(n,s+1,m);
	  dum3 = (m_r3*ZETA(n-1,abs(s)))*ROT(n,s,m);
	  ROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi+dum2*m_exphi+dum3)/ETA(n,m);
        }

  int m = m_p-2;
  for (int n = m+2; n < 2*m_p-m-1; n++)
    for (int s = -n+1; s < n; s++)
      {
	dum1 = (m_r1*ETA(n,s-1))*ROT(n,s-1,m);
	dum2 = (m_r2*ETA(n,-s-1))*ROT(n,s+1,m);
	dum3 = (m_r3*ZETA(n-1,abs(s)))*ROT(n,s,m);
	ROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi + dum2*m_exphi + dum3)/ETA(n,m);
      }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
void
CRotCoeff::computeIncGradCoeff()
{ 
  const CSHCoeff & SH = m_SH;
  Complex dum1, dum2, dum3;

  for (int n = 2*m_p-3; n < 2*m_p-1; n++)
    {
      dROT(n,0,0) = -sqrt((REAL)n*(n+1))*m_exiphi*SH(n,1);
      for (int s = 1; s < n; s++)
	{	 
	  dROT(n,-s,0) = (s*m_cott)*SH(n,s) -
	    sqrt((REAL)(n-s)*(n+s+1))*m_exiphi*SH(n,s+1);
	  dROT(n,s,0) = conj(dROT(n,-s,0));
	}
      
      dROT(n,-n,0) = (n*m_cott)*SH(n,n);
	  dROT(n,n,0) = conj(dROT(n,-n,0));
    }
  
  for (int m = 0; m < m_p-2; m++)
    for (int n = 2*m_p-m-3; n < 2*m_p-m-1; n++)
      for (int s = -n+1; s < n; s++)
	{
	  dum1 = ETA(n,s-1)*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
	  dum2 = ETA(n,-s-1)*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
	  dum3 = ZETA(n-1,abs(s))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
	  dROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi+dum2*m_exphi+dum3)/ETA(n,m);
        }
  
  int m = m_p-2;
  for (int n = m+2; n < 2*m_p-m-1; n++)
    for (int s = -n+1; s < n; s++)
      {
	dum1 = ETA(n,s-1)*(m_r1*dROT(n,s-1,m) + m_dr1*ROT(n,s-1,m));
	dum2 = ETA(n,-s-1)*(m_r2*dROT(n,s+1,m) + m_dr2*ROT(n,s+1,m));
	dum3 = ZETA(n-1,abs(s))*(m_r3*dROT(n,s,m) + m_dr3*ROT(n,s,m));
	dROT(n-1,s,m+1) = m_exxi*(dum1*m_exiphi+dum2*m_exphi+dum3)/ETA(n,m);
      }
}

/******************************************************************/
/******************************************************************//**
*  Computing Q coefficients, 
******************************************************************/
void
CRotCoeff::computeQCoeff()
{
  REAL SH[2*N_POLES-1];
  CSHCoeff::specialSH(SH, 2*N_POLES-1, 1.0);

  REAL temp;
  for (int n = 1; n < 2*N_POLES-1; n++)
    {
      m_Q[n][0][0] = 0.0;
      m_Q[n][0][1] = SH[n];
    }

  for (int n = 2; n < 2*N_POLES-1; n++)
    {
      m_Q[n-1][1][0] = (ETA(n,-1)*m_Q[n][0][1] +  ZETA(n-1,0))/ETA(n,0);
      
      if (n > 2)
				m_Q[n-1][1][1] = (ETA(n,1)/ETA(n,0))*m_Q[n][0][1];
      else
				m_Q[n-1][1][1] = 0.0;
    }

  for (int m = 1; m < N_POLES-1; m++)
    for (int n = m+2; n < 2*N_POLES-m-1; n++)
		{
			m_Q[n-1][m+1][0] = (ETA(n,m-1)*m_Q[n][m][0] +  ZETA(n-1,m))/ETA(n,m);

			if (n > m+2)
				m_Q[n-1][m+1][1] =  (ETA(n,m+1)/ETA(n,m))*m_Q[n][m][1];
			else
				m_Q[n-1][m+1][1] = 0.0;
      }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
void
CRotCoeff::incOrder()
{
  allocate();
  m_p++;

  if (!isSingular())
    {
      while (m_SH.getOrder() < 2*m_p-1)
	m_SH.inc();

      computeIncCoeff();
      if (m_bGrad)
	computeIncGradCoeff();
    }
}

/******************************************************************/
/******************************************************************//**
*  
******************************************************************/
void
CRotCoeff::outputRot(int p) const
{
  if (p > m_p)
    p = m_p;

  cout << "****ROT COEFFICIENTS****" << endl;
  for (int m = 0; m < m_p; m++)
    {
      cout << "\t---m = " << m << "---" << endl;
      for (int n = m; n < m_p; n++)
	{
	  for (int s = -n; s <= n; s++)
	    {
	      REAL r = fabs(ROT(n,s,m).real())>1e-15 ? 
		ROT(n,s,m).real() : 0;
	      REAL im = fabs(ROT(n,s,m).imag())>1e-15 ? 
		ROT(n,s,m).imag() : 0;
	      cout << "(" << r << "," << im << ") | ";
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
CRotCoeff::outputdRot(int p) const
{
  if (p > m_p)
    p = m_p;

  cout << "****dROT COEFFICIENTS****" << endl;
  for (int m = 0; m < m_p; m++)
    {
      cout << "\t---m = " << m << "---" << endl;
      for (int n = m; n < m_p; n++)
	{
	  for (int s = -n; s <= n; s++)
	    {
	      REAL r = fabs(dROT(n,s,m).real())>1e-15 ? 
		dROT(n,s,m).real() : 0;
	      REAL im = fabs(dROT(n,s,m).imag())>1e-15 ? 
		dROT(n,s,m).imag() : 0;
	      cout << "(" << r << "," << im << ") | ";
	    }
	  cout << endl;
	}
    }
}
