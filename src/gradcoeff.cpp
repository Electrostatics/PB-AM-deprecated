#include <cmath>

#include "gradcoeff.h"

/******************************************************************//**
* Creating Grad Coefficent class
******************************************************************/

CGradCoeff::CGradCoeff(const CMCoeff & M, const CSpPnt & c) 
{ 
  m_M[0] = M.dMdr(c.rho());							// Computing dMatrix/d rad
  m_M[1] = M.dMdt(c.theta(), c.phi());  // Computing dMatrix/d theta
  m_M[2] = M.dMdp();										// Computing dMatrix/d phi
}


/******************************************************************//**
* Printing out each partial derivative of the re-expansion operator
******************************************************************/
ostream & 
operator<<(ostream & out, const CGradCoeff & G)
{
  cout << "\t---dMdr---" << endl;
  cout << G.m_M[0] << endl;

  cout << "\t---dMdt---" << endl;
  cout << G.m_M[1] << endl;

  cout << "\t---dMdp---" << endl;
  cout << G.m_M[2] << endl;
  
  return out;
}

