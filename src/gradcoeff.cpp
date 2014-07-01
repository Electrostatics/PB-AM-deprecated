#include <cmath>

#include "gradcoeff.h"

/******************************************************************//**
* # File: gradcoeff.cpp
* #
* # Date: June 2014
* #
* # Description: This file contains the class GradCoeff and its functions
* #
* # Author: Lotan, Felberg
* #
* # Copyright ( c )
* #
******************************************************************/

/******************************************************************/
/******************************************************************//**
* Creating Grad Coefficent class
* 
******************************************************************/

CGradCoeff::CGradCoeff(const CMCoeff & M, const CSpPnt & c) 
{ 
  m_M[0] = M.dMdr(c.rho()); 
  m_M[1] = M.dMdt(c.theta(), c.phi());
  m_M[2] = M.dMdp();
}

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

