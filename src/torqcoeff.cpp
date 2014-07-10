#include "torqcoeff.h"
#include "shcoeff.h"

/******************************************************************/
/******************************************************************//**
* Initialize the torque coefficients.
******************************************************************/
CTorqCoeff::CTorqCoeff(const vector<REAL> & charges, const vector<CPnt> & pos, 
		       int p, REAL rad) 
{
  assert(p > 0);
  reset(p);
	
  REAL I[p], Ix[p], Iy[p], Iz[p];
  for (int i = 0; i < pos.size(); i++)			// for each charge
	{
		CSpPnt spos = CartToSph(pos[i]);				// convert cartesian to spherical coords
		CSHCoeff SH(p);													// Create a spherical harmonics series for this point
		SH.reset(spos.theta(), spos.phi(), p);	// Reset it WRT to current charge's theta and phi angles
		
		if (rad == 0)														// If the rad is set to zero
			CSHCoeff::besseli(I, p, CMCoeff::KAPPA*spos.rho());	// Create a besseli I for p poles and kappa*r of charge
		else
			CSHCoeff::besseli(I, p, 0.0);					// else, make besseli, I of current position of p poles
		
		CPnt k = charges[i]*pos[i];							// Multiply charge magnutude by X, Y and Z coords
		REAL r = spos.rho()*CMCoeff::IRS;				// Scale radius of charge by scaling factor 1/Ravg
		for (int n = 0; n < p; n++)							// For each level of poles
		{
			Ix[n] = I[n]*k.x();										
			Iy[n] = I[n]*k.y();
			Iz[n] = I[n]*k.z();
			k *= r;
		}
		
		m_M[0] += (SH * Ix);										// The M is a 3 component MCoeff class object
		m_M[1] += (SH * Iy);										// that is a component of the torque calculation
		m_M[2] += (SH * Iz);										// It is equation 42/43 of Lotan 2006 without the gamma part
	} 
} // end CTorqCoeff constructor
