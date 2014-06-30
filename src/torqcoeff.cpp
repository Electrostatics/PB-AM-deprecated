#include "torqcoeff.h"
#include "shcoeff.h"

CTorqCoeff::CTorqCoeff(const vector<REAL> & charges, const vector<CPnt> & pos, 
		       int p, REAL rad) 
{
  assert(p > 0);
  reset(p);

  REAL I[p], Ix[p], Iy[p], Iz[p];
  for (int i = 0; i < pos.size(); i++)
    {
      CSpPnt spos = CartToSph(pos[i]);
      CSHCoeff SH(p);
      SH.reset(spos.theta(), spos.phi(), p);
            
      if (rad == 0)
	CSHCoeff::besseli(I, p, CMCoeff::KAPPA*spos.rho());
      else
	CSHCoeff::besseli(I, p, 0.0);

      CPnt k = charges[i]*pos[i];
      REAL r = spos.rho()*CMCoeff::IRS;
      for (int n = 0; n < p; n++)
	{
	  Ix[n] = I[n]*k.x();
	  Iy[n] = I[n]*k.y();
	  Iz[n] = I[n]*k.z();
	  k *= r;
	}

      m_M[0] += (SH * Ix);
      m_M[1] += (SH * Iy);
      m_M[2] += (SH * Iz);
    } 
}
