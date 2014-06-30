#ifndef _ROTCOEFF_H
#define _ROTCOEFF_H

#include "mcoeff.h"
#include "tricoeff.h"
#include "shcoeff.h"

class CRotCoeff
{
 public:
  CRotCoeff(bool bGrad = true);

  static void initConstants();

  void reset(REAL theta, REAL phi, REAL xi, int p);
  void reset(const CPnt & axis, REAL ang, int p);
  void reset(const CQuat & Q, int p);

  void rotate(const CTriCoeff & Tin, CTriCoeff & Tout, bool bFor)
    { rotate(Tin, Tout, 1, m_p, bFor); }
  void rotate(const CTriCoeff & Tin, CTriCoeff & Tout, int p1, int p2, 
	      bool bFor);
  void rotate(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
    { rotate(Min, Mout, 1, m_p, bFor); }
  void rotate(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2, 
	      bool bFor);
  void dRotateT(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
    { dRotateT(Min, Mout, 1, m_p, bFor); }
  void dRotateT(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2, 
		bool bFor);
  void dRotateP(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
    { dRotateP(Min, Mout, 1, m_p, bFor); }
  void dRotateP(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2,
		bool bFor);

  void incOrder();
  void decOrder()
    { deallocate(), m_p--; }

  void saveUndo();
  void undo();
    
  void outputRot(int p) const;
  void outputdRot(int p) const;

  bool isSingular() const
    { return m_bSing; }
  void getParams(REAL & st, REAL & ct, REAL & sp, REAL & cp) const
    { st = m_sint; ct = m_cost; sp = m_exphi.imag(); cp = m_exphi.real(); } 
  
 private:
  void allocate();
  void deallocate();
  void reallocate(int p);

  void initParams(REAL theta, REAL phi, REAL xi);

  void dRotateTSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2);
  void dRotatePSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2);

  void computeCoeff();
  void computeIncCoeff();
  void computeGradCoeff();
  void computeIncGradCoeff();

  const Complex & ROT(int n, int s, int m) const;
  Complex & ROT(int n, int s, int m); 
  const Complex & dROT(int n, int s, int m) const;
  Complex & dROT(int n, int s, int m); 

  static Complex derv(const Complex & c, int s)
    { return Complex(-s*c.imag(), s*c.real()); }
  static void computeQCoeff();
  static REAL & ZETA(int n, int m) { return m_zeta[n][m]; }
  static REAL & ETA(int n, int m) { return m_eta[n][2*N_POLES-1+m]; }
  
  static REAL m_zeta[N_POLES*2][2*N_POLES-1];
  static REAL m_eta[N_POLES*2][4*N_POLES-1];
  static REAL m_Q[2*N_POLES-1][N_POLES][2];

  vector<Complex**> m_R, m_RU, m_dR, m_dRU;
  REAL m_sint, m_cost, m_cott;
  Complex m_exphi, m_exiphi, m_exxi;
  REAL m_r1, m_r2, m_r3, m_dr1, m_dr2, m_dr3;
  CSHCoeff m_SH;
  REAL m_theta, m_phi, m_xi;
  bool m_bSing, m_bGrad;
  CQuat m_Quat, m_QuatU;
  REAL m_rbu[12];
  Complex m_cbu[3];
  int m_p, m_pU;
};

inline const Complex &  
CRotCoeff::ROT(int n, int s, int m) const
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_R[n][n+s][m]; 
}

inline Complex &  
CRotCoeff::ROT(int n, int s, int m) 
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_R[n][n+s][m]; 
}

inline const Complex &  
CRotCoeff::dROT(int n, int s, int m) const
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_dR[n][n+s][m]; 
}

inline Complex &  
CRotCoeff::dROT(int n, int s, int m) 
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_dR[n][n+s][m]; 
}

inline void 
CRotCoeff::rotate(const CTriCoeff & Tin, CTriCoeff & Tout, int p1, int p2, 
		  bool bFor)
{ 
  rotate(Tin[0], Tout[0], p1, p2, bFor);
  rotate(Tin[1], Tout[1], p1, p2, bFor);
  rotate(Tin[2], Tout[2], p1, p2, bFor);
}

inline void
CRotCoeff::saveUndo()
{
  if (m_bGrad)
    {
      for (int n = 0; n < 2*m_p-1; n++)
	for (int s = -(2*m_p-1); s < 2*m_p- 1; s++)
	  for (int m = 0; m < 2*m_p - 1; m++)
	    {
	      m_RU[n][n+s][m] = m_R[n][n+s][m];
	      m_dRU[n][n+s][m] = m_dR[n][n+s][m];
	    }
    }
  else
    {
      for (int n = 0; n < 2*m_p-1; n++)
	for (int s = -(2*m_p-1); s < 2*m_p- 1; s++)
	  for (int m = 0; m < 2*m_p - 1; m++)
	    m_RU[n][n+s][m] = m_R[n][n+s][m];
    }

  m_rbu[0] = m_sint;
  m_rbu[1] = m_cost;
  m_rbu[2] = m_cott;
  m_rbu[3]  = m_r1;
  m_rbu[4]  = m_r2;
  m_rbu[5]  = m_r3;
  m_rbu[6]  = m_dr1;
  m_rbu[7]  = m_dr2;
  m_rbu[8]  = m_dr3;
  m_rbu[9] = m_theta;
  m_rbu[10] = m_phi;
  m_rbu[11] = m_xi;
  
  m_cbu[0] = m_exphi;
  m_cbu[1] = m_exiphi;
  m_cbu[2] = m_exxi;

  m_SH.saveUndo();
  m_QuatU = m_Quat;

  m_pU = m_p;
}

inline void
CRotCoeff::undo()
{
  if (m_bGrad)
    {
      for (int n = 0; n < 2*m_p-1; n++)
	for (int s = -(2*m_p-1); s < 2*m_p- 1; s++)
	  for (int m = 0; m < 2*m_p - 1; m++)
	    m_R[n][n+s][m] = m_RU[n][n+s][m];
    }
  else
    {
      for (int n = 0; n < 2*m_p-1; n++)
	for (int s = -(2*m_p-1); s < 2*m_p- 1; s++)
	  for (int m = 0; m < 2*m_p - 1; m++)
	    m_R[n][n+s][m] = m_RU[n][n+s][m];
    }

  m_sint = m_rbu[0];
  m_cost = m_rbu[1];
  m_cott = m_rbu[2];
  m_r1 = m_rbu[3];
  m_r2 = m_rbu[4];
  m_r3 = m_rbu[5];
  m_dr1 = m_rbu[6];
  m_dr2 = m_rbu[7];
  m_dr3 = m_rbu[8];
  m_theta = m_rbu[9];
  m_phi = m_rbu[10];
  m_xi = m_rbu[11];
  
  m_exphi = m_cbu[0];
  m_exiphi = m_cbu[1];
  m_exxi = m_cbu[2];

  m_SH.undo();

  m_Quat = m_QuatU;
  m_p = m_pU;
}


#endif
