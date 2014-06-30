#ifndef _TRANSCOEFF_H_
#define _TRANSCOEFF_H_

#include "mcoeff.h"
#include "shcoeff.h"

class CTransCoeff
{
 public:
  CTransCoeff(bool bGrad = true);
  void reset(REAL rho, int p);
  static void initConstants();

  void translate(const CMCoeff & Min, CMCoeff & Mout, bool tpose)
    { translate(Min, Mout, m_p, tpose); }
  void translate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose);
  void dTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose)
    { dTranslate(Min, Mout, m_p, tpose); }
  void dTranslate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose);
  void incTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose);

  REAL computeError(const CMCoeff & tM1, const CMCoeff & tM2, int p);
  int getOrder() const
    { return m_p; }

  void incOrder();
  void decOrder();

  void outputTrans(int p) const;
  void outputdTrans(int p) const;

  void saveUndo();
  void undo();
  
  void exportMat(double ** T);

 private:
  void allocate();
  void deallocate();
  void reallocate(int p);
  void initParams(REAL d);
  void computeCoeff(); 
  void computeCoeff_(vector<REAL**> & U);
  void computeIncCoeff(); 
  void computeIncCoeff_(vector<REAL**> & U);

  static REAL & ALPHA(int n, int m) { return m_alpha[n][m]; }
  static REAL & BETA(int n, int m) { return m_beta[n][m]; }
  static REAL & GAMMA(int n, int m) { return m_gamma[n][N_POLES-1+m]; }
  static REAL & DELTA(int n, int m) { return m_delta[n][N_POLES-1+m]; }
  static bool & EVEN(int n) { return m_even[n]; } 

  static REAL m_alpha[N_POLES*2][N_POLES];
  static REAL m_beta[N_POLES*2][N_POLES];
  static REAL m_gamma[N_POLES*2][2*N_POLES-1];
  static REAL m_delta[N_POLES*2][2*N_POLES-1];
  static bool m_even[4*N_POLES];

  
  REAL & TRANS(int l,int n, int m) { return m_T[l][n][m]; }
  REAL & dTRANS(int l, int n, int m) { return m_dT[l][n][m]; }
  REAL TRANS(int l,int n, int m) const { return m_T[l][n][m]; }
  REAL dTRANS(int l, int n, int m) const { return m_dT[l][n][m]; }

  vector<REAL**> m_T, m_TU, m_dT, m_dTU;
  REAL m_exkid, m_ir, m_d, m_kd;
  REAL m_K[N_POLES*2], m_KU[N_POLES*2];
  bool m_bGrad;
  REAL m_rbu[4];
  int m_p, m_pU;
};

inline void
CTransCoeff::incOrder()
{
  allocate();
  m_p++;
  CSHCoeff::incBesselk(m_K, 2*m_p-1, m_kd);
  CSHCoeff::incBesselk(m_K, 2*m_p, m_kd);
  computeIncCoeff();
}

inline void
CTransCoeff::decOrder()
{
  deallocate();
  m_p--;
}

inline void
CTransCoeff::saveUndo()
{
  if (m_bGrad)
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    {
	      m_TU[l][n][m] = m_T[l][n][m];
	      m_dTU[l][n][m] = m_dT[l][n][m];
	    }
    }
  else
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    m_TU[l][n][m] = m_T[l][n][m];
    }

  m_rbu[0] = m_exkid;
  m_rbu[1] = m_ir;
  m_rbu[2] = m_d;
  m_rbu[3] = m_kd;

  for (int i = 0; i < 2*m_p; i++)
    m_KU[i] = m_K[i];

  m_pU = m_p;
}

inline void
CTransCoeff::undo()
{
  if (m_bGrad)
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    {
	      m_T[l][n][m] = m_TU[l][n][m];
	      m_dT[l][n][m] = m_dTU[l][n][m];
	    }
    }
  else
    {
      for (int l = 0; l < 2*m_p-1; l++)
	for (int n = 0; n < m_p; n++)
	  for (int m = 0; m <= n; m++)
	    m_T[l][n][m] = m_TU[l][n][m];
    }

  m_exkid = m_rbu[0];
  m_ir = m_rbu[1];
  m_d = m_rbu[2];
  m_kd = m_rbu[3];

  for (int i = 0; i < 2*m_p; i++)
    m_K[i] = m_KU[i];

  m_p = m_pU;
}







#endif
