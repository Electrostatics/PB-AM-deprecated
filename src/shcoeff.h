#ifndef _SHCOEFF_H_
#define _SHCOEFF_H_

#include <complex>
#include <vector>
#include "util.h"
#include "mcoeff.h"

using namespace std;

class CSHCoeff : public CMCoeff
{
public:
  static void init(REAL rs, REAL kappa);
  CSHCoeff(int p = 0, int res = N_POLES);
  void reset(REAL theta, REAL phi, int p);

  void inc();
  static void besselk(REAL K[], int n, REAL val);
  static void incBesselk(REAL K[], int p, REAL val); 
  static void besseli(REAL I[], int n, REAL val);
  static void specialSH(REAL SH[], int n, REAL val);

private:
  static REAL CONST1[2*N_POLES][2*N_POLES], CONST2[2*N_POLES][2*N_POLES], 
    CONST3[2*N_POLES][2*N_POLES], CONST4[2*N_POLES], CONST5[2*N_POLES],
    CONST6[2*N_POLES];


  void legendre();
  void incLegendre();
  
  vector<REAL> m_P;

  REAL m_xval, m_negterm, m_sqroot, m_sqrtterm;
  vector<Complex> m_cis;
  REAL m_theta, m_phi;
};

#endif

  
