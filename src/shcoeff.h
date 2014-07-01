#ifndef _SHCOEFF_H_
#define _SHCOEFF_H_

#include <complex>
#include <vector>
#include "util.h"
#include "mcoeff.h"

using namespace std;
/******************************************************************//**
* # File: shcoeff.h
* #
* # Date: June 2014
* #
* # Description: This is the header file for the class SHCoeff
* #								SH = spherical harmonics
* #
* # Author: Lotan, Felberg
* #
* # Copyright ( c )
* #
******************************************************************/


	//!  The CSHCoeff class 
/*!
		The CSHCoeff class contains information about spherical
		harmonics.
*/
class CSHCoeff : public CMCoeff
{
public:
  static void init(REAL rs, REAL kappa);
	
		//!  The CSHCoeff class constructor
/*!
		Initialize a CSHCoeff object.
		\param p is an int that is a minimal number of poles. Default is zero.
		\param res is an in that describes the total number of residuals (poles) for the expansion. Default is N_POLES
		\return an object of the CSHCoeff class
*/
  CSHCoeff(int p = 0, int res = N_POLES);
  void reset(REAL theta, REAL phi, int p);
	
  void inc();
  static void besselk(REAL K[], int n, REAL val);
  static void incBesselk(REAL K[], int p, REAL val); 
  static void besseli(REAL I[], int n, REAL val);
  static void specialSH(REAL SH[], int n, REAL val);
	
private:
  static REAL CONST1[2*N_POLES][2*N_POLES];
	static REAL CONST2[2*N_POLES][2*N_POLES]; 
	static REAL CONST3[2*N_POLES][2*N_POLES];		//<! sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
	static REAL CONST4[2*N_POLES];
	static REAL CONST5[2*N_POLES];
	static REAL CONST6[2*N_POLES];
	
  void legendre();
  void incLegendre();
  
  vector<REAL> m_P;
	
  REAL m_xval;
	REAL m_negterm;
	REAL m_sqroot;
	REAL m_sqrtterm;
  vector<Complex> m_cis;
  REAL m_theta;																//<! Theta angle for spherical harmonics comp.
	REAL m_phi;																	//<! Phi angle for spherical harmonics comp.
}; // end CSHCoeff

#endif

  
