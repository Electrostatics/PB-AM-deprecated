#ifndef _SHCOEFF_H_
#define _SHCOEFF_H_

#include <complex>
#include <vector>
#include "util.h"
#include "mcoeff.h"

using namespace std;

//!  The CSHCoeff class 
/*! The CSHCoeff class contains information about spherical
 harmonics. Contains computations for EQ(1) of Lotan 2006
 as well as bessel function calculations for equations 2 and 3 */
class CSHCoeff : public CMCoeff
{
public:
  static void init(REAL rs, REAL kappa);
	
	//!  The CSHCoeff class constructor
	/*!
	 Initialize a CSHCoeff object. 
	 \param p is an int that is a minimal number of poles. Default is zero.
	 \param res is an in that describes the total number of residuals (poles) for the expansion. Default is N_POLES
	 \return an object of the CSHCoeff class. */
  CSHCoeff(int p = 0, int res = N_POLES);
  void reset(REAL theta, REAL phi, int p);
	
	//!  The CSHCoeff inc function
	/*! Function to increase SH for NPOLES=m_p+1 */	
  void inc();
	
	//!  The CSHCoeff besselk function
	/*! Function to create the i modified spherical bessel function. 		
	 kn(z) = sqrt(pi/(2z)) * K n+0.5 (z)
	 
	 There	is a recursion relationship as follows:
	 k n+1(z) = k n(z) + ( k n-1(z) * z^2 )/( (2n+1)*(2n-1) )
	 
	 \param K a vector of floating points to store besselK functions of size p
	 \param p an integer indicating the number of bessel functions to compute
	 \param val a floating point with the value at which to evaluate the bessel
	 function.  Given as z in the Lotan 2006 paper. */
  static void besselk(REAL K[], int n, REAL val);
  static void incBesselk(REAL K[], int p, REAL val); 
	
	//!  The CSHCoeff besseli function
	/*! Function to create the i modified spherical bessel function. 		
	 in(z) = sqrt(pi/(2z)) * I n+0.5 (z)
	 
	 There	is a recursion relationship as follows:
	 i n+1(z) = (2n+1)*(2n+3)*( i n-1(z) - i n(z) )/ z^2
	 
	 But the alternate method used is:
	 
	 i n(z) = 1 + sum_(j=1)^L t_j^n(z^2/2)
	 where: t_j^n(y) = (1/j)*t_(j-1)^n(y)*(y/(2n+2j+3))
	 
	 \param I a vector of floating points to store besseli functions of size p
	 \param p an integer indicating the number of bessel functions to compute
	 \param val a floating point with the value at which to evaluate the bessel
	 function.  Given as z in the Lotan 2006 paper. */
  static void besseli(REAL I[], int n, REAL val);
	
	
  static void specialSH(REAL SH[], int n, REAL val);
	
private:
  static REAL CONST1[2*N_POLES][2*N_POLES];		//!< (2l-1)/(l-m) for use in legendre computation
	static REAL CONST2[2*N_POLES][2*N_POLES];		//!< (l+m-1)/(l-m) for use in legendre computation
	static REAL CONST3[2*N_POLES][2*N_POLES];		//!< sqrt((n-m)!/(n+m)!) in EQ1, Lotan 2006
	static REAL CONST4[2*N_POLES];							//!< 1.0/((2*n-1)*(2*n-3)), Used for besselK recursion EQ3 in Lotan 2006
	//	static REAL CONST5[2*N_POLES];						// unused
	static REAL CONST6[2*N_POLES];							//!< (2l-1)!! double factorial, for use in legendre recursion
	
	//!  The CSHCoeff legendre function
	/*!
	 Function to create legendre polynomials for a given theta and phi.  
	 Created with p poles and with recursion functions. 
	 */
  void legendre();
	
	//!  The CSHCoeff incLegendre function
	/*!
	 Function to increase legendre polynomials for a given theta and phi.  
	 */
  void incLegendre();
  
  vector<REAL> m_P;														//!< Vector of legendre polynomials of NPOLES * NPOLES 
	
  REAL m_xval;																//!< xval a floating point of cos(theta)
	REAL m_negterm;															//!< Stored negative term of legendre polynomial
	REAL m_sqroot;															//!< store sqrt( 1 - cos(theta) )
	REAL m_sqrtterm;														//!< Related to m_sqroot above
  vector<Complex> m_cis;											//!< Complex number = exp( i*m*phi ) in EQ(1)
  REAL m_theta;																//!< Theta angle for spherical harmonics comp.
	REAL m_phi;																	//!< Phi angle for spherical harmonics comp.
}; // end CSHCoeff

#endif


