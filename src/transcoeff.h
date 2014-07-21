#ifndef _TRANSCOEFF_H_
#define _TRANSCOEFF_H_

#include "mcoeff.h"
#include "shcoeff.h"

//!  The transCoeff class
/*!	The class that contains all details for a translation coefficients, as
	described by the appendix of Lotan 2006.   */
class CTransCoeff
{
 public:
 
//!  The transCoeff class constructor
/*!	Constructing a transcoeff class, using input of whether or 
		not to compute the gradient
		\param bGrad is a boolean indicating of whether or not to compute gradients
						the default is yes  */ 
  CTransCoeff(bool bGrad = true);
	
//! TransCoeff reset function
/*!  Function to reset the translation coefficients and reallocate space
			\param rho a floating point with the new desired distance, d
			\param p an int containing the new number of poles */
  void reset(REAL rho, int p);

//! TransCoeff initConstants function
/*!  Function to initialize constants used in trans. coeff recursion */
  static void initConstants();

//! TransCoeff translate function
/*!  Function to apply the tranlation operator to the MP coefficents
			\param Min matrix coefficient input
			\param Mout matrix coefficient output
			\param tpose boolean to indicate whether or not to transpose matrix */
  void translate(const CMCoeff & Min, CMCoeff & Mout, bool tpose)
    { translate(Min, Mout, m_p, tpose); }

//! TransCoeff translate function
/*!  Function to apply the tranlation operator to the MP coefficents
			\param Min matrix coefficient input
			\param Mout matrix coefficient output
			\param p an integer indicating the number of poles
			\param tpose boolean to indicate whether or not to transpose matrix */
  void translate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose);
	
//! TransCoeff dTranslate function
/*!  Function to apply the derivative tranlation operator to the MP coefficents
			\param Min matrix coefficient input
			\param Mout matrix coefficient output
			\param tpose boolean to indicate whether or not to transpose matrix */
  void dTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose)
    { dTranslate(Min, Mout, m_p, tpose); }

//! TransCoeff dTranslate function
/*!  Function to apply the derivative tranlation operator to the MP coefficents
			\param Min matrix coefficient input
			\param Mout matrix coefficient output
			\param p an integer indicating the number of poles
			\param tpose boolean to indicate whether or not to transpose matrix */
  void dTranslate(const CMCoeff & Min, CMCoeff & Mout, int p, bool tpose);

//! TransCoeff dTranslate function
/*! Apply the translation operator to transform the top level of MP coeeffs 
assuming translation based upon p-1 levels already done. */
  void incTranslate(const CMCoeff & Min, CMCoeff & Mout, bool tpose);

//! TransCoeff computeError function
/*!  Function to compute the error between two matrix coefficients
			\param Min matrix coefficient input
			\param Mout matrix coefficient output
			\param p an integer indicating the number of poles */
  REAL computeError(const CMCoeff & tM1, const CMCoeff & tM2, int p);

//! Get the number of poles in the transfer coefficient object
  int getOrder() const
    { return m_p; }

  void incOrder();
  void decOrder();

// Print out functions
  void outputTrans(int p) const;
  void outputdTrans(int p) const;

// Saving and undoing tranlation matrices
  void saveUndo();
  void undo();

//! TransCoeff exportMat function
/*!  Export the translation matrix as a 2X2 matrix */  
  void exportMat(double ** T);

 private:
 
//! TransCoeff allocate function
/*!  Function to allocate space in T and dT matrices */
  void allocate();
//! TransCoeff deallocate function
/*!  Function to deallocate space in T and dT matrices */
  void deallocate();
//! TransCoeff reallocate function
/*!  Function to reallocate space in T and dT matrices
			calls either de or allocate according to m_p  
			\param p an int of desired number of poles */
  void reallocate(int p);
	
//! TransCoeff initParams function
/*!  Function to initialize parameters
			\param d is a floating point number of the magnitude of the 
							distance for translation */	
  void initParams(REAL d);

//! TransCoeff computeCoeff function
/*!  Function to compute translation coefficients.  Related to section of Appendix A.1
			of Lotan 2006.  */
  void computeCoeff(); 	
//! TransCoeff computeCoeff_ function
/*!  Function to compute translation coefficients for l!=0.  Related to section of Appendix A.1
			of Lotan 2006.  */
  void computeCoeff_(vector<REAL**> & U);

//! TransCoeff computeIncCoeff function
/*!  Function to compute addition translation coefficients when the 
			number of poles is increased.  */
  void computeIncCoeff(); 
//! TransCoeff computeIncCoeff_ function
/*!  Function to compute addition translation coefficients for l!=0 when the 
			number of poles is increased.  */
  void computeIncCoeff_(vector<REAL**> & U);

  static REAL & ALPHA(int n, int m) { return m_alpha[n][m]; }
  static REAL & BETA(int n, int m) { return m_beta[n][m]; }
  static REAL & GAMMA(int n, int m) { return m_gamma[n][N_POLES-1+m]; }
  static REAL & DELTA(int n, int m) { return m_delta[n][N_POLES-1+m]; }
  static bool & EVEN(int n) { return m_even[n]; } 

  static REAL m_alpha[N_POLES*2][N_POLES];			//!< Given as alpha(n,m) in Lotan, 2006 (eq 1.9)
  static REAL m_beta[N_POLES*2][N_POLES];				//!< Given as beta(n,m) in Lotan, 2006 (eq 1.9)
  static REAL m_gamma[N_POLES*2][2*N_POLES-1];	//!< Given as eta(n,m) in Lotan, 2006 (eq 1.9)
  static REAL m_delta[N_POLES*2][2*N_POLES-1];	//!< Given as mu(n,m) in Lotan, 2006 (eq 1.9)
  static bool m_even[4*N_POLES];								//!< Flag for whether number is even or not

//! TransCoeff TRANS function
/*!  Function to return the l,n,m coefficient of the translation coefficient array */	  
  REAL & TRANS(int l,int n, int m) { return m_T[l][n][m]; }
//! TransCoeff dTRANS function
/*!  Function to return the l,n,m coefficient of the translation derivative coefficient array */	 
  REAL & dTRANS(int l, int n, int m) { return m_dT[l][n][m]; }
//! TransCoeff TRANS function
/*!  Function to return the l,n,m coefficient of the translation coefficient array */	 
  REAL TRANS(int l,int n, int m) const { return m_T[l][n][m]; }
//! TransCoeff dTRANS function
/*!  Function to return the l,n,m coefficient of the translation derivative coefficient array */	 
  REAL dTRANS(int l, int n, int m) const { return m_dT[l][n][m]; }

  vector<REAL**> m_T;							//!< vector of translation coefficients
	vector<REAL**> m_TU;						//!< vector of translation coefficients stored temporarily incase of undo
	vector<REAL**> m_dT;						//!< vector of translation derivatives
	vector<REAL**> m_dTU;						//!< vector of translation derivatives stored temporarily incase of undo
  REAL m_exkid;										//!< exponential of k/d = exp(-k*d)/d, part of addition theorem (EQ5)
	REAL m_ir;											//!< Inverse of distance = 1/d
	REAL m_d;												//!< Distance of interest
	REAL m_kd;											//!< Distance times the inverse debye length (k * d)	
	REAL m_rbu[4];									//!< vector of exkid, ir, d and kd stored temporarily incase of undo
	
  REAL m_K[N_POLES*2];						//!< vector of modified spherical bessel functions of the second kind
	REAL m_KU[N_POLES*2];						//!< vector of MSBF stored temporarily incase of undo
  bool m_bGrad;										//!< A boolean indicating whether or not to compute the gradient

  int m_p;												//!< An int of the poles
	int m_pU;												//!< An int of the poles stored temporarily incase of undo
};  // end CTransCoeff

/////////////////////////////////////
////// Inline functions

//!  CTransCoeff incOrder
/*!	Increase the number of poles of the calculation.  Memory must
			be allocated, the bessel functions increased and the coefficients
			added too
*/
inline void
CTransCoeff::incOrder()
{
  allocate();
  m_p++;
  CSHCoeff::incBesselk(m_K, 2*m_p-1, m_kd);
  CSHCoeff::incBesselk(m_K, 2*m_p, m_kd);
  computeIncCoeff();
}

//!  CTransCoeff decOrder
/*!	Decrease the number of poles of the calculation.  Memory must
			be deallocated etc
*/
inline void
CTransCoeff::decOrder()
{
  deallocate();
  m_p--;
}

//!  CTransCoeff saveUndo
/*!	Save the translation coefficients and derivatives
		as an undone move
*/
inline void
CTransCoeff::saveUndo()
{
  if (m_bGrad)							// if computing the gradient, save both T and derivatives
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
    m_KU[i] = m_K[i];			// Save K's 
	
  m_pU = m_p;							// Save the number of poles
} // end saveUndo

//!  CTransCoeff undo
/*!	Revert to old tranlation coefficients using saved 
			components
*/
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
} // end undo

#endif
