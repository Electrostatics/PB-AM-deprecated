#ifndef _ROTCOEFF_H
#define _ROTCOEFF_H

#include "mcoeff.h"
#include "tricoeff.h"
#include "shcoeff.h"

//!  The rotation coefficient class
/*!
		The class that contains all details for rotation coefficients,
		which are described in the Appendix, section A.1 of Lotan 2006. 
*/

class CRotCoeff
{
public:

//!  The rotation coefficient class constructor
/*!
		This creates an object of the RotCoeff class, using input
		of the user.  Sets bSing to false and NPOLES to 1
		\param bGrad a boolean that identifies whether or not the
							gradients of rotation coefficients should be computed
							or not.  Default is yes.  */
  CRotCoeff(bool bGrad = true);

//!  The rotation coefficient initConstants function
/*! This function initializes the rotation coefficients
	as described in A.1 of Lotan 2006.  It initializes a and b
	constants as described in the appendix. Calls computeQCoeff */	
  static void initConstants();
	
  void reset(REAL theta, REAL phi, REAL xi, int p);
  void reset(const CPnt & axis, REAL ang, int p);
  void reset(const CQuat & Q, int p);

//! CRotCoeff rotate function
/*!  Function to apply the rotation operator to the tricoefficents
			\param Tin tri coefficient input
			\param Tout tri coefficient output
			\param bFor boolean to indicate to compute the original 
					coefficient or the conjugate transpose  */	
  void rotate(const CTriCoeff & Tin, CTriCoeff & Tout, bool bFor)
	{ rotate(Tin, Tout, 1, m_p, bFor); }
//! CRotCoeff rotate function
/*!  Function to apply the rotation operator to the tricoefficents
			\param Tin tri coefficient input
			\param Tout tri coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void rotate(const CTriCoeff & Tin, CTriCoeff & Tout, int p1, int p2, 
							bool bFor);							
//! CRotCoeff rotate function
/*!  Function to apply the rotation operator to the MP coefficents
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void rotate(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
	{ rotate(Min, Mout, 1, m_p, bFor); }
//! CRotCoeff rotate function
/*!  Function to apply the rotation operator to the MP coefficents
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */	
  void rotate(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2, 
							bool bFor);
//! CRotCoeff dRotateT function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
		with respect to THETA
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void dRotateT(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
	{ dRotateT(Min, Mout, 1, m_p, bFor); }
//! CRotCoeff dRotateT function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
		with respect to THETA
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void dRotateT(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2, 
								bool bFor);
//! CRotCoeff dRotateP function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
		with respect to PHI 
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the number of poles
			\param p2 an integer indicating the number of poles
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void dRotateP(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
	{ dRotateP(Min, Mout, 1, m_p, bFor); }
//! CRotCoeff dRotateP function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
	with respect to PHI 
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles
			\param bFor boolean to indicate to compute the original 
					coefficient (true) or the conjugate transpose  */
  void dRotateP(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2,
								bool bFor);
	
  void incOrder();
  void decOrder()
	{ deallocate(), m_p--; }

// Print out functions	
  void outputRot(int p) const;
  void outputdRot(int p) const;
	
// Saving and undoing rotation coefficient matrices		
  void saveUndo();
  void undo();

// Determine if the sint values approch 0		
  bool isSingular() const
	{ return m_bSing; }

// Return parameters of rotation coefficient matrix, sin theta, cos theta, etc.
  void getParams(REAL & st, REAL & ct, REAL & sp, REAL & cp) const
	{ st = m_sint; ct = m_cost; sp = m_exphi.imag(); cp = m_exphi.real(); } 
  
private:
//! RotCoeff allocate function
/*!  Function to allocate space in R and dR matrices */
  void allocate();
//! RotCoeff deallocate function
/*!  Function to deallocate space in R and dR matrices */
  void deallocate();
//! RotCoeff reallocate function
/*!  Function to reallocate space in R and dR matrices.
			will call either re or allocate depending on size
			of input WRT to npoles
			\param p an int of new number of poles   */
  void reallocate(int p);
	
//! RotCoeff initParams function
/*!  Function to initialize parameters
			\param theta is a floating point number of theta angle
			\param phi is a floating point number of phi angle
			\param xi is a floating point number of xi */
  void initParams(REAL theta, REAL phi, REAL xi);

//! CRotCoeff dRotateTSing function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
		with respect to THETA for sint -> 0, when near singlularity
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles */	
  void dRotateTSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2);
//! CRotCoeff dRotatePSing function
/*!  Function to apply the gradient of rotation operator to the MP coefficents
		with respect to PHI for sint -> 0, when near singlularity
			\param Min mp coefficient input
			\param Mout mp coefficient output
			\param p1 an integer indicating the lowest number of poles
			\param p2 an integer indicating the highest number of poles */	
  void dRotatePSing(const CMCoeff & Min, CMCoeff & Mout, int p1, int p2);

//! RotCoeff computeCoeff function
/*!  Function to compute rotation coefficients.  Related to section of Appendix A.1
			of Lotan 2006.  */	
  void computeCoeff();
//! RotCoeff computeIncCoeff function
/*!  Function to compute additional rotation coefficients when the 
			number of poles is increased.  */
  void computeIncCoeff();
//! RotCoeff computeGradCoeff function
/*!  Function to compute gradient rotation coefficients.  Related to section of Appendix A.1
			of Lotan 2006.  */	
  void computeGradCoeff();
//! RotCoeff computeIncGradCoeff function
/*!  Function to compute additional gradient rotation coefficients when the 
			number of poles is increased.  */
  void computeIncGradCoeff();

// Retrieving components of rotation and gradient coeff matrices	
  const Complex & ROT(int n, int s, int m) const;
  Complex & ROT(int n, int s, int m); 
  const Complex & dROT(int n, int s, int m) const;
  Complex & dROT(int n, int s, int m); 
	
  static Complex derv(const Complex & c, int s)
	{ return Complex(-s*c.imag(), s*c.real()); }
	
//! RotCoeff computeQCoeff function
/*!  Function to compute special spherical harmonics,
		used only when the function is near a singularity, e.g. sint -> 0 */
  static void computeQCoeff();

// Printing/Retrieving a and b coefficients	
  static REAL & ZETA(int n, int m) { return m_zeta[n][m]; }
  static REAL & ETA(int n, int m) { return m_eta[n][2*N_POLES-1+m]; }
  
  static REAL m_zeta[N_POLES*2][2*N_POLES-1];	//!< Given as a(n,m) in Lotan, 2006 (eq 1.2)
  static REAL m_eta[N_POLES*2][4*N_POLES-1];	//!< Given as b(n,m) in Lotan, 2006 (eq 1.2)
  static REAL m_Q[2*N_POLES-1][N_POLES][2];		//!< Q coefficients used in computing dR for singularities
	
  vector<Complex**> m_R;					//!< vector of rotation coefficients
	vector<Complex**> m_RU;					//!< vector of rotation coefficients stored temporarily incase of undo
	vector<Complex**> m_dR;					//!< vector of gradient rotation coeff
	vector<Complex**> m_dRU;				//!< vector of gradient rotation coeff stored temporarily incase of undo
  REAL m_sint;										//!< The sin of theta
	REAL m_cost;										//!< The cos of theta
	REAL m_cott;										//!< The cotangent of theta
  Complex m_exphi;								//!< The complex evaluation of e^(i*phi) = cos(phi) + i*sin(phi)
	Complex m_exiphi;								//!< The complex evaluation of e^(-i*phi) = cos(phi) - i*sin(phi) = conj(exphi)
	Complex m_exxi;									//!< The complex evaluation of e^(i*xi) = cos(xi) + i*sin(xi) 
  REAL m_r1;											//!< -0.5*(1 + m_cost) component of EQ 1.1
	REAL m_r2;											//!< 0.5*(1 - m_cost) component of EQ 1.1
	REAL m_r3;											//!< -sint component of EQ 1.1
	REAL m_dr1;											//!< 0.5*m_sint component of EQ 1.4
	REAL m_dr2;											//!< 0.5*m_sint component of EQ 1.4 
	REAL m_dr3;											//!< -m_cost component of EQ 1.4 
  CSHCoeff m_SH;									//!< A spherical harmonic object
  REAL m_theta;										//!< The angle theta in radians
	REAL m_phi;											//!< The angle phi in radians
	REAL m_xi;											//!< Evaluation of e^(i*xi) = cos(xi) + i*sin(xi)
  bool m_bSing;										//!< A bool indicating whether the deriv approaches a singularity (sint -> 0)
	bool m_bGrad;										//!< A bool indicating whether or not to compute the gradient
  CQuat m_Quat;										//!< A quat of molecule's orientation
	CQuat m_QuatU;									//!< A quat of molecule's orientation stored temporarily incase of undo
  REAL m_rbu[12];									//!< vector of many real values (sint etc) stored temporarily incase of undo
  Complex m_cbu[3];								//!< vector of exphi, exiphi, exxi and kd stored temporarily incase of undo
  int m_p;												//!< An int of the poles
	int m_pU;												//!< An int of the poles stored temporarily incase of undo
}; // end CRotCoeff class

/////////////////////////////////////
////// Inline functions

//!  CRotCoeff ROT
/*!	Printing/retrieving the n,s,mth index of rot coeff matrix */
inline const Complex &  
CRotCoeff::ROT(int n, int s, int m) const
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_R[n][n+s][m]; 
}
//!  CRotCoeff ROT
/*!	Printing/retrieving the n,s,mth index of rot coeff matrix   */
inline Complex &  
CRotCoeff::ROT(int n, int s, int m) 
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_R[n][n+s][m]; 
}
//!  CRotCoeff dROT
/*!	Printing/retrieving the n,s,mth index of dRot coeff matrix  */
inline const Complex &  
CRotCoeff::dROT(int n, int s, int m) const
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_dR[n][n+s][m]; 
}
//!  CRotCoeff dROT
/*!	Printing/retrieving the n,s,mth index of dRot coeff matrix   */
inline Complex &  
CRotCoeff::dROT(int n, int s, int m) 
{ 
  assert(n < 2*m_p - 1 && n >= 0);
  assert(abs(s) < 2*m_p - 1);
  assert(m >= 0 && m < 2*m_p - 1);
  return m_dR[n][n+s][m]; 
}

//!  CRotCoeff rotate
/*!	Rotate three matrices of tricoeff class
			\param Tin a TriCoeff object for input 
			\param Tout a TriCoeff object for output
			\param p1 an int 
			\param p2 an int 
			\param bFor a boolean of  
*/
inline void 
CRotCoeff::rotate(const CTriCoeff & Tin, CTriCoeff & Tout, int p1, int p2, 
									bool bFor)
{ 
  rotate(Tin[0], Tout[0], p1, p2, bFor);
  rotate(Tin[1], Tout[1], p1, p2, bFor);
  rotate(Tin[2], Tout[2], p1, p2, bFor);
}

//!  CRotCoeff saveUndo
/*!	Save the rotation coefficients and derivatives
		as an undone move
*/
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
} // end saveUndo


//!  CRotCoeff undo
/*!	Revert to old rotation coefficients using saved 
			components
*/
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
} // end undo


#endif
