#ifndef _MPE_H_
#define _MPE_H_

#include <complex>
#include <cmath>
#include <vector>
#include "util.h"
#include "mcoeff.h"
#include "gradcoeff.h"
#include "torqcoeff.h"
#include "xform.h"
#include "rotcoeff.h"

using namespace std;
class CProtein;

//!  The multipole expansion class
/*!
		The class that contains all details for a multipole expansion.  
*/

class CMPE
{
public:
  static CMPE * tp;
	
	//! CMPE class constructor
	/*! 
			\param mol a pointer to the protein object for which the MPE will be created
			\return an object of the CMPE class
	*/
  CMPE(const CProtein & mol);
	
	//! CMPE class constructor
	/*! 
			\param charges a pointer to a vector of charges
			\param pos a pointer to a vector of positions for the aforementioned charges in xyz coords
			\param rad a floating point number of the radius of the CG sphere containing all charges
			\param id an int with numerical ID of the CG sphere
			\param p an int with the number of poles ??
			\return an object of the CMPE class
	*/
  CMPE(const vector<REAL> & charges, const vector<CPnt> & pos, 
       REAL rad, int id, int p);
	
  void reset(int p, const CQuat & Q = CQuat());
	
  static void initConstants(REAL kappa, REAL dielp, REAL diels, int nmol, 
														REAL rs);
  static void undoXForms();
  static void initXForms(const vector<CMPE*> & mpe);
  static void computeForce(vector<CMPE*> & mpe, 
													 const vector<CPnt*> & cen, vector<REAL> & pot,
													 vector<CPnt> & force, vector<CPnt> & torque);
  static void solve(vector<CMPE*> & mpe, const vector<CPnt*> & cen, bool bPot);
  static void updateSolve(vector<CMPE*> & mpe, const vector<CPnt*> & cen);
  static REAL computeForceAt(const vector<CMPE*> & mpe,
														 const vector<CPnt*> & cen,
														 const CPnt & P, CPnt & force,
														 CPnt & torque, int p);
  static REAL computePotAt(const vector<CMPE*> & mpe, 
													 const vector<CPnt*> & cen,
													 const CPnt & P);
													 
	//! Function for mutual polarization 
	/*!
		Function of MPE class that computes the mutual polarization of
		many bodies in a 
	*/
  static void polarize(vector<CMPE*> & mpe, bool bPot);
  static void updateXForms(const vector<CPnt*> & cen, vector<CMPE*> & mpe);
  static void reexpand(const vector<CMPE*> & mpe);
  static void computePairPot(const vector<CMPE*> & mpe, int i, int j, 
														 REAL & p1, REAL & p2);
  void setOrient(const CQuat & Q);
  const CQuat & getOrient() const
	{ return m_orient; }
  int getOrder() const
	{ return m_p; }
  void setOrder(int p);
	
  REAL getRad()
	{ return m_rad; }
	
  void saveUndo();
  void undo();
  
  static int npol;
  static REAL npol_t;
  static bool m_bInfinite;
  static int m_unit;
  int ngpol;
  REAL ngpol_t;
	
private:
  void initCD();
	
	//!  The MPE initialize function
/*!
		Private function that initializes an MPE from a 
		collection of charges and their positions
		\param charges a vector of floating point charges
		\param pos a vector of cartesian coordinates, one per charge
		\param p an integer that generally represents the number of poles
*/
  void initialize(const vector<REAL> & charges, const vector<CPnt> & pos,
									int p);
  REAL recomputeGrad(const vector<CMPE*> & mpe, int i, int j);
	
		//!  The MPE recompute function
/*!
		Recomputes the MP coefficients at molecule i
		\param mpe a vector of MPE objects
		\param i an integer designating molecule of interest
		\return floating point error between MP coeff of iteration n - iteration (n+1)
*/
  REAL recompute(const vector<CMPE*> & mpe, int i);
  void computeForceOn(vector<CMPE*> & mpe, CPnt & force, CPnt & torque, int i);
  void incRotate();
  void incOrder();
	
  void reexpand(const vector<CMPE*> & mpe, int i);
  void reexpandGrad(const vector<CMPE*> & mpe, int i);
  
  static void prepareDTA(const vector<CMPE*> & mpe, int j);
  static void sphToCartMat(const CPnt & P, CPnt * R);
  static CPnt sphToCart(const CPnt * R, const CPnt & P);
  
  static CXForm & XFS(int i, int j) 
	{ assert(i < j); return *m_xfs[IDX[i]+(j-i)-1]; }
	
  static REAL KAPPA, DIEL_S, DIEL_P;
  static int N_MOL;
  static CGradCoeff ** m_tmpG;
  static CMCoeff * m_tmpM;
  static CMCoeff m_tM;					//!< A temporary coefficient object for the multipole expansion
  static CGradCoeff m_tG1, m_tG2, *m_tG;
  static CXForm ** m_xfs;
  static int * IDX;							//!<
  static int m_total;						//!< Sum of ((number of poles per molecule)*(N_POLE+1)/2)
	
  REAL m_C[N_POLES];						//!< Part of GAMMA vector, dielectric boundary crossing operator, coefficients for eq(19) in paper
	REAL m_CD[N_POLES];						//!< Part of DELTA vector, cavity polarization operator. Coefficients for eq(20) in paper
  REAL m_rad;										//!< Radius of the sphere for this object of CMPE
  CMCoeff m_M;									//!< 									
	CMCoeff m_rM;									//!< 
	CMCoeff m_pM;									//!< A coefficient object for the current multipole expansion
	CMCoeff m_L;									//!< A coefficient object for the local expansion 
  CGradCoeff * m_pG;						//!< A vector of gradient coefficients 
	CGradCoeff * m_dL;
  CTorqCoeff m_T;								//!< A coefficient matrix for torques of system.
	CTorqCoeff m_rT;
  CRotCoeff m_rot;
  CQuat m_orient;
	CQuat m_orientU;
  int m_id;											//!< Numerical ID of given MPE
  int m_p;											//!< Has something to do with the poles ??
	int m_pU; 
}; //end CMPE



//!  The MPE setOrder function
/*!
		Function that sets the order ??  
*/
inline void
CMPE::setOrder(int p)
{
  assert(p <= N_POLES);
  if (m_p < p)  
    while (m_p < p)
      incOrder();
  else if (m_p > p)
    {
      for (int k = m_p; k > p; k--)
	m_rot.decOrder(); 
      m_rT.setOrder(p);
      m_rM.setOrder(p);
      m_pM.setOrder(p);
      
      if (m_pG)
	for (int j = 0; j < N_MOL; j++)
	  m_pG[j].setOrder(p);

      m_p = p;
    }
} // end incOrder

//!  The MPE incOrder function
/*!
		Function that sets the order ??  
*/
inline void
CMPE::incOrder()
{
  m_p++;
  m_rot.incOrder();
  incRotate();
  m_pM.copy_p(m_rM, m_p);

  if (m_pG)
    for (int j = 0; j < N_MOL; j++)
      m_pG[j].setOrder(m_p);
}

//!  The MPE sphToCartMat function
/*!
		Function converts spherical coords ??
*/
inline void
CMPE::sphToCartMat(const CPnt & P, CPnt * R)
{
  CSpPnt d = CartToSph(P);
  REAL sin_t = sin(d.theta());
  REAL cos_t = cos(d.theta());
  REAL sin_p = sin(d.phi());
  REAL cos_p = cos(d.phi());
  REAL ir = 1/d.rho();
  REAL irt = ir/sin_t;

  R[0] = CPnt(sin_t*cos_p, ir*cos_t*cos_p, -irt*sin_p);
  R[1] = CPnt(sin_t*sin_p, ir*cos_t*sin_p, irt*cos_p);
  R[2] = CPnt(cos_t, -ir*sin_t, 0.0);
}

//!  The MPE saveUndo function
/*!
		Function  ??
*/
inline void
CMPE::saveUndo()
{
  m_rM.saveUndo();
  m_pM.saveUndo();
  m_L.saveUndo();
  for (int i = 0; i < N_MOL; i++)
    m_pG[i].saveUndo(); 
  m_dL.saveUndo();
  m_rT.saveUndo();
  m_rot.saveUndo();
  m_orientU = m_orient;
 
  m_pU = m_p;
}

//!  The MPE undo function
/*!
		Function  ??
*/
inline void
CMPE::undo()
{
  m_rM.undo();
  m_pM.undo();
  m_L.undo();
  for (int i = 0; i < N_MOL; i++)
    m_pG[i].undo(); 
  m_dL.undo();
  m_rT.undo();
  m_rot.undo();
  m_orient = m_orientU;
 
  m_p = m_pU; 
}

//!  The MPE undoXForms function
/*!
		Function  ??
*/
inline void
CMPE::undoXForms()
{
  for (int i = 0; i < N_MOL; i++)
    for (int j = i+1; j < N_MOL; j++)
      XFS(i,j).undo();
}


#endif


