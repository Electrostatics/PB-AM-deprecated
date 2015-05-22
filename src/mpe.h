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
/*! The class that contains all details for a multipole expansion.  */

class CMPE
{
public:
  static CMPE * tp;
	
	//! CMPE class constructor
	/*! \param mol a pointer to the protein object for which the MPE will be created
			\return an object of the CMPE class */
  CMPE(const CProtein & mol);
	
	//! CMPE class constructor
	/*! Constructs an MPE object for the computation of self-polarization
			\param charges a pointer to a vector of charges
			\param pos a pointer to a vector of positions for the aforementioned charges in xyz coords
			\param rad a floating point number of the radius of the CG sphere containing all charges
			\param id an int with numerical ID of the CG sphere
			\param p an int with the number of poles ??
			\return an object of the CMPE class 	*/
  CMPE(const vector<REAL> & charges, const vector<CPnt> & pos, REAL rad, int id, int p);

	//! CMPE reset function
	/*!  CMPE function to reset both the number of poles and 
		the orientation of the multipole expansion
		\param p is an integer number of poles
		\param Q is a CQuat object for rotating the MPE  	*/
  void reset(int p, const CQuat & Q = CQuat());
	
		//! CMPE initConstants function
	/*!  CMPE function to initialize constants of MPE functions
		\param kappa a floating point of the inverse debye length
		\param dielp a floating point of the dielectric of the protein
		\param diels a floating point of the dielectric of the solvent
		\param nmol an integer of the number of molecules in the system
		\param rs a floating point of a scaling factor for length 	*/
  static void initConstants(REAL kappa, REAL dielp, REAL diels, int nmol, REAL rs);

		//! CMPE initXForms function
	/*!  CMPE function to initialize the transforms of the MPE functions
		\param mpe a vector of MPE objects */	
  static void initXForms(const vector<CMPE*> & mpe);

		//! CMPE computeForce function
	/*!  CMPE function to compute forces and torques on a system
		\param mpe a vector of MPE objects 
		\param cen a vector of CG object centers in XYZ coordinates
		\param pot a vector of potentials, one for each molecule 
		\param force a vector of forces, one for each molecule 
		\param torque a vector of torques, one for each molecule */	
  static void computeForce(vector<CMPE*> & mpe, const vector<CPnt*> & cen, vector<REAL> & pot,
													 vector<CPnt> & force, vector<CPnt> & torque);
											
 // static REAL computeForceAt(const vector<CMPE*> & mpe, const vector<CPnt*> & cen,
//														 const CPnt & P, CPnt & force, CPnt & torque, int p);

		//! CMPE computePot function
	/*!  CMPE function to compute potential on a system at a given cartesian point
		\param mpe a vector of MPE objects 
		\param cen a vector of CG object centers in XYZ coordinates
		\param P a vector of cartesian coordinates to compute the potential at
		\return a floating point of the potential at point P  */	
  static REAL computePotAt(const vector<CMPE*> & mpe, const vector<CPnt*> & cen, const CPnt & P);
	
		//! CMPE solve function
	/*!  A function to run mutual polarization on the system
		\param mpe a vector of MPE objects    */		
  static void solve(vector<CMPE*> & mpe, const vector<CPnt*> & cen, bool bPot);

	//! Function for reexpansion of MPE 
	/*! Function of MPE class that reexpands the input matrix pM by
		the matrix transform, as in EQ 46 of paper.
	\param mpe a vector of MPE objects to use for polarization	*/
  static void reexpand(const vector<CMPE*> & mpe);

	//! Function for updating the solutions of the mutual MPE 
	/*! Update solution to multipole expansion.  Calls Transforms and 
				polarize schemes.
	\param mpe a vector of MPE objects to use for polarization	
	\param cen a vector of centers of CG spheres */	
	static void updateSolve(vector<CMPE*> & mpe, const vector<CPnt*> & cen);
													 
	//! Function for mutual polarization 
	/*! Function of MPE class that computes the mutual polarization of
	many bodies in a system, given a vector of the MPE objects.
	\param mpe a vector of MPE objects to use for polarization
	\param bPot a boolean of whether or not to only compute potential
					true for only potential, false for gradients as well 	*/
  static void polarize(vector<CMPE*> & mpe, bool bPot);

	//! Function for updateXForms 
	/*! Function of MPE class that updates transforms. 
	\param cen a vector of centers of CG spheres
	\param mpe a vector of MPE objects to use for polarization 	*/		
  static void updateXForms(const vector<CPnt*> & cen, vector<CMPE*> & mpe);

		//!  The MPE computePairPot function
	/*! Function to compute pairwise interaction of two molecules in a system.
		\param mpe a vector of MPEs, one for each molecule in the system
		\param i an integer of the index of the first molecule of interest
		\param j an integer of the index of the second molecule of interest
		\param p1 a floating point of the potential on molecule 1
		\param p2 a floating point of the potential on molecule 2 */
  static void computePairPot(const vector<CMPE*> & mpe, int i, int j, REAL & p1, REAL & p2);

		//!  The MPE setOrient function
	/*! Function that sets the orientation of the MPE object given a quaternion
		\param Q a quaternion object to rotate the MPE */
  void setOrient(const CQuat & Q);

// Functions inline
	
		//!  The MPE getOrient function
	/*! Function that returns the orientation of the MPE object */
  const CQuat & getOrient() const
	{ return m_orient; }
		//!  The MPE getOrder function
	/*! Function that returns number of poles of the MPE object */
  int getOrder() const
	{ return m_p; }
		//!  The MPE getRad function
	/*! Function that returns the radius of the MPE object */	
  REAL getRad()
	{ return m_rad; }

// Inline functions, below	
  void setOrder(int p);
  void saveUndo();
  void undo();
  static void undoXForms();

  static bool m_bInfinite;	//!< Indicates whether simulation for an infinite grid is being performed
  static int m_unit;				//!< Unit needed for an infinite grid.  Set to 1.
  int ngpol;
  REAL ngpol_t;
	
private:

		//!  The MPE initCD function
	/*! 		Computes Gamma and Delta coefficients for the MPE, using equations
			(19) and (20) from the 2006 paper. */
  void initCD();
	
	//!  The MPE initialize function
	/*! Private function that initializes an MPE from a 
		collection of charges and their positions
		\param charges a vector of floating point charges
		\param pos a vector of cartesian coordinates, one per charge
		\param p an integer that generally represents the number of poles */
  void initialize(const vector<REAL> & charges, const vector<CPnt> & pos, int p);

		//!  The MPE recomputeGrad function
	/*! Function that computes the gradient of A matrix through iterative computation
		\param mpe a vector of MPE objects to compute gradients
		\param i an integer of the ith molecule 
		\param j an integer of the jth molecule  */
  REAL recomputeGrad(const vector<CMPE*> & mpe, int i, int j);
	
		//!  The MPE recompute function
	/*! Recomputes the MP coefficients at molecule i
		\param mpe a vector of MPE objects
		\param i an integer designating molecule of interest
		\return floating point error between MP coeff of iteration n - iteration (n+1) */
  REAL recompute(const vector<CMPE*> & mpe, int i);

		//!  The MPE reexpandGrad function
	/*! Recomputes the gradient MP coefficients at molecule i, in accordance with EQ 53
		\param mpe a vector of MPE objects
		\param i an integer designating molecule of interest */
  void reexpandGrad(const vector<CMPE*> & mpe, int i);	
	
		//!  The MPE reexpand function
	/*! Recomputes the MP coefficients at molecule i, in accordance with EQ 51
		\param mpe a vector of MPE objects
		\param i an integer designating molecule of interest */
  void reexpand(const vector<CMPE*> & mpe, int i);

		//!  The MPE computeForceOn function
	/*! Computes the force on molecule i
			\param mpe a vector of MPEs, one for each molecule in the system
			\param force a xyz vector of forces
			\param torque a xyz vector of torques
			\param i an integer of the molecule index   */
  void computeForceOn(vector<CMPE*> & mpe, CPnt & force, CPnt & torque, int i);

		//!  The MPE incRotate function
	/*! Increases the number of poles in rotation coefficients */ 
  void incRotate();

		//!  The MPE recompute function
	/*! Precomputes the sum of the products of dT(i,j)*A(j) at molecule j
		\param mpe a vector of MPE objects
		\param i an integer designating molecule of interest*/  
  static void prepareDTA(const vector<CMPE*> & mpe, int j);

	// Converting spherical to cartesian coordinates
  static CPnt sphToCart(const CPnt * R, const CPnt & P);
  static void sphToCartMat(const CPnt & P, CPnt * R);

		//!  The MPE XFS function
	/*! 		Returns a pointer to the index for the transform
		between molecule i and j. */  
  static CXForm & XFS(int i, int j) 
	{ assert(i < j); return *m_xfs[IDX[i]+(j-i)-1]; }

// Inline functions
  void incOrder();
	
  static REAL KAPPA; 						//!< A floating point number of the inverse debye length
	static REAL DIEL_S;						//!< A floating point number of the dielectric of the solvent
	static REAL DIEL_P;						//!< A floating point number of the dielectric of the protein
  static int N_MOL;							//!< An integer of the number of molecules in the system
	
// Gradient matrix coefficients
  static CGradCoeff ** m_tmpG;
  static CMCoeff * m_tmpM;

  static CGradCoeff m_tG1;			//!< A temporary gradient coefficient, equal to sum( delT*A ) in EQ 53
	static CGradCoeff m_tG2;			//!< A temporary gradient coefficient, equal to ( delT*A ) in EQ 53 
  CGradCoeff * m_pG;						//!< A vector of gradient coefficients, equivalent do dA 
	static CGradCoeff *m_tG;			//!< A vector of transformed gradient transforms, of length NMOL
	
  static CXForm ** m_xfs;				//!< An array of transforms that holds one for each pair of mols in the system
  static int * IDX;							/** An integer array of summed interactions.  For NMOL=4, 
																		 IDX[0]=0, IDX[1] = 3, IDX[2] = 3+2 = 5, IDX[3] = 3+2+1 = 6  */
  static int m_total;						//!< Sum of ((number of poles per molecule)*(N_POLE+1)/2)
	
	// Form factors, gamma and delta 
  REAL m_C[N_POLES];						/** Part of GAMMA vector, dielectric boundary crossing operator, 
																										coefficients for eq(19) in paper */
	REAL m_CD[N_POLES];						/** Part of DELTA vector, cavity polarization operator. 
																										Coefficients for eq(20) in paper */
																										
  REAL m_rad;										//!< Radius of the sphere for this object of CMPE

// Matrix coefficients for EQ 51, used for solving for A matrix of system
  static CMCoeff m_tM;					//!< A transform coefficient object for the multipole expansion, L in paper
  CMCoeff m_M;									//!< A matrix coefficient object transformed by Gamma 									
	CMCoeff m_rM;									//!< Coefficients multipole expansion, E in paper
	CMCoeff m_pM;									//!< A gradient coefficient object for the current multipole expansion, A in paper
	CMCoeff m_L;									//!< A coefficient object for the local expansion 

	CGradCoeff m_dL;							//!< A gradient of the local expansion 
  CTorqCoeff m_T;								//!< A coefficient matrix for torques of system.
	CTorqCoeff m_rT;							//!< A rotated coefficient matrix for torques of system, the H matrix of EQ 41
  CRotCoeff m_rot;							//!< A rotation coefficient object for MPE
  CQuat m_orient;								//!< A quaternion describing the orientation of the MPE
	CQuat m_orientU;;							//!< A quaternion describing the orientation of the MPE, saved incase of undo
  int m_id;											//!< Numerical ID of given MPE
  int m_p;											//!< Number of poles for this MPE object
	int m_pU; 										//!< Number of poles for this MPE object, saved incase of undo
}; //end CMPE

///////////////////////////////////////////
/////// Inline functions

//!  The MPE setOrder function
/*! Sets the order (number of poles) for the MPE object
	\param p an int describing the number of poles desired.
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
/*! Function that increases the order (number of poles)
	used for setting the number of poles in the MPE object
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
/*! Function that converts cartesian coordinates given
	as input to a saved matrix used for converting spherical
	derivatives to cartesian derivatives, as in EQ 48.
	\param P a cartesian coordinate object for deriv. evaluation
	\param R a vector of cartesian coordinates used for conversion
					in EQ 48.  */
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
/*!  Function that stores everthing in the MPE object 
		to U parameters incase of future undo */

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
/*!  Function that reverts MPE object to saved parameters
	using everything stored in the objects within the class */
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
/*! That reverts all transforms between each pair of
molecules in the system.  */
inline void
CMPE::undoXForms()
{
  for (int i = 0; i < N_MOL; i++)
    for (int j = i+1; j < N_MOL; j++)
      XFS(i,j).undo();
}


#endif


