#ifndef _XFORM_H_
#define _XFORM_H_

#include <vector>
#include "mcoeff.h"
#include "gradcoeff.h"
#include "torqcoeff.h"
#include "shcoeff.h"
#include "rotcoeff.h"
#include "transcoeff.h"
#include "util.h"

#define MAX_REL_ERROR 1e-3							//!< The maximum relative order allowed for mutual polarization
#define MIN_REL_ERROR (MAX_REL_ERROR)		//!< The minimum relative order desired for mutual polarization

//!  The XForm class
/*!
		The class that contains all details for the 
		re-expansion operator (section 4.1 of the Lotan 2006 paper) 
*/
class CXForm
{
public:
  CXForm() : m_pM1(NULL), m_pM2(NULL), m_rot(false) {}
	
//!  The XForm class constructor
/*!
		The constructor creates an XForm class from two
		matrix coefficients M1 and M2
		\param M1 is a matrix coefficient from the first molecule
		\param M2 is a matrix coefficient from the 2nd molecule
*/
  CXForm(const CMCoeff & M1, const CMCoeff & M2) : 
	m_pM1(&M1), m_pM2(&M2), m_rot(true) {}
  void init(const CPnt & P, int p)
	{ reset(P, p); }
  
  static void initConstants();

//!  The XForm class xform function
/*!
		The function takes a Matrix coefficient class and computes
		the gradient of it.
		\param Min is a matrix coefficient for transformation
		\param Gout is a gradient coefficient object made from function
		\param bFor is a boolean to indicate whether the input should be transposed
									for translation coefficients */
  void xform(const CMCoeff & Min, CGradCoeff & Gout, bool bFor);
//!  The XForm class xform function
/*!
		The function takes a Matrix coefficient class and computes
		the re-expansion of it, Z= T dot X
		\param Min is a matrix coefficient for transformation
		\param Mout is a re-expanded coefficient object made from function
		\param bFor is a boolean to indicate whether the input should be transposed
									for translation coefficients */
  void xform(const CMCoeff & Min, CMCoeff & Mout, bool bFor);
//!  The XForm class xform function
/*!
		The function takes a triplet matrix coefficient class and computes
		the re-expansion of the three matrices, Z= T dot X
		\param Gin is a triplet matrix coefficient for transformation
		\param Gout is a re-expanded triplet coefficient object made from function
		\param bFor is a boolean to indicate whether the input should be transposed
									for translation coefficients */	
  void xform(const CTriCoeff & Gin, CTriCoeff & Gout, bool bFor);
	
//!  The XForm class reset function
/*! Reset the XForm object by the given XYZ coords and p number of poles  */	
  void reset(const CPnt & P, int p);

//!  The XForm class sphToCart function
/*! Convert a partial deriv in spherical coordinates to cartesian: 
		[dF/dr, dF/dt, dF,dphi] to [dF/dx, dF/dy, dF,dz]
		\param p a 3 coordinate function F to convert	*/		
  void sphToCart(CPnt & p);
//!  The XForm class sphToCart function
/*! Convert a partial deriv of GradCoeff object in spherical to cartesian: 
		\param Gin a GradCoeff object to convert	*/	
  void sphToCart(CGradCoeff & Gin);

// Return the number of poles	
  int getOrder() const
	{ return m_p; }
  REAL getError() const
	{ return m_relError; }

// increase or decrease number of poles
  REAL incOrder();
  REAL decOrder();
  
  bool isInc()
	{ return (m_relError >= MAX_REL_ERROR); }
  bool isDec()
	{ return (m_relError < MIN_REL_ERROR); }
	
  void saveUndo();
  void undo();
	
private:
  void compRelError()
	{
		REAL a = fabs(m_resid1) + fabs(m_resid2);
		m_relError = (m_p > 2 ? a/fabs(m_base) : a);
	}
  
  CRotCoeff m_rot;				//!< A rotCoeff object used for computing rotation coefficients
  CTransCoeff m_trans;		//!< A TransCoeff object used for computing translation coefficients 
  int m_p;								//!< The number of poles
	int	m_pU;								//!< The number of poles, stored incase of undo
  const CMCoeff * m_pM1;  //!< The matrix coefficient of the first molecule
	const CMCoeff * m_pM2;  //!< The matrix coefficient of the second molecule
  REAL m_relError;				//!< A floating point of the relative error of the two residues over the basis
	REAL m_resid1;					//!< Error in S*pM1 VS S*pM2 for p-1 poles
	REAL m_resid2;					//!< Error in S*pM1 VS S*pM2 for p poles
	REAL m_base;						//!< Inner product of the first matrix R*pM1 and second S*R*pm2
  REAL m_rbu[4];					//!< A vector of stored class coefficients stored for undo case
  CMCoeff m_tM1;					//!< A temporary matrix used for computing rotation/translation transforms
	CMCoeff m_tM2;					//!< A temporary matrix used for computing rotation/translation transforms
	CMCoeff m_tM3;					//!< A temporary matrix used for computing rotation/translation transforms
  CGradCoeff m_tG1;				//!< A temporary matrix used for computing rotation/translation gradient transforms
	CGradCoeff m_tG2;				//!< A temporary matrix used for computing rotation/translation gradient transforms
  CPnt m_R[3];						//!< The SphToCart conversion matrix for dz/dSp at r, t, phi 
	CPnt m_RU[3];						//!< The SphToCart conversion matrix for dz/dSp at r, t, phi stored for undo case 
};

///////////////////////////////////////////
/////// Inline functions

//!  The XForm saveUndo function
/*! 		Function that saves the current step, for cases of undo later */
inline void
CXForm::saveUndo()
{
  m_rot.saveUndo();
  m_trans.saveUndo();

  m_rbu[0] = m_relError;
  m_rbu[1] = m_resid1;
  m_rbu[2] = m_resid2;
  m_rbu[3] = m_base;

  m_RU[0] = m_R[0];
  m_RU[1] = m_R[1];
  m_RU[2] = m_R[2];

  m_pU = m_p;
}

//!  The XForm undo function
/*!
		Function that replaces current XForm values with
		saved values to revert to an old state.
*/
inline void
CXForm::undo()
{
  m_rot.undo();
  m_trans.undo();

  m_relError = m_rbu[0];
  m_resid1 = m_rbu[1];
  m_resid2 = m_rbu[2];
  m_base = m_rbu[3];

  m_R[0] = m_RU[0];
  m_R[1] = m_RU[1];
  m_R[2] = m_RU[2];

  m_p = m_pU;
}

#endif
