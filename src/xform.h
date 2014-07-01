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

/******************************************************************//**
* # File: xform.h
* #
* # Date: June 2014
* #
* # Description: This file is the header file for the class XForm
* #								Xform = transform
* #
* # Author: Lotan, Felberg
* #
* # Copyright ( c )
* #
******************************************************************/


#define MAX_REL_ERROR 1e-3							//!< The maximum relative order allowed for mutual polarization
#define MIN_REL_ERROR (MAX_REL_ERROR)		//!< The minimum relative order desired for mutual polarization

//!  The XForm class
/*!
		The class that contains all details for a transform ??  
*/
class CXForm
{
 public:
  CXForm() : m_pM1(NULL), m_pM2(NULL), m_rot(false) {}
  CXForm(const CMCoeff & M1, const CMCoeff & M2) : 
    m_pM1(&M1), m_pM2(&M2), m_rot(true) {}
  void init(const CPnt & P, int p)
    { reset(P, p); }
  
  static void initConstants();
    
  void xform(const CMCoeff & Min, CGradCoeff & Gout, bool bFor);
  void xform(const CMCoeff & Min, CMCoeff & Mout, bool bFor);
  void xform(const CTriCoeff & Gin, CTriCoeff & Gout, bool bFor);
  void reset(const CPnt & P, int p);
 
  void sphToCart(CPnt & p);
  void sphToCart(CGradCoeff & Gin);

  int getOrder() const
    { return m_p; }
  REAL getError() const
    { return m_relError; }
  
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
  
  CRotCoeff m_rot;
  CTransCoeff m_trans;
  int m_p, m_pU;
  const CMCoeff * m_pM1, *m_pM2;
  REAL m_relError, m_resid1, m_resid2, m_base;
  REAL m_rbu[4];
  CMCoeff m_tM1, m_tM2, m_tM3;
  CGradCoeff m_tG1, m_tG2;
  CPnt m_R[3], m_RU[3];
};

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

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

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

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
