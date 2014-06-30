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

class CMPE
{
public:
  static CMPE * tp;
  CMPE(const CProtein & mol);
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
  //  static void polarize(vector<CMPE*> & mpe, int p);
  static REAL computeForceAt(const vector<CMPE*> & mpe,
			     const vector<CPnt*> & cen,
			     const CPnt & P, CPnt & force,
			     CPnt & torque, int p);
  static REAL computePotAt(const vector<CMPE*> & mpe, 
			   const vector<CPnt*> & cen,
			   const CPnt & P);
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
  void initialize(const vector<REAL> & charges, const vector<CPnt> & pos,
		  int p);
  REAL recomputeGrad(const vector<CMPE*> & mpe, int i, int j);
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
  static CMCoeff m_tM;
  static CGradCoeff m_tG1, m_tG2, *m_tG;
  static CXForm ** m_xfs;
  static int * IDX;
  static int m_total;

  REAL m_C[N_POLES], m_CD[N_POLES];
  REAL m_rad;
  CMCoeff m_M, m_rM, m_pM, m_L;
  CGradCoeff *  m_pG, m_dL;
  CTorqCoeff m_T, m_rT;
  CRotCoeff m_rot;
  CQuat m_orient, m_orientU;
  int m_id;
  int m_p, m_pU; 
};

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
}

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

inline void
CMPE::undoXForms()
{
  for (int i = 0; i < N_MOL; i++)
    for (int j = i+1; j < N_MOL; j++)
      XFS(i,j).undo();
}


#endif


