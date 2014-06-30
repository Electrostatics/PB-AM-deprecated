#include <cmath>
#include "xform.h"
#include "gradcoeff.h"

void
CXForm::initConstants()
{
  CRotCoeff::initConstants();
  CTransCoeff::initConstants();
}

// Initialize to a MtoL transform by the vector P.
void 
CXForm::reset(const CPnt & P, int p)
{
  assert(p > 0 && p <= N_POLES);

  m_p = p;
  CSpPnt S = CartToSph(P);
  m_rot.reset(S.theta(), S.phi(), M_PI, p);
  m_trans.reset(S.rho(), p);

  m_rot.rotate(*m_pM1, m_tM1, 1, p, true);
  m_rot.rotate(*m_pM2, m_tM2, 1, p, true);
  m_trans.translate(m_tM1, m_tM3, m_p, false);
  m_base = inprod(m_tM2, m_tM3);	  
  
  m_resid1 = (p > 1 ? m_trans.computeError(m_tM1, m_tM2, p-1) : 0.0);
  m_resid2 = m_trans.computeError(m_tM1, m_tM2, p);
  compRelError();

  REAL sint, cost, sinp ,cosp;
  m_rot.getParams(sint, cost, sinp, cosp);
  REAL ir = 1.0/S.rho();

  if (m_rot.isSingular())
    {
      m_R[0] = CPnt(0.0, ir*cost, 0.0);
      m_R[1] = CPnt(0.0, 0.0, ir);
      m_R[2] = CPnt(cost, 0.0, 0.0);
    }
  else 
    {
      REAL irst = ir/sint, irct = ir*cost;
            
      m_R[0] = CPnt(sint*cosp, irct*cosp, -irst*sinp);
      m_R[1] = CPnt(sint*sinp, irct*sinp, irst*cosp);
      m_R[2] = CPnt(cost, -ir*sint, 0.0);
    }
}

// Compute the derivatives of the transformed MP coeffs
void 
CXForm::xform(const CMCoeff & Min, CGradCoeff & Gout, bool bFor)
{
  assert(Min.getOrder() >= m_p);

  // Derivative with respect to rho
  m_rot.rotate(Min, m_tM1, true);
  m_trans.dTranslate(m_tM1, m_tM2, !bFor); 
  m_rot.rotate(m_tM2, Gout[dRHO], false);
  
  // Derivative with respect to theta and phi (first part)
  m_trans.translate(m_tM1, m_tM2, !bFor);
  m_rot.dRotateT(m_tM2, Gout[dTHETA], false);
  m_rot.dRotateP(m_tM2, Gout[dPHI], false);
  
  // Derivative with respect to theta (second part)
  m_rot.dRotateT(Min, m_tM1, true);
  m_trans.translate(m_tM1, m_tM2, !bFor);
  m_rot.rotate(m_tM2, m_tM1, false);
  Gout[dTHETA] += m_tM1;
  
  // Derivative with respect to phi (second part)
  m_rot.dRotateP(Min, m_tM1, true);
  m_trans.translate(m_tM1, m_tM2, !bFor);
  m_rot.rotate(m_tM2, m_tM1, false);
  Gout[dPHI] += m_tM1;

  sphToCart(Gout);
}

// Transform the MP coeffs
void 
CXForm::xform(const CMCoeff & Min, CMCoeff & Mout, bool bFor)
{
  assert(Min.getOrder() >= m_p);

  m_rot.rotate(Min, m_tM1, true);
  m_trans.translate(m_tM1, m_tM2, !bFor); 
  m_rot.rotate(m_tM2, Mout, false);
}

// Transform the MP coeff triplet.
void 
CXForm::xform(const CTriCoeff & Gin, CTriCoeff & Gout, bool bFor)
{
  xform(Gin[0], Gout[0], bFor);
  xform(Gin[1], Gout[1], bFor);
  xform(Gin[2], Gout[2], bFor);
}

REAL
CXForm::incOrder()
{
  assert(m_p < N_POLES);

  m_p++;
  m_trans.incOrder();
  m_rot.incOrder();

  if (m_tM1.getOrder() < m_p)
    {
      m_rot.rotate(*m_pM1, m_tM1, m_p, m_p, true);
      m_rot.rotate(*m_pM2, m_tM2, m_p, m_p, true);
    }
 
  m_resid1 = m_resid2;
  m_resid2 = m_trans.computeError(m_tM1, m_tM2, m_p);
  m_base += m_resid2;
  compRelError();
     
  return m_relError;
}

REAL
CXForm::decOrder()
{
  assert (m_p > 1);
    
  m_p--;

  m_trans.decOrder();
  m_rot.decOrder();

  m_base -= m_resid2;
  m_resid2 = m_resid1;
  m_resid1 = (m_p > 1 ? m_trans.computeError(m_tM1, m_tM2, m_p-1) : 0.0);
  compRelError();
  
  return m_relError;
}

void
CXForm::sphToCart(CPnt & p)
{
  p = CPnt(dot(m_R[0],p), dot(m_R[1],p), dot(m_R[2],p));
}

void
CXForm::sphToCart(CGradCoeff & G)
{
  for (int n = 0; n < G.getOrder(); n++)
    for (int m = 0; m <= n; m++)
      {
	CPnt ar(G[dRHO](n,m).real(), G[dTHETA](n,m).real(), 
		G[dPHI](n,m).real());
	CPnt ai(G[dRHO](n,m).imag(), G[dTHETA](n,m).imag(), 
		G[dPHI](n,m).imag());

	CPnt br(dot(m_R[0],ar), dot(m_R[1],ar), dot(m_R[2],ar));
	CPnt bi(dot(m_R[0],ai), dot(m_R[1],ai), dot(m_R[2],ai));

	G[dX](n,m) = Complex(br.x(), bi.x());
	G[dY](n,m) = Complex(br.y(), bi.y());
	G[dZ](n,m) = Complex(br.z(), bi.z());
      }
}
