#include <sys/time.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>

#include "mpe.h"
#include "protein.h"


#define MAX_POLAR_DEV (1e-2)															//!< Max polarization error desired
#define MAX_POLAR_DEV_SQR (MAX_POLAR_DEV*MAX_POLAR_DEV)		//!< Square of max polarization error
#define MAX_POL_ROUNDS 20																	//!< Maximum desired polarization rounds

#define CLKS_PER_MSEC (0.001*CLOCKS_PER_SEC)

// MPE class members described in the mpe.h file
CMPE * CMPE::tp;
REAL CMPE::DIEL_P = 1.0;
REAL CMPE::DIEL_S = 1.0;
int CMPE::N_MOL = 0;
CMCoeff CMPE::m_tM;
CGradCoeff CMPE::m_tG1;
CGradCoeff CMPE::m_tG2;
CGradCoeff *CMPE::m_tG;
CXForm ** CMPE::m_xfs;
int * CMPE::IDX;		
int CMPE::m_total;
bool CMPE::m_bInfinite = false;		
int CMPE::m_unit = 1;							 

/******************************************************************/
/******************************************************************//**
* Initializing parameters for multipole expansion
* Inputs:  kappa, dielectric of the protein, water 
*				dielectric, the nummber of molecules, average scaling factor
* 
******************************************************************/
void
CMPE::initConstants(REAL kappa, REAL dielp, REAL diels, int nmol, REAL rs)
{
  DIEL_P = dielp;		// Dielectric of protein
  DIEL_S = diels;		// Dielectric of solvent
  N_MOL = nmol;		// Number of molecules 
  
	// Initializing other classes: SHCoeff and Transform
  CSHCoeff::init(rs, kappa);
  CXForm::initConstants();

  IDX = new int[N_MOL];		// Create vector of size(number of molecules)
  IDX[0] = 0;							// List starts with 0. 
	// For each subsequent value, tells how many inter-molecule interactions are given
	//  by the previous molecule's interactions plus it's own
	
	// If NMOL=2 	// IDX[0]=0; IDX[1]= = IDX[0] + 1 own interaction = 1
	// If NMOL=3	// IDX[0]=0, IDX[1] = 2, IDX[2] = IDX[1] + 1 own = 3 
	// If NMOL=4	// IDX[0]=0, IDX[1] = 3, IDX[2] = 3+2 = 5, IDX[3] = 3+2+1 = 6 
  for (int i = N_MOL-1; i > 0; i--)
    IDX[N_MOL-i] = IDX[N_MOL-i-1] + i;

  if (m_bInfinite)
    m_xfs = new CXForm*[m_unit * (N_MOL-1)];
  else
    m_xfs = new CXForm*[IDX[N_MOL-1]];
  m_tG = new CGradCoeff[N_MOL];
}

/******************************************************************/
/******************************************************************//**
* Initialize transforms.  Create a transform for each pair of proteins
*   in the system, using their matrix coefficients.
******************************************************************/
void
CMPE::initXForms(const vector<CMPE*> & mpe)
{
  int i_max = (m_bInfinite ? m_unit : N_MOL);		

  for (int i = 0; i < i_max; i++)
    for (int j = i+1; j < N_MOL; j++)
      m_xfs[IDX[i]+(j-i)-1] = new CXForm(mpe[j]->m_pM, mpe[i]->m_pM); 
}

/******************************************************************/
/******************************************************************//**
* Function that sets the orientation of the MPE object given a quaternion
******************************************************************/
void
CMPE::setOrient(const CQuat & rot)
{
  m_rot.reset(rot, m_p);

  m_rot.rotate(m_M, m_rM, true);
  m_rot.rotate(m_T, m_rT, true);
  m_rT.rotate(rot, m_p);
  m_orient = rot;
}


/******************************************************************/
/******************************************************************//**
* Rotate the top level of the expansions associated with this molecule.
******************************************************************/
void
CMPE::incRotate()
{
  assert(m_p <= N_POLES);
 
  m_rot.rotate(m_M, m_rM, m_p, m_p, true);
  m_rot.rotate(m_T, m_rT, m_p, m_p, true);
  m_rT.incRotate(m_orient);
}

/******************************************************************/
/******************************************************************//**
* Initialize a multipole expansion with a protein object
******************************************************************/
CMPE::CMPE(const CProtein & mol) : m_pM(0, N_POLES), m_rot(false) 
{
  m_rad = mol.getRadius();
  m_id = mol.getID();
  initCD();													// Initializing conformational constants for the given protein
	
  int num = mol.getNumCharges();
  vector<CPnt> pos(num);
  vector<REAL> ch(num);
  for (int i = 0; i < num; i++)
	{
		pos[i] = mol.getChargePos(i);
		ch[i] = mol.getCharge(i);
	}
	
  initialize(ch, pos, 1);
  reset(1);
}

/******************************************************************/
/******************************************************************//**
* Create a multipole expansion from a collection of charges.
******************************************************************/
CMPE::CMPE(const vector<REAL> & charges, const vector<CPnt> & pos, 
	   REAL rad, int id, int p) : 
  m_rad(rad),  m_id(id), m_pM(0, N_POLES), m_rot(false)
{ 
  initCD(); 
  initialize(charges, pos, p); 
  reset(p); 
}

/******************************************************************/
/******************************************************************//**
* Initialize many multipole expansion coefficients
******************************************************************/
void
CMPE::initialize(const vector<REAL> & charges, const vector<CPnt> & pos,
		 int p)
{
  m_M = CMCoeff(charges, pos, N_POLES, m_rad);			// Create matrix expansion coeff object
  m_M *= m_C;																				// Transform MCoeff by the Gamma operator

  m_pG = new CGradCoeff[N_MOL];											// Create gradient coefficients for 0-N_POLES
  m_T = CTorqCoeff(charges, pos, N_POLES, m_rad);		// Create torque coefficients
  m_T *= m_C;																				// Transform torque coeffs by the Gamma operator
}

/******************************************************************/
/******************************************************************//**
* CMPE function to reset both the number of poles and 
		the orientation of the multipole expansion
******************************************************************/
void
CMPE::reset(int p, const CQuat & Q)
{
  if (p > N_POLES)
    p = N_POLES;
  m_p = p;

  setOrient(Q);
  m_pM.copy(m_rM, p);

  if (m_pG)
    for (int j = 0; j < N_MOL; j++)
      m_pG[j].reset(p);
}

/******************************************************************/
/******************************************************************//**
* Initializing the surface charge distribution, using equations
*			(19) and (20) from the 2006 paper.
******************************************************************/
void
CMPE::initCD()
{
  if (m_rad == 0.0)
	{
		for (int n = 0; n < N_POLES; n++)
		{
			m_C[n] = 1.0;
			m_CD[n] = 0.0;
		}
	}
	else
	{
		REAL K[N_POLES+1], I[N_POLES+1];
		CSHCoeff::besselk(K, N_POLES+1, CMCoeff::KAPPA * m_rad);
		CSHCoeff::besseli(I, N_POLES+1, CMCoeff::KAPPA * m_rad);
		
		REAL eps = DIEL_P/DIEL_S - 1;
		REAL ex = exp(CMCoeff::KAPPA * m_rad);
		REAL r = m_rad;
		REAL k = CMCoeff::KAPPA*CMCoeff::KAPPA*m_rad*m_rad;
		
		m_C[0] = ex/K[1];															// Coefficients for eq(19) in paper
		m_CD[0] = (r*ex)/K[1] * k*I[1]/3;							// Coefficients for eq(20) in paper
		
		for (int n = 1; n < N_POLES; n++)
		{
			r *= (m_rad*m_rad*CMCoeff::IRS*CMCoeff::IRS);
			m_C[n] = ex*(2*n+1)/((2*n+1)*K[n+1] + n*K[n]*eps);
			m_CD[n] = (r*ex)/((2*n+1)*K[n+1] + n*K[n]*eps) *
			(k*I[n+1]/(2*n+3) - n*I[n]*eps);
		}
	}
}

/******************************************************************/
/******************************************************************//**
* Function to update transforms.  Inputs are protein centers and 
* multipole expansions for each protein.  
******************************************************************/
void 
CMPE::updateXForms(const vector<CPnt*> & cen, vector<CMPE*> & mpe)
{
  int max[N_MOL];
  for (int i = 0; i < N_MOL; i++)
    max[i] = 0;
	
  int i_max = (m_bInfinite ? m_unit : N_MOL);		// If infinite grid, use unit (1) else, use NMOL
  for (int i = 0; i < i_max; i++)
    for (int j = i+1; j < N_MOL; j++)
		{
			int minp = (mpe[i]->getOrder() > mpe[j]->getOrder() ?			// Find the lowest pole order of
									mpe[j]->getOrder() : mpe[i]->getOrder());			// the two proteins
			
			XFS(i,j).init(*(cen[i]) - *(cen[j]), minp);								// reset the dist from 1 to 2 and 
																																//the min # of poles
			
			while (XFS(i,j).isDec() && XFS(i,j).getOrder() > 1)				// If there is less than the tolerated 
				XFS(i,j).decOrder();																		// error in the mut polarization, reduce 
																																// the number of poles for many factors
			
			while (XFS(i,j).isInc() && XFS(i,j).getOrder() < N_POLES)	// Increase number of poles if needed
			{
				assert(mpe[i]->getOrder() >= XFS(i,j).getOrder());
				if (mpe[i]->getOrder() == XFS(i,j).getOrder())
					mpe[i]->incOrder();
				
				assert(mpe[j]->getOrder() >= XFS(i,j).getOrder());
				if (mpe[j]->getOrder() == XFS(i,j).getOrder())
					mpe[j]->incOrder();
				
				XFS(i,j).incOrder();
			} 
			
			int p = XFS(i,j).getOrder();																// Set number of poles
			assert(p <= mpe[i]->getOrder() && p <= mpe[j]->getOrder());	// make sure it is less than or 
			//equal to the poles set for the MPE of each prot
			
			if (p > max[i])
				max[i] = p;
			if (p > max[j])
				max[j] = p;
		}
	
  m_total = 0;
  if (N_MOL > 1)																						//With more than one molecule
	{
		for (int i = 0; i < i_max; i++)
		{
			mpe[i]->setOrder(max[i]);															// Reset the pole order
			m_total += max[i]*(max[i]+1)/2;												// Create a count of total poles
		}
	}
  else
    m_total = mpe[0]->getOrder()*(mpe[0]->getOrder()+1)/2;
} // end updateXForms

/******************************************************************/
/******************************************************************//**
* Routine for mutual polarization.  First computes the total number
* of molecules in the system, then computes the initial mutual 
* polarization and its error.
******************************************************************/
void
CMPE::polarize(vector<CMPE*> & mpe, bool bPot)
{
  if (N_MOL == 1)
    return;
	
  int i_max = (m_bInfinite ? m_unit : N_MOL);		// If infinite grid, use unit (1) else, use NMOL
  REAL itot = 1.0/m_total;											// 1/(N_POLES^2), weighting factor for mutual polarization
  REAL d[i_max], dev = 0.0;
  timeval t1, t2;
	
  for (int i = 0; i < i_max; i++)							
	{
		d[i] = mpe[i]->recompute(mpe, i);
		dev += d[i];
	}
	
  int ct = i_max;
  int i = 0;
  while (dev*itot > MAX_POLAR_DEV_SQR)					// While the mutual polarization error is greater than tol.
	{
		dev -= d[i];
		d[i] = mpe[i]->recompute(mpe, i);
		dev += d[i];
		
		i = (i+1) % i_max;
		if (ct > MAX_POL_ROUNDS*i_max)
		{
			cout << "Polarization does not converge!!! dev=" 
			<< dev << " " << ct << endl;
			exit(0);
		}
		ct++;
	}
	
  if (bPot)																			// If we only want to compute potential, exit
    return;
	
  for (int j = 0; j < i_max; j++)								// Else, compute the gradient for force, torque computations
	{
		dev = 0.0;
		prepareDTA(mpe, j);
		for (int i = 0; i < i_max; i++)
		{
			d[i] = mpe[i]->recomputeGrad(mpe, i, j);
			dev += d[i];
		}
		
		ct = i_max;
		int i = 0;
		while (dev*((1.0/3.0)*itot) > MAX_POLAR_DEV_SQR)
		{
			dev -= d[i];
			d[i] = mpe[i]->recomputeGrad(mpe, i, j);
			dev += d[i];
			
			i = (i+1) % i_max;
			if (ct > MAX_POL_ROUNDS*i_max)
	    {
	      cout << "Gradient polarization does not converge!!! dev=" 
				<< dev << " (wrt " << j << ")" << endl;
	      exit(0);
	    }
			ct++;
		}
		
	}
}  // end polarize

/******************************************************************/
/******************************************************************//**
* Function to perform part of iteration to solve for A matrix
in the paper. Of EQ 51, L = SUM( T * A )
******************************************************************/
void
CMPE::reexpand(const vector<CMPE*> & mpe, int i)
{
  m_L.reset(m_p);
  
  for (int j = 0; j < i; j++)
	{
		XFS(j,i).xform(mpe[j]->m_pM, m_tM, false);			// performing part of EQ 51
		m_L += m_tM;
	}
  for (int j = i+1; j < N_MOL; j++)
	{
		XFS(i,j).xform(mpe[j]->m_pM, m_tM, true);			// performing part of EQ 51
		m_L += m_tM;
	}
}

/******************************************************************/
/******************************************************************//**
* Recompute the MP coeffs at mol i
******************************************************************/
REAL
CMPE::recompute(const vector<CMPE*> & mpe, int i)
{
  if (N_MOL == 1)	// if there is only 1 molecule in the system,
    return 0.0;		// we do not need to recompute the MPEs

  reexpand(mpe, i);

  m_tM = m_L;						// Sum(T * A) of EQ 51
  m_tM *= mpe[i]->m_CD;	// Delta*Sum(T * A) of EQ 51
  m_tM += m_rM;					// Delta*Sum(T * A) + E of EQ 51

  REAL dev = CMCoeff::computeDev(m_pM, m_tM);  // Compute change of EQ 52
  m_pM = m_tM;
  return dev;
}

/******************************************************************/
/******************************************************************//**
* Function to reexpand the gradient of molecule i
******************************************************************/
void
CMPE::reexpandGrad(const vector<CMPE*> & mpe, int i)
{
  m_dL.reset(getOrder());
	
  for (int k = 0; k < i; k++)
	{
		XFS(k,i).xform(mpe[k]->m_pM, m_tG2, false);
		m_dL -= m_tG2;
	}
  for (int k = i+1; k < N_MOL; k++)
	{
		XFS(i,k).xform(mpe[k]->m_pM, m_tG2, true);
		m_dL += m_tG2;
	}
	
  for (int k = 0; k < i; k++)
	{
		XFS(k,i).xform(mpe[k]->m_pG[i], m_tG2, false);
		m_dL += m_tG2;
	}
  for (int k = i+1; k < N_MOL; k++)
	{
		XFS(i,k).xform(mpe[k]->m_pG[i], m_tG2, true);
		m_dL += m_tG2;
	}
}

/******************************************************************/
/******************************************************************//**
* Recompute the gradient of the MP coeffs at mol i respect to the position
* of mol j. Computes EQ 53 in the Lotan 2006 paper.
******************************************************************/
REAL
CMPE::recomputeGrad(const vector<CMPE*> & mpe, int i, int j)
{
  if (N_MOL == 1)
    return 0.0;
	
  m_tG1 = m_tG[i];
  
  for (int k = 0; k < i; k++)					// Computing T*d(A)
	{	
		XFS(k,i).xform(mpe[k]->m_pG[j], m_tG2, false);  
		m_tG1 += m_tG2;										// Adding it to precomputed dT*A
	}
  for (int k = i+1; k < N_MOL; k++)
	{
		XFS(i,k).xform(mpe[k]->m_pG[j], m_tG2, true);
		m_tG1 += m_tG2;
	}
	
  if (i == j)
    m_dL = m_tG1;
	
  m_tG1 *= m_CD;
  REAL dev = CGradCoeff::computeDev(m_pG[j], m_tG1);
	
  m_pG[j] = m_tG1;
  return dev;
}

/******************************************************************/
/******************************************************************//**
* Precompute the sum of the products of dT(i,j)*A(i)
******************************************************************/
void
CMPE::prepareDTA(const vector<CMPE*> & mpe, int j)
{
  m_tG[j].reset(mpe[j]->getOrder());
	
  int i_max = (m_bInfinite ? m_unit : N_MOL);
  for (int i = 0; i < N_MOL; i++)
	{
		if (i == j)
		{
			for (int k = 0; k < i; k++)
	    {
	      XFS(k,i).xform(mpe[k]->m_pM, m_tG2, false);
	      m_tG[i] -= m_tG2;
	    }
			for (int k = i+1; k < N_MOL; k++)
	    {
	      XFS(i,k).xform(mpe[k]->m_pM, m_tG2, true);
	      m_tG[i] += m_tG2;
	    }
		}
		else if (j < i)
		{
			XFS(j,i).xform(mpe[j]->m_pM, m_tG[i], false);
		}
		else
		{
			XFS(i,j).xform(mpe[j]->m_pM, m_tG[i], true);
			m_tG[i].recip();
		}
	}
}

/******************************************************************/
/******************************************************************//**
* Computes the force on molecule i, contains equations of the Lotan
2006 paper, both EQ 37 and EQ 41
******************************************************************/
void
CMPE::computeForceOn(vector<CMPE*> & mpe, CPnt & force, CPnt & torque, int i)
{
  force = -inprod(m_dL, (m_pM+m_pM)-m_rM); // EQ 37: <dL,A> + <L, dA>
  //force = inprod(m_pM, m_dL) + inprod((*m_pG), m_L); // EQ 37: <dL,A> + <L, dA>
  torque = cross(m_dL, m_rT);							// EQ 41: H X dL
}

/******************************************************************/
/******************************************************************//**
* Function to compute pairwise interaction of two molecules in a system. 
******************************************************************/
void
CMPE::computePairPot(const vector<CMPE*> & mpe, int i, int j, 
		     REAL & p1, REAL & p2)
{
  if (i < j)
	{
		XFS(i,j).xform(mpe[j]->m_pM, m_tM, true);
		p1 = inprod(m_tM,mpe[i]->m_pM);
		XFS(i,j).xform(mpe[i]->m_pM, m_tM, false);
		p2 = inprod(m_tM, mpe[j]->m_pM);
	}
  else if (j < i)
	{
		XFS(j,i).xform(mpe[j]->m_pM, m_tM, false);
		p1 = inprod(m_tM,mpe[i]->m_pM);
		XFS(j,i).xform(mpe[i]->m_pM, m_tM, true);
		p2 = inprod(m_tM, mpe[j]->m_pM);
	}
}

/******************************************************************/
/******************************************************************//**
*  CMPE function to compute forces and torques on a system 
******************************************************************/
void
CMPE::computeForce(vector<CMPE*> & mpe, const vector<CPnt*> & cen,
		   vector<REAL> & pot, vector<CPnt> & force, 
		   vector<CPnt> & torque)
{
  int i_max = (m_bInfinite ? m_unit : N_MOL);
	
  pot.clear(); pot.resize(i_max);
  force.clear(); force.resize(i_max);
  torque.clear(); torque.resize(i_max);
  
  for (int j = 0; j < i_max; j++)
	{
		pot[j] = inprod(mpe[j]->m_L, mpe[j]->m_pM);  // Omega in paper EQ 28
		mpe[j]->computeForceOn(mpe, force[j], torque[j], j);
	}
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/
/*REAL
CMPE::computeForceAt(const vector<CMPE*> & mpe, const vector<CPnt*> & cen,
		     const CPnt & P, CPnt & force, CPnt & torque, int p)
{
  force.zero();
  torque.zero();
  REAL pot = 0.0;
	
  CPnt R[3];
  for (int i = 0; i < N_MOL; i++)
	{
		CSpPnt c = CartToSph(P - (*cen[i]));
		CMCoeff L(1.0, P - (*cen[i]), p, CMCoeff::LOCAL_K);
		CGradCoeff dL(L, c);
		
		CPnt f = inprod(mpe[i]->m_pM, dL);
		
		CMPE::sphToCartMat(P - (*cen[i]), R);
		force += CPnt(f.x()*R[0] + f.y()*R[1] + f.z()*R[2]);
		
		pot += inprod(mpe[i]->m_pM, L);
	}
	
  return pot;
}
*/

/******************************************************************/
/******************************************************************//**
* Function to compute potential on a system at a given cartesian point
******************************************************************/
REAL
CMPE::computePotAt(const vector<CMPE*> & mpe, const vector<CPnt*> & cen,
		   const CPnt & P)
{
  REAL pot = 0.0;
	
  for (int i = 0; i < mpe.size(); i++)
	{
		CMCoeff L(1.0, P - (*cen[i]), mpe[i]->getOrder(), CMCoeff::LOCAL_K);
		pot += inprod(mpe[i]->m_pM, L);
	}
	
  return pot;
}

/******************************************************************/
/******************************************************************//**
* A function to run mutual polarization on the system
******************************************************************/
void 
CMPE::solve(vector<CMPE*> & mpe, const vector<CPnt*> & cen, bool bPot)
{ 
  updateXForms(cen, mpe); 
  polarize(mpe, bPot);
  updateXForms(cen, mpe); 
  polarize(mpe, bPot);
}

/******************************************************************/
/******************************************************************//**
* Function of MPE class that reexpands the input matrix pM by
		the matrix transform, as in EQ 46 of paper.
******************************************************************/
void
CMPE::reexpand(const vector<CMPE*> & mpe)
{
  int i_max = (m_bInfinite ? m_unit : N_MOL);
  for (int i = 0; i < i_max; i++)
	{
		mpe[i]->reexpand(mpe, i);
		mpe[i]->reexpandGrad(mpe, i);
	}
}

/******************************************************************/
/******************************************************************//**
* Update solution to multipole expansion.  Calls Transforms and 
*   polarize schemes.  
******************************************************************/
void
CMPE::updateSolve(vector<CMPE*> & mpe, const vector<CPnt*> & cen)
{
  updateXForms(cen, mpe); 
  polarize(mpe, false);
}
