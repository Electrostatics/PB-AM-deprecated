#ifndef _PROTEIN_H_
#define _PROTEIN_H_

#include <map>

#include "pdb.h"
#include "util.h"
#include "mpe.h"

using namespace std;

class CProtein
{
public:
  static void initParameters(REAL kappa, REAL dielp, REAL diels, 
			     int nmol, REAL rs);
  static void computeForces(vector<CPnt> & force, vector<CPnt> & torque);
  static bool isInitialized()
    {return m_initialized; }
  static REAL computeDistance(const CProtein & p1, const CProtein & p2)
    { return (p1.getPosition() - p2.getPosition()).norm(); }
  static void getInterfaceAtoms(const CProtein & P1, const CProtein & P2,
				vector <CPnt> & pos);
  
  CProtein(const char * fname);
  void computeTransformTo(const CProtein & target, CQuat & q, 
			  CPnt & trans) const;
 
  bool inCollision(const CProtein & mol);
  void translate(const CPnt & trans)
    { m_undoCenter = m_center; m_center += trans; }
  void setPos(const CPnt & trans)
    { m_center = trans; }
  void rotate(const CQuat & rot)
    {  m_orient = rot * m_orient; m_mpe->setOrient(m_orient); }
  void setOrient(const CQuat & rot)
    { m_orient = rot;  m_mpe->setOrient(m_orient); }
  void PDBoutput(ostream & fout, int sindex, const char * chainid,
		 const CQuat & rot, const CPnt & trans); 
  void QCDoutput(ostream & fout); 
  void untransform()
    { m_center = m_undoCenter; }
  int getNumCharges() const
    { return m_chargedAtoms.size(); }
  REAL getCharge(int i) const
    { return m_chargedAtoms[i]->getCharge(); }
  const CPnt & getChargePos(int i) const
    { return m_chargedAtoms[i]->getPos(); }
  
  REAL getRadius() const
    {return m_rad; }
  REAL getSumCharge() const
    { return m_sumCharge; }
  const CQuat & getOrientation() const
    { return m_orient; }
  const CPnt & getPosition() const
    {return m_center; }
  const vector<CAtom> & getAtoms() const
    { return m_atoms; }
  const CAtom * getAtom(int resnum, int acode) const;
  int getID() const
    { return m_id; }
   
  static map<int, REAL> CHARGES[NUM_AAS];
 private:
  static void loadChargeMap();
  REAL getMaxAtomRad();
  CPnt computeCenter();
  void computeRadius();
  
  static bool m_initialized;
  static vector<CMPE*> m_exps;
  static vector<CPnt*> m_cens;
  static vector<CProtein*> m_mols;
  static bool m_bFirst;
  
  vector<CAtom> m_atoms;
  vector<const CAtom*> m_chargedAtoms;

  CMPE * m_mpe;
  CPnt m_center, m_undoCenter;
  CQuat m_orient, m_undoOrient;
  REAL m_rad;
  REAL m_sumCharge;

  int m_id;
};

#endif
