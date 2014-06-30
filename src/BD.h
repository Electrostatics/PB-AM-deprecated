#ifndef _BD_H_
#define _BD_H_

#include <vector>

#include "util.h"
#include "pdb.h"
#include "protein.h"

using namespace std;

struct CONTACT
{
  CONTACT(const CAtom * a1_, const CAtom * a2_, REAL dist_) :
    a1(a1_), a2(a2_), dist(dist_) {}
  const CAtom * a1, * a2;
  REAL dist;
};


class CBD
{
 public:
  CBD(const char * fname1, const char * fname2, REAL kappa);
  enum STATUS {ESCAPED, DOCKED, RUNNING};
  static const char STATUS_STRING[3][10];
  CBD::STATUS run();
  REAL computeRate(int numTraj, int nDocked);
  static CPnt getRandVec(REAL std)
    { return std*CPnt(normRand(), normRand(), normRand()); }

 private:
  void saveState();
  bool isDocked() const;
  bool escaped(REAL dist) const;
  void loadContacts(const char * fconts);
  REAL compute_dt() const; 
  CBD::STATUS makeMove(const CPnt & dR2, const CPnt & dO1, 
		const CPnt & dO2, REAL dt);
  void computeInterfaceVectors(const char * fname1, const char * fname2);
  
  
  CProtein * m_mol1, * m_mol2;
  REAL m_Dtr;
  REAL m_Dr1, m_Dr2;
  vector<CONTACT> m_contacts;
  REAL m_maxContDist;
  int m_fcnt;
  CPnt m_patch1, m_patch2, m_orth1, m_orth2;
  CQuat m_rot1, m_rot2;
};

inline bool
CBD::escaped(REAL dist) const
{
  return (CProtein::computeDistance(*m_mol1, *m_mol2) > dist);
}

inline REAL 
CBD::compute_dt() const
{
  REAL d = CProtein::computeDistance(*m_mol1, *m_mol2) - 
    m_maxContDist;

  if (d-5 > 0)
    return 1.0 + (d-5)/15;
  else
    return 1.0;
}


#endif
