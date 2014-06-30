#ifndef _PDB_H_
#define _PDB_H_

#include <vector>

#include "util.h"

#define NUM_ELEM 5
#define NUM_PLACES 10
#define NUM_AAS 22

using namespace std;

class AA;
class CAtom {
 public:
  CAtom(int acode_, int rcode_, int resnum_, REAL x, REAL y, REAL z);
  CAtom(const CAtom & a) :
    m_acode(a.m_acode), m_rcode(a.m_rcode), m_resnum(a.m_resnum),
    m_pos(a.m_pos), m_charge(a.m_charge), m_rad(a.m_rad) {}
  const CAtom & operator=(const CAtom & a)
    { 
      m_acode = a.m_acode; m_rcode = a.m_rcode; m_resnum = a.m_resnum;
      m_pos = a.m_pos; m_charge = a.m_charge; m_rad = a.m_rad;
    }

  static const char ELEMSYM[NUM_ELEM];
  enum ELEMENT {H, C, N, O, S};

  static const char PLACESYM[NUM_PLACES];
  enum PLACE {ALPHA = 1, BETA = 2, GAMMA = 3, DELTA = 4, EPSILON = 5, ZETA = 6, 
	      ETA = 7, NITRO = 8, TERM = 9};
  enum BRANCH {BR1 = 1, BR2 = 2, BR3 = 3};

  static const int ELEM_MASK, PLACE_MASK, BRANCH1_MASK, BRANCH2_MASK;
  static const int ELEM_SHIFT, PLACE_SHIFT, BRANCH1_SHIFT, BRANCH2_SHIFT;

  static int getAtomCode(const char * aname);
  static void generateName(int acode, char * aname);

  void setPos(const CPnt & p)
    { m_pos = p; }

  int getCode() const
    { return m_acode; }
  int getResCode() const
    { return m_rcode; }
  int getResNum() const
    { return m_resnum; }
  REAL getCharge() const
    { return m_charge; }
  const CPnt & getPos() const
    { return m_pos; }
  REAL x() const
    { return m_pos.x(); }
  REAL y() const
    { return m_pos.y(); }
  REAL z() const
    { return m_pos.z(); }
  REAL getRadius() const
    { return m_rad; }
  ELEMENT getElement() const
    { return (ELEMENT) ((m_acode & ELEM_MASK) >> ELEM_SHIFT); } 
  bool isHeavy() const
    { return getElement() != CAtom::H; }

 private:
  int m_acode;
  int m_rcode;
  int m_resnum;
  CPnt m_pos;
  REAL m_charge;
  REAL m_rad;
};


class CPDB;

class AA {
 public:
  const CAtom * operator[](const char * aname) const;
  const CAtom & operator[](int i) const
    { return m_atoms[i]; }
 
  enum AACODE {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, 
	       MET, PHE, PRO, SER, THR, TRP, TYR, VAL, NTR, CTR};
  static const char AANAME[NUM_AAS][4];
  static const char AALETTER[NUM_AAS][2];
  static AACODE getAACode(const char * aa);
  static const char * getName(AACODE code) 
    { return AANAME[(int) code]; }

  AACODE getType() const
    { return m_type; }
  int getNumAtoms() const
    { return m_atoms.size(); }
  void setType(AACODE type)
    { m_type = type; }
  void insertAtom(CAtom & atom)
    { m_atoms.push_back(atom); }
  void clear()
    { m_atoms.clear(); }

 private:
  vector<CAtom> m_atoms;
  AACODE m_type;
};

class CPDB
{
 public:
  static void loadFromPDB(const char * fname, vector<AA> & aas);
  static void writeToPDB(const char * fname, const vector<AA> & aas);

  static void writeLine(ostream & fout, int index, const char * chainid,
		 const CAtom & atom, const CQuat & rot, const CPnt & trans);
  static CAtom readline(const char * buf);
};



#endif
