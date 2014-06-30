#include <iostream>
#include <fstream>

#include "protein.h"
#include "pdb.h"

#define PARAM_FILE "charges_OPLS"
#define MAX_MPE_ERROR 0.02

map<int, REAL> CProtein::CHARGES[NUM_AAS];
bool CProtein::m_initialized = false;
vector<CMPE*> CProtein::m_exps;
vector<CPnt*> CProtein::m_cens;
vector<CProtein*> CProtein::m_mols;
bool CProtein::m_bFirst = true;

void
CProtein::initParameters(REAL kappa, REAL dielp, REAL diels, int nmol, REAL rs)
{
  loadChargeMap();
  CMPE::initConstants(kappa, dielp, diels, nmol, rs);

   m_initialized = true;
}

void
CProtein::loadChargeMap()
{
  ifstream fin;
  fin.open(PARAM_FILE);
   if (!fin.is_open())
    {
      cout << "Could not open input file " << PARAM_FILE << endl;
      exit(0);
    }
  
  char buf[200];
  char rname[4], aname[5];
  float ch;
  while (!fin.eof())
    {
      fin.getline(buf, 199);
      if (buf[0] == '#')
	continue;

      int r = sscanf(buf, "%s %s %g", &rname, &aname, &ch);
      if (r != 3)
	{
	  cout << "Bad input line in file " << PARAM_FILE 
	       << " : " << buf << endl;
	  exit(0);
	}

      CHARGES[AA::getAACode(rname)][CAtom::getAtomCode(aname)] = (REAL)ch;
    }
}

CProtein::CProtein(const char * fname)
{
  vector<AA> aas;
  CPDB::loadFromPDB(fname, aas);

  for (int i = 0; i < aas.size(); i++)
    for (int j = 0; j < aas[i].getNumAtoms(); j++)
      m_atoms.push_back(aas[i][j]);

  for (int i = 0; i < m_atoms.size(); i++)
    if (m_atoms[i].getCharge() != 0.0)
      m_chargedAtoms.push_back(&(m_atoms[i]));

  m_center = computeCenter();
  for (int i = 0; i < m_atoms.size(); i++)
    m_atoms[i].setPos(m_atoms[i].getPos() - m_center);

  computeRadius();

  m_sumCharge = 0.0;
  for (int i = 0; i < getNumCharges(); i++)
    m_sumCharge += getCharge(i);

  REAL s = 0.0;
  for (int i = 0; i < getNumCharges(); i++)
    s += (getCharge(i));
  cout << "sum = " << s <<  " rad = " << m_rad << endl;

  m_mpe = new CMPE(*this);

  m_exps.push_back(m_mpe);
  m_cens.push_back(&m_center);
  m_mols.push_back(this);
}

CPnt
CProtein::computeCenter()
{
  CPnt p;

  for (int i = 0; i < m_atoms.size(); i++)
    p += m_atoms[i].getPos();
	
  p *= 1.0/m_atoms.size();
  return p;
}

void
CProtein::computeRadius()
{
  REAL max = 0.0;
  for (int i = 0; i < m_atoms.size(); i++)
    {
      REAL dist = m_atoms[i].getPos().norm() + m_atoms[i].getRadius();
      if (dist > max)
	max = dist;
    }
  
  m_rad = max;
}


// never used.
void 
CProtein::computeTransformTo(const CProtein & target, CQuat & q, 
			     CPnt & trans) const
{
  trans = conj(getOrientation())*(target.getPosition() - getPosition());
  q =  conj(getOrientation()) * target.getOrientation();
}

bool 
CProtein::inCollision(const CProtein & mol)
{
  if (CProtein::computeDistance(*this, mol) > 
      getRadius() + mol.getRadius())
    return false;
  else
    return true;
}

REAL 
CProtein::getMaxAtomRad()
{
  REAL max = 0.0;
  for (int i = 0; i < m_atoms.size(); i++)
    if (m_atoms[i].getRadius() > max)
      max = m_atoms[i].getRadius();

  return max;
}

void
CProtein::PDBoutput(ostream & fout, int sindex, const char * chainid,
		    const CQuat & rot, const CPnt & trans)
{
  for (int i = 0; i < m_atoms.size(); i++)
    CPDB::writeLine(fout, sindex+i, chainid, m_atoms[i], rot, trans);
}

void
CProtein::QCDoutput(ostream & fout)
{
  char buf[100];
  sprintf(buf, "ATOM\t1\tLTL\tTST\t%8.3f %8.3f %8.3f %5.3f %5.2f",
	  0.0, 0.0, 0.0, 0.0, getRadius());
  fout << buf << endl;

  for (int i = 0; i < m_chargedAtoms.size(); i++)
    {
      const CAtom * a = m_chargedAtoms[i];
      CPnt p = a->getPos();
      //CPnt p = m_orient *  m_atoms[i].getPos() + m_center;
      sprintf(buf, "ATOM\t1\tLTL\tTST\t%8.3f %8.3f %8.3f %5.3f %5.2f",
	      p.x(), p.y(), p.z(),  a->getCharge(), 0.0);
      fout << buf << endl;
    }
}

const CAtom *
CProtein::getAtom(int resnum, int acode) const
{
  for (int i = 0; i < m_atoms.size(); i++)
    if (m_atoms[i].getResNum() == resnum && m_atoms[i].getCode() == acode)
      return &(m_atoms[i]);

  return NULL;
}

void
CProtein::getInterfaceAtoms(const CProtein & P1, const CProtein & P2,
			    vector <CPnt> & pos)
{
  const vector<CAtom> A1 = P1.getAtoms();
  const vector<CAtom> A2 = P2.getAtoms();

  vector<CAtom>::const_iterator it1 = A1.begin();
  vector<CAtom>::const_iterator it2;

  CPnt d = P1.getPosition() - P2.getPosition();
  for (; it1 != A1.end(); it1++)
    {
      if (!it1->isHeavy())
      	continue;

      it2 = A2.begin();
      for (; it2 != A2.end(); it2++)
	{
	  if (!it2->isHeavy())
	    continue;
	  //cout << it1->getPos() << it2->getPos() << endl;
	  if (((it1->getPos() - it2->getPos()) + d).norm() <= 4.0)
	    {
	      pos.push_back(it1->getPos() + P1.getPosition() );
	      break;
	    }
	}
    }

  
  it2 = A2.begin();
  for (; it2 != A2.end(); it2++)
    {
      if (!it2->isHeavy())
	continue;

      it1 = A1.begin();
      for (; it1 != A1.end(); it1++)
	{
	  if (!it1->isHeavy())
	    continue;
	  
	  if ((it2->getPos() - it1->getPos() - d).norm() <= 4.0)
	    {
	      pos.push_back(it2->getPos() + P2.getPosition());
	      break;
	    }
	}
    }
}

void
CProtein::computeForces(vector<CPnt> & force, vector<CPnt> & torque) 
{
  vector<REAL> pot;

  if (m_bFirst)
    {
      CMPE::initXForms(m_exps);
      m_bFirst = false;
    }

  CMPE::updateSolve(m_exps, m_cens);
  CMPE::computeForce(m_exps, m_cens, pot, force, torque); 
}
