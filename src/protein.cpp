#include <iostream>
#include <fstream>

#include "protein.h"
#include "pdb.h"

#define PARAM_FILE "charges_OPLS"				//<! Hard link to the file that contains a list of all OPLS charges for each AA
#define MAX_MPE_ERROR 0.02							//!< Maximum error allowed for the multipole expansions

map<int, REAL> CProtein::CHARGES[NUM_AAS];
bool CProtein::m_initialized = false;
vector<CMPE*> CProtein::m_exps;
vector<CPnt*> CProtein::m_cens;
vector<CProtein*> CProtein::m_mols;
bool CProtein::m_bFirst = true;				

/******************************************************************/
/******************************************************************//**
* Initializing parameters for protein file
* Inputs:  kappa, dielectric of the protein, water 
*				dielectric, the nummber of molecules, average scaling factor
*				RS~avg molecular radius
******************************************************************/
void
CProtein::initParameters(REAL kappa, REAL dielp, REAL diels, int nmol, REAL rs)
{
	// Loading charge map for barnase and barstar
  loadChargeMap();
	// Initialize constants for the multipole expansions
  CMPE::initConstants(kappa, dielp, diels, nmol, rs);

	// Now the system is initialized
  m_initialized = true;
}

/******************************************************************/
/******************************************************************//**
* Function for loading the charge map
******************************************************************/
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

			// Load charges for each amino acid type
      CHARGES[AA::getAACode(rname)][CAtom::getAtomCode(aname)] = (REAL)ch;
    }
}

/******************************************************************/
/******************************************************************//**
* Initializing protein class
******************************************************************/
CProtein::CProtein(const char * fname)
{
  vector<AA> aas;
	string ss = fname;

	if (ss.find("pdb") != std::string::npos)
	{	
		CPDB::loadFromPDB(fname, aas);
	} else if (ss.find("pqr") != std::string::npos)
	{
		CPDB::loadFromPQR(fname, aas);
	}
		
	for (int i = 0; i < aas.size(); i++)							// Add each atom in each amino acid of 
		for (int j = 0; j < aas[i].getNumAtoms(); j++)	// protein to matrix of atoms
			m_atoms.push_back(aas[i][j]);
	
	for (int i = 0; i < m_atoms.size(); i++)				// Add charge on atom in each amino acid to vector of charges
	{
		if (m_atoms[i].getCharge() != 0.0)
				m_chargedAtoms.push_back(&(m_atoms[i]));
	}
	
  m_center = computeCenter();											// compute center of geometry if the atom
  for (int i = 0; i < m_atoms.size(); i++)				// reposition each atom WRT to the center of mass of protein
    m_atoms[i].setPos(m_atoms[i].getPos() - m_center);
	
  computeRadius();
	
  m_sumCharge = 0.0;															// Total charge of protein
  for (int i = 0; i < getNumCharges(); i++)			
    m_sumCharge += getCharge(i);
	
  REAL s = 0.0;																		// Total charge of protein to print as output
  for (int i = 0; i < getNumCharges(); i++)
    s += (getCharge(i));
  cout << "sum = " << s <<  " rad = " << m_rad << endl;
	
  m_mpe = new CMPE(*this);														// Creating a multipole expansion class for the protein
	
  m_exps.push_back(m_mpe);														// Pushing the MPE to an array
  m_cens.push_back(&m_center);												// Adding the protein center to an array of centers
  m_mols.push_back(this);															// Adding the protein to an array of proteins
}

/******************************************************************/
/******************************************************************//**
* Compute center of geometry of the protein of interest
******************************************************************/
CPnt
CProtein::computeCenter()
{
  CPnt p;
	
  for (int i = 0; i < m_atoms.size(); i++)
    p += m_atoms[i].getPos();
	
  p *= 1.0/m_atoms.size();
  return p;
}

/******************************************************************/
/******************************************************************//**
* Calculate radius of protein
******************************************************************/
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

/******************************************************************/
/******************************************************************//**
* Function to determine whether two CG molecules are collided or not
******************************************************************/
bool 
CProtein::inCollision(const CProtein & mol)
{
  if (CProtein::computeDistance(*this, mol) > 
      getRadius() + mol.getRadius())
    return false;
  else
    return true;
}

/******************************************************************/
/******************************************************************//**
* Determine which of all of the atoms in the CG simulation has the
			largest radius.
******************************************************************/
REAL 
CProtein::getMaxAtomRad()
{
  REAL max = 0.0;
  for (int i = 0; i < m_atoms.size(); i++)
    if (m_atoms[i].getRadius() > max)
      max = m_atoms[i].getRadius();

  return max;
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/
/*void
CProtein::PDBoutput(ostream & fout, int sindex, const char * chainid,
		    const CQuat & rot, const CPnt & trans)
{
  for (int i = 0; i < m_atoms.size(); i++)
    CPDB::writeLine(fout, sindex+i, chainid, m_atoms[i], rot, trans);
}*/

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/
/*void
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
}*/


/******************************************************************/
/******************************************************************//**
* Return atom of residue number resnum and atom code acode
******************************************************************/
const CAtom *
CProtein::getAtom(int resnum, int acode) const
{
  for (int i = 0; i < m_atoms.size(); i++)
    if (m_atoms[i].getResNum() == resnum && m_atoms[i].getCode() == acode)
      return &(m_atoms[i]);

  return NULL;
}

/******************************************************************/
/******************************************************************//**
* Function that returns a list of atoms at the interface
******************************************************************/
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

/******************************************************************/
/******************************************************************//**
* Compute forces between molecules.
* Input: vector of forces and torques acting on each molecule in the 
*					system
******************************************************************/
void
CProtein::computeForces(vector<CPnt> & force, vector<CPnt> & torque) 
{
  vector<REAL> pot;
	
  if (m_bFirst)
	{
		CMPE::initXForms(m_exps);						// Initialize transforms if this is the first 
		m_bFirst = false;										// force calc of the simulation
	}
	
  CMPE::updateSolve(m_exps, m_cens);			// Update the solution to mutual polarization.
  CMPE::computeForce(m_exps, m_cens, pot, force, torque); 
}
