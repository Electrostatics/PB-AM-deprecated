#include <fstream>
#include <iostream>

#include "pdb.h"
#include "protein.h"

/******************************************************************//**
 * #
 * # File: pdb.cpp
 * #
 * # Date: June 2014
 * #
 * # Description: This file details the class Atom and its functions
 * #
 * # Author: Lotan, Felberg
 * #
 * # Copyright ( c )
 * #
******************************************************************/

const char CAtom::ELEMSYM[NUM_ELEM] = {'H', 'C', 'N', 'O', 'S'};					//!< Atomic symbols
const char CAtom::PLACESYM[NUM_PLACES] = {' ', 'A', 'B', 'G', 'D', 'E', 
					     'Z', 'H', 'N', 'T'};																				//!< Letter codes for location of atoms

const int CAtom::ELEM_SHIFT = 8;
const int CAtom::PLACE_SHIFT = 4;
const int CAtom::BRANCH1_SHIFT = 2;
const int CAtom::BRANCH2_SHIFT = 0;

const int CAtom::ELEM_MASK = (7 << CAtom::ELEM_SHIFT);
const int CAtom::PLACE_MASK = (15 <<CAtom::PLACE_SHIFT); 
const int CAtom::BRANCH1_MASK = (3 << CAtom::BRANCH1_SHIFT);
const int CAtom::BRANCH2_MASK = (3 << CAtom::BRANCH2_SHIFT);

const char AA::AANAME[NUM_AAS][4] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", 
				 "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", 
				 "MET", "PHE", "PRO", "SER", "THR", "TRP", 
				 "TYR", "VAL", "NTR", "CTR"};																//!< 3 Letter codes for each amino acid
const char AA::AALETTER[NUM_AAS][2] = {"A", "R", "N", "D", "C", "Q", "E", "G", 
				   "H", "I", "L", "K", "M", "F", "P", "S", 
				   "T", "W", "Y", "V", "X", "Z"};														//!< 1 Letter codes for each amino acid
					 
/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

CAtom::CAtom(int acode_, int rcode_, int resnum_, REAL x, REAL y, REAL z) :
  m_acode(acode_), m_rcode(rcode_), m_resnum(resnum_), m_pos(x, y, z), 
  m_charge(0.0), m_rad(0.0)
{
  map<int,REAL>::const_iterator it = CProtein::CHARGES[m_rcode].find(m_acode);
  if (it != CProtein::CHARGES[m_rcode].end() && it->second != 0.0)
    m_charge = CProtein::CProtein::CHARGES[m_rcode][m_acode];

  m_rad = 2.0;
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

void 
CAtom::generateName(int acode, char * aname)
{
  int branch2 = (acode & BRANCH2_MASK) >> BRANCH2_SHIFT;
  if (branch2)
    aname[0] = '0' + branch2;
  else 
    aname[0] = ' ';

  aname[1] = ELEMSYM[(acode & ELEM_MASK) >> ELEM_SHIFT];
  aname[2] = PLACESYM[(acode & PLACE_MASK) >> PLACE_SHIFT];

  int branch1 = (acode & BRANCH1_MASK) >> BRANCH1_SHIFT;
  if (branch1)
    aname[3] = '0' + branch1;
  else 
    aname[3] = ' ';

  aname[4] = '\0';
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

int 
CAtom::getAtomCode(const char * aname)
{
  int elem, place = 0, branch1 = 0, branch2 = 0;

  int n_elem, n_place = 0, n_branch1, n_branch2;
  if (strlen(aname) == 1)
    n_elem = 0;
  else if (strlen(aname) == 2)
    {
      n_elem = 0;
      n_place = 1;
    }
  else if (strlen(aname) == 3)
    {
      if (aname[0] <= '9' && aname[0] >= 0)
	{
	  branch2 = (int)(aname[0] - '0');
	  n_elem = 1;
	  n_place = 2;
	}
      else
	{
	  n_elem = 0;
	  n_place = 1;
	  branch1 = (int)(aname[2] - '0');
	}
    }

  else if (strlen(aname) == 4)
    {
      if (aname[0] <= '9' && aname[0] >= 0)
	{
	  branch1 = (int)(aname[3] - '0');
	  branch2 = (int)(aname[0] - '0');

	  n_elem = 1;
	  n_place = 2;
	}
      else
	{ 
	  branch1 = (int)(aname[2] - '0');
	  branch2 = (int)(aname[3] - '0');

	  n_elem = 0;
	  n_place = 1;
	}
    }
  
  for (int i = 0; i < NUM_ELEM; i++)
    if (aname[n_elem] == ELEMSYM[i]) 
      {
	elem = i;
	break;
      }
  
  if (n_place != 0)
    for (int i = 0; i < NUM_PLACES; i++)
      if (aname[n_place] == PLACESYM[i]) 
	{
	  place = i;
	  break;
	}
  
  return ((elem << ELEM_SHIFT) |
	  (place << PLACE_SHIFT) | 
	  (branch1 << BRANCH1_SHIFT) |
	  (branch2 << BRANCH2_SHIFT));
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/

const CAtom *
AA::operator[](const char * aname) const
{
  int acode = CAtom::getAtomCode(aname);
  vector<CAtom>::const_iterator it = m_atoms.begin();
  for (; it != m_atoms.end(); it++)
    if (it->getCode() == acode)
      return &(*it);
  
  return NULL;
}

/******************************************************************/
/******************************************************************//**
* loadFromPDB
* Inputs: pdb filename for protein, empty vector for amino acid seq
******************************************************************/

void CPDB::loadFromPDB(const char * fname, vector<AA> & aas)
{
  aas.clear();

  ifstream fin(fname);
  if (!fin.is_open())		//!< Checking for file
    {
      cout << "Could not open input file " << fname << endl;
      exit(0);
    }

  char buf[200];
  char aname[5], rname[5];
  int curr_res = -1;
  float x,y,z;
  int i = 0;
  AA aa;
  bool first = true;
  while (!fin.eof())
    {
      fin.getline(buf, 199);
			
			if (strncmp(buf, "ATOM", 4) != 0)  //<! If there is an atom in the current line
				continue;

      CAtom a = readline(buf);					//<! read it in with readLine function, return atom class with name, xyz etc
     
      if (a.getResNum() != curr_res)
			{
				if (!first)
					aas.push_back(aa);
				else 
					first = false;

				aa.clear();
				aa.setType((AA::AACODE) a.getResCode());
				curr_res = a.getResNum();
			}

      aa.insertAtom(a);
    }
  
  aas.push_back(aa);

  /*
  int p = 0;
  vector<AA>::iterator it = aas.begin();
  for (; it != aas.end(); it++)
    {
      vector<ATOM_>::iterator at = it->atoms.begin();
      for (; at != it->atoms.end(); at++)
	cout << at->name << " " << AA_NAMES[it->type] << " " << p << " " 
	     << at->pos[0] << " " << at->pos[1] << " " << at->pos[2] << endl;

	  p++;
    }
  */
}

/******************************************************************/
/******************************************************************//**
* readLine: read a line from a pdb file to obtain information about 
* an atom
******************************************************************/

CAtom
CPDB::readline(const char * buf)
{
  int i = 12;
  while (buf[i] == ' ')
    i++;

  int j = 0;
  char aname[5];
  while (buf[i] != ' ')
    aname[j++] = buf[i++];

  aname[j] = '\0';
  int acode = CAtom::getAtomCode(aname);

  char rname[4];
  strncpy(rname, &(buf[17]),3);
  rname[3] = '\0';
  int rcode = AA::getAACode(rname);

  i = 22;
  while (buf[i] == ' ')
    i++;

  int resnum = atoi(&(buf[i]));
  
  i = 30;
  while (buf[i] == ' ')
    i++;
  float x = atof(&(buf[i]));

  i = 38;
  while (buf[i] == ' ')
    i++;
  float y = atof(&(buf[i]));

  i = 46;
  while (buf[i] == ' ')
    i++;
  float z = atof(&(buf[i]));

  return CAtom(acode, rcode, resnum, x, y, z);
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/
void CPDB::writeToPDB(const char * fname, const vector<AA> & aas)
{
  ofstream fout(fname);
  if (!fout.is_open())
    {
      cout << "Could not output coordinates to: " << fname << endl;
      return;
    }

  char buf[100];
  sprintf(buf,"HEADER    %-70s", fname);
  fout << buf << endl;
  sprintf(buf, "%-80s", "COMPND");
  fout << buf << endl;
  sprintf(buf, "%-80s", "SOURCE");
  fout << buf << endl;
}

/******************************************************************/
/******************************************************************//**
* 
******************************************************************/
		   
void 
CPDB::writeLine(ostream & fout, int index, const char * chainid,
		const CAtom & atom, const CQuat & rot, const CPnt & trans)
{
  fout << "ATOM  ";

  fout.width(5);
  //fout << right;
  fout.setf(ios::right, ios::adjustfield);
  fout << index;

  fout << ' ';

  fout.setf(ios::left, ios::adjustfield);
  char aname[5];
  CAtom::generateName(atom.getCode(), aname);

  fout << aname << ' ';
  fout << AA::getName((AA::AACODE) atom.getResCode()) << ' ';
  fout << chainid;

  fout.setf(ios::right, ios::adjustfield);
  fout.width(4);
  fout << atom.getResNum();

  fout << "    ";
  
  char buf[50];
  CPnt p = rot * atom.getPos() + trans;
  sprintf(buf, "%8.3f%8.3f%8.3f%26s", p.x(), p.y(), p.z()," ");
  fout << buf << endl;;
}

/******************************************************************/
/******************************************************************//**
* Using the amino acid code name (3 or 1 letter), identify the given
* amino acid 
******************************************************************/

AA::AACODE 
AA::getAACode(const char * aa)
{
  for (int i = 0; i < NUM_AAS; i++)
    {
      if (strcmp(aa, AANAME[i]) == 0 ||
						strcmp(aa, AALETTER[i]) == 0)
				return (AACODE) i;
    }

  cout << "bad AA name " << aa << endl;
  exit(0);
  return AA::GLY;
}

