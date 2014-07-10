#ifndef _PDB_H_
#define _PDB_H_

#include <vector>
#include "util.h"

#define NUM_ELEM 5
#define NUM_PLACES 10
#define NUM_AAS 22

using namespace std;

class AA;

//!  The Atom class
/*!
		The class that contains all details for a atom object
*/
class CAtom {
 public:
 
 //!  The Atom class constructor
/*!
		Creates an Atom class. using multiple inputs.  Maps the given atom to
		it's OPLS partial charge.
		\param acode an integer containing the atom code (atomic symbol)
		\param rcode an integer containing the residue name
		\param resnum an integer containing the residue number of atom
		\param x a floating point containing the x coordinate of the atom
		\param y a floating point containing the y coordinate of the atom
		\param z a floating point containing the z coordinate of the atom
*/
  CAtom(int acode_, int rcode_, int resnum_, REAL x, REAL y, REAL z);
	
	 //!  The Atom class constructor
/*!
		Creates an Atom class. using another atom object
		\param a a CAtom object 
*/
  CAtom(const CAtom & a) :
	m_acode(a.m_acode), m_rcode(a.m_rcode), m_resnum(a.m_resnum),
	m_pos(a.m_pos), m_charge(a.m_charge), m_rad(a.m_rad) {}
	
	//!  The Atom class =
  const CAtom & operator=(const CAtom & a)
	{ 
		m_acode = a.m_acode; m_rcode = a.m_rcode; m_resnum = a.m_resnum;
		m_pos = a.m_pos; m_charge = a.m_charge; m_rad = a.m_rad;
	}
	
	static const char ELEMSYM[NUM_ELEM];
  enum ELEMENT {H, C, N, O, S};							//!< Elements
	
  static const char PLACESYM[NUM_PLACES];

		 //!  The Atom class enum PLACE
/*!
		An enumerator for the various positions on a protein that an atom
		may be found.
*/
  enum PLACE {ALPHA = 1, BETA = 2, GAMMA = 3, DELTA = 4, EPSILON = 5, ZETA = 6, 
		ETA = 7, NITRO = 8, TERM = 9};					
  enum BRANCH {BR1 = 1, BR2 = 2, BR3 = 3};
	
  static const int ELEM_MASK, PLACE_MASK, BRANCH1_MASK, BRANCH2_MASK;
  static const int ELEM_SHIFT, PLACE_SHIFT, BRANCH1_SHIFT, BRANCH2_SHIFT;
	
		//!  The Atom class getAtomCode function
/*!
		A function that returns a code from the character string
		\param aname a character string of 4 that specifies atom name
*/
  static int getAtomCode(const char * aname);
	
	//!  The Atom class generateName function
/*!
		A function that generates a 4 character string containing 
		the name of the atom
		\param acode is the integer code of the atom
		\param aname a pointer to a character string that is used
							to identify the atom
*/
  static void generateName(int acode, char * aname);
	
		//!  The Atom class setPos function
/*!
		A function that sets the atom position to xyz coords given
		\param p a CPnt object
*/
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
  int m_acode;				//!< Atom element code
  int m_rcode;				//!< Residue code				
  int m_resnum;				//!< Residue ID number
  CPnt m_pos;					//!< XYZ coordinates of the atom
  REAL m_charge;			//!< Charge of the atom
  REAL m_rad;					//!< Radius of the atom
};  // end CAtom


class CPDB;

//!  The AA class
/*!
		The class that contains all details for an amino acid object
*/
class AA {
 public:
 
  //!  The AA [] operator
	/*!
		Given a character string of the aname, return the atom of interest
		\param aname a four character pointer for atom of interest
		\return a pointer to the CAtom object of interest
*/
  const CAtom * operator[](const char * aname) const;
	
	  //!  The AA [] operator
	/*!
		Given the atom index, return the atom of interest
		\param i an integer of the atom index
		\return the CAtom object of interest
*/
  const CAtom & operator[](int i) const
	{ return m_atoms[i]; }
	
//!  The AA class enum AACODE
/*!
		An enumerator for the amino acids given by their 3 letter codes
*/
  enum AACODE {ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, 
								MET, PHE, PRO, SER, THR, TRP, TYR, VAL, NTR, CTR};


  static const char AANAME[NUM_AAS][4];		//!< 3 Letter codes for each amino acid
  static const char AALETTER[NUM_AAS][2]; //!< 1 Letter codes for each amino acid
	
	//!  The AA class function getAACode
/*!
		Return the enum code for a given amino acid
		\param aa a string of amino acid name
		\return an enum for the amino acid in question
*/
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
  vector<CAtom> m_atoms;				//!< A vector of CAtom objects contained in the given amino acid
  AACODE m_type;								//!< The enum telling what amino acid the object is
}; // end AA

//!  The PDB class
/*!
		The class that contains all details for a PDB object
*/
class CPDB
{
 public:
 
 //!  The PDB function loadFromPDB
/*!
		The function that opens a PDB and reads in the atoms in it, storing
		each in their respective amino acid groups
		\param fname a character string of the filename handle
		\param a vector of amino acids to store information from PDB
*/
  static void loadFromPDB(const char * fname, vector<AA> & aas);
 
 //!  The PDB function writeToPDB
/*!
		The function that writes out a group of amino acids to a PDB
		\param fname a character string of the filename handle
		\param a vector of amino acids to print to PDB
*/	
	static void writeToPDB(const char * fname, const vector<AA> & aas);
	
	
 //!  The PDB function writeLine
/*!
		The function called by writeToPDB to write out each line
		\param fout a character string of the filename handle
		\param int an index of the atom number
		\param chainId a character string of the chain ID
		\param atom a CAtom object of the atom to print
		\param rot a CQuat object of the atom's orientation
		\param trans a CPnt object of the atoms translation
*/
  static void writeLine(ostream & fout, int index, const char * chainid,
													const CAtom & atom, const CQuat & rot, const CPnt & trans);
  
 //!  The PDB function readLine
/*!
		The function called by loadFromPDB to read in a line and
		store as a CAtom class object
		\param buf a string from the PDB file
		\return CAtom object of the atom information contained on the line
*/	
	static CAtom readline(const char * buf);
};  //end CPDB



#endif
