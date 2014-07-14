#ifndef _BD_H_
#define _BD_H_

#include <vector>

#include "util.h"
#include "pdb.h"
#include "protein.h"

using namespace std;

//!  The Brownian dynamics class
/*! The class that contains all details for a BD run, where 1 protein is fixed
		at the center while the other is allowed to translate and rotate, until it has
		escaped the system or it has come into contact with the fixed protein. */

class CBD
{
 public:
 
 //! CBD class constructor
 /*! 
			\param fname1 a constant character pointer to the filename of first protein
			\param fname2 a constant character pointer  to the filename of second protein
			\param kappa a real number of the inverse debye length of the system
			\return an object of the CBD class
	*/
  CBD(const char * fname1, const char * fname2, REAL kappa);
	
	//! An enum of the status of the BD system.
  enum STATUS {
								ESCAPED,  //!< Enum value for escaped protein in BD run
								DOCKED,		//!< Enum value for docked protein in BD run
								RUNNING		//!< Enum value for simulation that isn't complete in BD
	};
	
  static const char STATUS_STRING[3][10]; //!< A string detailing the possible statuses of the system

 //! CBD run function
 /*! 
			Requires no inputs, performs a run on the BD class object
			\return an enum of the BD run's status, whether it be Escaped, docked or running
	*/
  CBD::STATUS run();
	
  REAL computeRate(int numTraj, int nDocked);
	
	//! CBD getRandVec function
	/*! 
			\param std a floating point of random vector length
			\return an CPnt object of random orientation and lenght std
	*/
  static CPnt getRandVec(REAL std)
    { return std*CPnt(normRand(), normRand(), normRand()); }

 private:
  void saveState();
  bool isDocked() const;
  bool escaped(REAL dist) const;
  REAL compute_dt() const; 
  CBD::STATUS makeMove(const CPnt & dR2, const CPnt & dO1, 
		const CPnt & dO2, REAL dt);
  void computeInterfaceVectors(const char * fname1, const char * fname2);
  
  
  CProtein * m_mol1, * m_mol2;			//!< Pointers to the two proteins required for the BD run
  REAL m_Dtr;												//!< Translational diffusion coefficient for the second protein
  REAL m_Dr1;												//!< Rotational diffusion coefficient for first protein
	REAL m_Dr2;												//!< Rotational diffusion coefficient for second protein
  REAL m_maxContDist;								//!< Closest distance that the proteins may be to each other, considered contact distance
  int m_fcnt;
  CPnt m_patch1, m_patch2, m_orth1, m_orth2;
  CQuat m_rot1, m_rot2;
}; // end CBD


///////////////////////////////////////////
///// Inline functions

//!  Escaped function
/*!
		The function that determines whether the second protein has 
		escaped the system or not. Calls on the Protein member function computeDistance
		\param dist a floating point number for the max distance allowed between two proteins
		\return bool of whether the protein has escaped or not, true for escaped
*/
inline bool
CBD::escaped(REAL dist) const
{
  return (CProtein::computeDistance(*m_mol1, *m_mol2) > dist);
}

//!  compute_dt function
/*!
		The function that determines simulation timestep. If the proteins are close enough,
		timestep is shortened by a factor of their distance, otherwise it is defaulted to 1 ps.
		Calls on the Protein member function computeDistance
		\return floating point of timestep in picoseconds.
*/
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

// UNUSED
/*

struct CONTACT
{
  CONTACT(const CAtom * a1_, const CAtom * a2_, REAL dist_) :
    a1(a1_), a2(a2_), dist(dist_) {}
  const CAtom * a1, * a2;
  REAL dist;
};

// In the private portion of the BD class:

  vector<CONTACT> m_contacts;
	void loadContacts(const char * fconts);
	
*/