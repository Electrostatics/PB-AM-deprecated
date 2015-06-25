#ifndef _PROTEIN_H_
#define _PROTEIN_H_

#include <map>

#include "pdb.h"
#include "util.h"
#include "mpe.h"

using namespace std;

//!  The protein class
/*! The class that contains all details for a protein and
its multipole expansions, etc */

class CProtein {
public:
	//!  The protein initParameters function
	/*! Function that initializes parameters for a protein class
	        \param kappa the inverse debye length
	        \param dielp the protein dielectric constant
	        \param diels the solvent dielectric constant
	        \param nmol the number of molecules in the system
	        \param rs a scaling factor of the molecule  */
	static void initParameters(REAL kappa, REAL dielp, REAL diels, int nmol,
			REAL rs);

	//!  The protein class constructor
	/*! CProtein class constructor
	        \param fname a character string of filename for input PDB */
	CProtein(const char *fname);

	//!  The protein inCollision function
	/*! Function that determines whether two molecules have collided or not
	        \param mol a vector of CProtein objects
	        \return a boolean indicating whether the molecules have collided
	                (true) */
	bool inCollision(const CProtein &mol);

	// void PDBoutput(ostream & fout, int sindex, const char * chainid,
	//							 const CQuat & rot, const CPnt & trans);
	// void QCDoutput(ostream & fout);

	//!  The protein computeForces function
	/*! Function that computes forces and torques on a group of CProtein
	    objects
	        \param force a vector of forces
	        \param torque a vector of torques */
	static void computeForces(vector<CPnt> &force, vector<CPnt> &torque);

	//!  The protein getAtom function
	/*! Function that returns the desired atom
	        \param resnum an int of the residue number of the desired atom
	        \param acode an int of the atom code of the desired atom */
	const CAtom *getAtom(int resnum, int acode) const;

	//!  The protein getInterfaceAtoms function
	/*! Function that returns a list of atoms at the interface
	        \param P1 a CProtein object
	        \param P2 a second CProtein object that will dock with P1
	        \param pos a vector of positions of atoms at the P!-P2 interface */
	static void getInterfaceAtoms(const CProtein &P1, const CProtein &P2,
			vector<CPnt> &pos);

	//  void computeTransformTo(const CProtein & target, CQuat & q,
	//													CPnt & trans) const;

	//!  The protein isInitialized
	/*! CProtein isInitialized function, check if object's constants
	    have been initialized
	        \return a boolean of whether or not the system has been initialized
	*/
	static bool isInitialized() { return m_initialized; }

	//!  The protein computeDistance
	/*! CProtein computeDistance function
	        \param p1 a CG molecule object
	        \param p2 a second CG molecule object
	        \return a floating point of the distance between the two objects */
	static REAL computeDistance(const CProtein &p1, const CProtein &p2) {
		return (p1.getPosition() - p2.getPosition()).norm();
	}
	//!  The protein translate
	/*! CProtein translate function
	        \param trans a vector of XYZ coordinates to translate the CProtein
	               object by */
	void translate(const CPnt &trans) {
		m_undoCenter = m_center;
		m_center += trans;
	}
	//!  The protein setPos
	/*! CProtein setPos function
	        \param trans a vector of XYZ coordinates to set the center at */
	void setPos(const CPnt &trans) { m_center = trans; }
	//!  The protein rotate
	/*! CProtein rotate function
	        \param rot a quaternion to rotate the CProtein object by */
	void rotate(const CQuat &rot) {
		m_orient = rot * m_orient;
		m_mpe->setOrient(m_orient);
	}
	//!  The protein setOrient
	/*! CProtein setOrient function
	        \param rot a quaterion to set the orientation of the molecule with
	   */
	void setOrient(const CQuat &rot) {
		m_orient = rot;
		m_mpe->setOrient(m_orient);
	}
	//!  The protein untransform
	/*! CProtein untransform function, move protein center back to saved
	    location */
	void untransform() { m_center = m_undoCenter; }

	//!  The protein getNumCharges
	/*! CProtein getNumCharges function
	        \return the number of charges in the protein */
	int getNumCharges() const { return m_chargedAtoms.size(); }
	//!  The protein getCharge
	/*! CProtein getCharge function
	            \param i the ith atom in the protein for which to obtain charge
	            \return a floating point of the charge at atom i */
	REAL getCharge(int i) const { return m_chargedAtoms[i]->getCharge(); }
	//!  The protein getChargePos
	/*! CProtein getChargePos function
	            \param i the ith atom in the protein for which to obtain charge
	                     position
	            \return a CPnt of the position of atom i */
	const CPnt &getChargePos(int i) const {
		return m_chargedAtoms[i]->getPos();
	}

	REAL getRadius() const { return m_rad; }
	REAL getSumCharge() const { return m_sumCharge; }
	const CQuat &getOrientation() const { return m_orient; }
	const CPnt &getPosition() const { return m_center; }
	const vector<CAtom> &getAtoms() const { return m_atoms; }

	int getID() const { return m_id; }

	//! Array of charges for each amino acid in the protein.
	static map<int, REAL> CHARGES[NUM_AAS];
private:
	//!  The protein loadChargeMap function
	/*! Function loads the charges from the file charges_OPLS  */
	static void loadChargeMap();

	//!  The protein computeCenter function
	/*! Compute center of geometry of the protein of interest  */
	CPnt computeCenter();

	//!  The protein computeRadius function
	/*! Compute radius of the coarse-grained protein of interest  */
	void computeRadius();

	//!  The protein getMaxAtomRad function
	/*! Determine which of all of the atoms in the CG simulation has the
	        largest radius.  */
	REAL getMaxAtomRad();

	// Variables
	//! Operator that indicates whether the protein constants are initialized
	static bool m_initialized;
	//! A vector of expansions that are kept for each protein in the system
	static vector<CMPE *> m_exps;
	//! A vector of positions for each protein in the system
	static vector<CPnt *> m_cens;
	//! A vector of protein classes, one for each protein in the system.
	static vector<CProtein *> m_mols;
	//! Operator, designates whether a force calc has been performed  (true=no)
	static bool m_bFirst;

	vector<CAtom> m_atoms; //!< A vector of all the atoms in the protein
	//! A vector of all the charged atoms in the protein
	vector<const CAtom *> m_chargedAtoms;

	CMPE *m_mpe;   //!< An MPE object of the multipole expansion of the object.
	CPnt m_center; //!< A cartesian coordinate object of current position
	//! A saved molecule coordinate, used for undoing a translation
	CPnt m_undoCenter;
	CQuat m_orient;     //!< A quaternion of the current orientation
	CQuat m_undoOrient; //!< A saved orientation, used for undoing a rotation
	REAL m_rad;         //!< The radius of the CProtein object
	REAL m_sumCharge;   //!< The total charge of the CProtein object

	int m_id; //!< The numerical ID of the molecule
};

#endif
