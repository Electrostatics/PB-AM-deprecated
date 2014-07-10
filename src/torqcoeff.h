#ifndef _TORQCOEFF_H_
#define _TORQCOEFF_H_

#define X 0
#define Y 1
#define Z 2

#include "tricoeff.h"

	//!  The CTorqCoeff class 
/*!
		The CTorqCoeff class contains information about 
		torque coefficients, related to TriCoeff, and used for
		computing torque coefficients in a molecule, the H variables  
*/
class CTorqCoeff : public CTriCoeff
{
 public:
  CTorqCoeff(int p = 0, int res = 0) : CTriCoeff(p, res) {} 
  CTorqCoeff(const CTriCoeff & G) : CTriCoeff(G) {} 
	
//!  The CTorqCoeff class constructor 
/*!
		The CTorqCoeff construction using the inputs:
		\param charges a vector of charges within the molecule
		\param pos a vector of charge positions in cartesian coords
		\param p an integer of number of poles
		\param rad a floating point number of the radius of the CG sphere. 
*/
  CTorqCoeff(const vector<REAL> & charges, const vector<CPnt>&  pos, 
	     int p, REAL rad);

};

#endif
