#include <fstream>
#include <vector>

#include "BD.h"
#include "protein.h"
#include "mcoeff.h"

/******************************************************************/
/******************************************************************//**
* Initializing constants
******************************************************************/
#define b_DIST 100.0											//!< Initial distance between 2 proteins for BD run
#define q_DIST 500.0											//!< Distance for molecules to be considered escaped
#define f_DIST 100.0											//!< Cutoff for protein force interactions

#define DIELECTRIC_WATER 78.0							//!< The dielectric constant of water
#define DIELECTRIC_PROT 4.0								//!< Dielectric constant of the protein

#define SALT_CONCENTRATION  0.0500				//!< [ Molar ]  
#define KAPPA (sqrt(SALT_CONCENTRATION)/3.04)	//!< Inverse debye length 

#define COULOMB_CONSTANT (8.988e9)				//!< [ N*m^2/C^2 ]
#define ELECTRON_CHARGE (1.60217733e-19)	//!<  [ coulombs ]
#define AVOGADRO_NUM (6.02209e23)				
#define KCAL 4184.0												//!<  [ 1 kCal = 4184 Joules ]
#define ANGSTROM (1e-10)									//!<  [ 1A = 1e-10 Meters ]
#define PICO_SEC (1e-12)									//!<  [ 1 ps = 1e-12 s ]	
#define COUL_K 332.141144

#define GAMMA (0.01*100/1000*ANGSTROM*PICO_SEC)
#define Kb (1.380658e-23)									//!<  [ m^2 kg/ s^2 / K ] 
#define TEMP 298.0												//!<  [ Kelvin ]
#define KbT (1.98718E-03*TEMP)						//!< (TEMP*Kb*AVOGADRO_NUM/KCAL)
#define IKbT (1.0/KbT)										//!<  1/KbT

#define TOL 2.5

#define PATCH_ANGLE 6.0
#define ROTATE_ANGLE 20.0
#define PATCH_SIZE (cos(PATCH_ANGLE*M_PI/180.0))
#define ROTATE_SIZE (cos(ROTATE_ANGLE*M_PI/180.0))

const char CBD::STATUS_STRING[3][10] = {"ESCAPED", "DOCKED", "RUNNING"};
int scount;
char temp_file[100];


/******************************************************************/
/******************************************************************//**
* Initialize CBD class
*****************************************************************/
CBD::CBD(const char * fname1, const char * fname2, REAL kappa) : m_fcnt(0)
{
  if (!CProtein::isInitialized())
		// Initialize the system parameters, w/ 2 molecules, and RS=lambda in paper=30
    CProtein::initParameters(kappa, DIELECTRIC_PROT, DIELECTRIC_WATER, 2, 30);

  m_mol1 = new CProtein(fname1); // MOLECULE 1 = Barnase
  m_mol2 = new CProtein(fname2); // MOLECULE 2 = Barstar
  m_maxContDist = m_mol1->getRadius() + m_mol2->getRadius() + 1.0;	// Closest distance that the proteins may be to each other
 
  m_Dtr = 0.030;												// Translational diffusion coefficient
  m_Dr1 = 4e-5; 												// Rotational diffusion coefficient barnase
	m_Dr2 = 4.5e-5;												// Rotational diffusion coefficient barstar

  computeInterfaceVectors(fname1, fname2);

}

/******************************************************************/
/******************************************************************//**
* run: Runs a BD simulation for 2 proteins.  Initializes one 
*			at the origin and the other at some distance b_DIST away with
*			random orientation. While the status is RUNNING, the system
*			computes forces and then moves the 2nd protein and rotates 
*     both while keeping the first fixed at the origin.  It prints 
*			out the position of each molecule, their distance and the 
*     force and torques of the system every 1000 steps.  
******************************************************************/
CBD::STATUS
CBD::run()
{
  ofstream fout(temp_file);
  m_mol1->setOrient(CQuat::chooseRandom());
  m_mol1->setPos(CPnt());
  m_mol2->setPos(b_DIST*randOrient());
  m_mol2->setOrient(CQuat::chooseRandom());
  
  vector<CPnt> force(2), torque(2);
  vector<REAL> pot;
  REAL fact = IKbT*COUL_K/DIELECTRIC_WATER;										// Conversion from internal units
  REAL srad = m_mol1->getRadius() + m_mol2->getRadius();			// Sum of protein radii
  STATUS status = RUNNING;
  scount = 0;
  CPnt dR2, dO1, dO2; 
  while (status == RUNNING)																		// 2 steps to BD run, compute dt and force computation
    {
      REAL dt = compute_dt();     
      REAL dist = CProtein::computeDistance(*m_mol1, *m_mol2);
      if (dist < f_DIST)																		 // If two proteins are within cutoff, compute forces
			{
				CProtein::computeForces(force, torque);
	  
				dR2 = (m_Dtr*dt*fact)*force[1];											// Move 2nd protein
				dO1 = (m_Dr1*dt*IKbT*COUL_K)*torque[0];							// Rotate both proteins
				dO2 = (m_Dr2*dt*IKbT*COUL_K)*torque[1];
		}
    else
		{
			dR2.zero();
			dO1.zero();
			dO2.zero();
		}

    if (scount % 1000 == 0)																	// Print out details every 1000 steps
		{
			CPnt t =  m_mol2->getPosition() - m_mol1->getPosition();
			fout << scount << ") " << m_mol1->getPosition() 
					<< m_mol2->getPosition()	<< "-> " << dist << endl;
			fout << force[1] << " " << torque[1] << endl;
		}

    status = makeMove(dR2, dO1, dO2, dt);										// Move system with given dr, d angle and dt
    scount++;
	}

  fout << CBD::STATUS_STRING[status] << endl;
  return status;
} // end CBD::run

/******************************************************************/
/******************************************************************//**
* makeMove
******************************************************************/
CBD::STATUS
CBD::makeMove(const CPnt & dR2, const CPnt & dO1, 
							const CPnt & dO2, REAL dt) 
{
  REAL bdist = CProtein::computeDistance(*m_mol1, *m_mol2);
  int c = 0;
  while(true)
	{
		c++;
		if (c > 500)
			cout << "stuck inside loop" << endl;
		CPnt dR2_ = dR2 + CBD::getRandVec(sqrt(2*dt*m_Dtr));
		m_mol2->translate(dR2_);
		
		if (!m_mol1->inCollision(*m_mol2))
			break;
		else
			m_mol2->untransform();
	}
	
  if (escaped(q_DIST))
    return ESCAPED;
	
  CPnt dO1_ = dO1 + CBD::getRandVec(sqrt(2*dt*m_Dr1));
  CPnt dO2_ = dO2 + CBD::getRandVec(sqrt(2*dt*m_Dr2));
  
  CQuat Q1(dO1_, dO1_.norm());
  CQuat Q2(dO2_, dO2_.norm());
  
  REAL adist = CProtein::computeDistance(*m_mol1, *m_mol2);
  if (adist > f_DIST)
	{
		m_rot1 = Q1*m_rot1;
		m_rot2 = Q2*m_rot2;
	}
  else
	{
		if (bdist > f_DIST)
		{
			m_mol1->rotate(m_rot1);
			m_mol2->rotate(m_rot2);
			
			m_orth1 = m_rot1 * m_orth1;
			m_patch1 = m_rot1 * m_patch1;
			m_orth2 = m_rot2 * m_orth2;
			m_patch2 = m_rot2 * m_patch2;
			
			m_rot1.identity();
			m_rot2.identity();
		}
		else
		{
			m_mol1->rotate(Q1);
			m_mol2->rotate(Q2);
			
			m_orth1 = Q1 * m_orth1;
			m_patch1 = Q1 * m_patch1;
			m_orth2 = Q2 * m_orth2;
			m_patch2 = Q2 * m_patch2;
		}
	}
	
  if (isDocked())
    return DOCKED;
  else 
    return RUNNING;
}

/******************************************************************/
/******************************************************************//**
* isDocked function to determine whether or not the moving protein has 
docked on the other
******************************************************************/
bool 
CBD::isDocked() const
{
  REAL d = CProtein::computeDistance(*m_mol1, *m_mol2);
  //cout << d << " " << m_maxContDist << endl;
  if (d > m_maxContDist)
    return false;

  CPnt v = (m_mol2->getPosition() - m_mol1->getPosition()).normalize();
  REAL a = dot(v, m_patch1);
  //cout << a << " " << PATCH_SIZE << endl;
  if (a < PATCH_SIZE)
    return false;

  a = -dot(v, m_patch2);
  //cout << a << " " << PATCH_SIZE << endl;
  if (a < PATCH_SIZE)
    return false;

  a = dot(m_orth1, m_orth2);
  //cout << a << " " << ROTATE_SIZE << endl;
  if (a < ROTATE_SIZE)
    return false;

  return true;
}


/******************************************************************/
/******************************************************************//**
* computeInterfaceVectors: compute patches - vector between proteins
* and compute orthogonal vectors to the two vectors
******************************************************************/
void 
CBD::computeInterfaceVectors(const char * fname1, const char * fname2)
{
  CPnt mean;
  vector<CPnt> pnts;
  
  CPnt d = (m_mol2->getPosition() - m_mol1->getPosition()).normalize();  //<! Normalized vector between 2 mols
  m_patch1 = d; 
  m_patch2 = -d;   
  m_orth1 = (CPnt(-d.z(), 0.0, d.x())).normalize(); // (cross(m_patch1, m_patch2)).normalize();
  m_orth2 = m_orth1;

}

/******************************************************************/
/******************************************************************//**
* computeRate
******************************************************************/
REAL
CBD::computeRate(int numTraj, int nDocked)
{
  REAL cnst = 4.0*M_PI*m_Dtr; 
  REAL maxR = q_DIST;
  REAL dR = 10.0;

  REAL k_b = 4*M_PI*m_Dtr*b_DIST ;
  REAL k_q = 4*M_PI*m_Dtr*q_DIST;
  cout << k_b << " " << k_q << endl;

  REAL delta = ((REAL)nDocked)/numTraj;
  REAL gamma = k_b/k_q;

  REAL K = k_b*delta/(1 - (1 - delta)*gamma);  

  cout << "The relative rate is: " << K << endl;
  return K;
}

/******************************************************************/
/******************************************************************//**
* approxBBall function to approximate.
\param V a vector of xyz coordinates for charges
\param cent an xyz coord for the CG sphere center
\param rad a floating point of the radius of the CG sphere
******************************************************************/
void 
approxBBall(const vector<CPnt> & V, CPnt & cen, REAL & rad)
{
  int n = V.size();
  REAL rad2; // radius squared
  REAL xmin, xmax, ymin, ymax, zmin, zmax;  // bounding box extremes   
  int  Pxmin, Pxmax, Pymin, Pymax, Pzmin, Pzmax;  // index of V[] at box extreme
	
  // find a large diameter to start with.    
  // first get the bounding box and V[] extreme points for it   
  xmin = xmax = V[0].x(); ymin = ymax = V[0].y();  zmin = zmax = V[0].z(); 
  Pxmin = Pxmax = Pymin = Pymax = Pzmin = Pzmax = 0;
	
  for (int i=1; i<n; i++) 
	{       
		if (V[i].x() < xmin) 
		{ xmin = V[i].x();  Pxmin = i; }   
		else if (V[i].x() > xmax) 
		{ xmax = V[i].x(); Pxmax = i; }
		
		if (V[i].y() < ymin) 
		{ ymin = V[i].y(); Pymin = i; }
		else if (V[i].y() > ymax) 
		{ ymax = V[i].y(); Pymax = i; }
		
		if (V[i].z() < zmin) 
		{ zmin = V[i].z(); Pzmin = i; }
		else if (V[i].z() > zmax) 
		{ zmax = V[i].z(); Pzmax = i; }
	} 
	
  // select the largest extent as an initial diameter for the ball
  CPnt dVx = V[Pxmax] - V[Pxmin]; // diff of Vx max and min   
  CPnt dVy = V[Pymax] - V[Pymin]; // diff of Vy max and min
  CPnt dVz = V[Pzmax] - V[Pzmin]; // diff of Vz max and min
  REAL dx2 = dVx.normsq(), dy2 = dVy.normsq(), dz2 = dVz.normsq(); // diffs squared
	
  if (dx2 >= dy2 && dx2 >= dz2) // x direction is largest extent
	{         
		cen = V[Pxmin] + (0.5*dVx);  // Center = midpoint of extremes
		rad2 = (V[Pxmax] - cen).normsq();   // radius squared   
	}
  else if (dy2 >= dz2) // y direction is largest extent
	{                                
		cen = V[Pymin] + (0.5*dVy);   // Center = midpoint of extremes
		rad2 = (V[Pymax] - cen).normsq();    // radius squared
	} 
  else  // z direction is largest extent
	{                                
		cen = V[Pzmin] + (0.5*dVz);   // Center = midpoint of extremes
		rad2 = (V[Pzmax] - cen).normsq();    // radius squared
	} 
	
  rad = sqrt(rad2);
	
  // now check that all points V[i] are in the ball    
  // and if not, expand the ball just enough to include them
  CPnt dV;   
  REAL dist, dist2;
  for (int i=0; i<n; i++) 
	{
		dV = V[i] - cen;
		dist2 = dV.normsq();
		if (dist2 <= rad2)    // V[i] is inside the ball already
			continue;
		
		// V[i] not in ball, so expand ball to include it
		dist = sqrt(dist2); 
		rad = (rad + dist) / 2.0;   // enlarge radius just enough
		rad2 = rad * rad;
		cen = cen + ((dist-rad)/dist) * dV;   // shift Center toward V[i]
	}
  
  cout << "Ball: " << cen << " " << rad << endl;
  return;
} // end approxBBall

/******************************************************************/
/******************************************************************//**
* readqcd
******************************************************************/
void readqcd(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     REAL & rad)
{
  ifstream fin(fname);
  char buf[100], temp;
  REAL sum = 0.0;
	
  pnt.clear();
  ch.clear();
  bool first = true;
  fin.getline(buf,99);
  while (!fin.eof())
	{
		double x,y,z,c,r;
		sscanf(buf, "ATOM\t1\tLTL\tTST\t%lf %lf %lf %lf %lf", 
					 &x, &y, &z, &c, &r);
		
		if (first)
		{
			rad = r;
			first = false;
		}
		else
		{
			pnt.push_back(CPnt(x,y,z));
			ch.push_back(c);
			sum += c;
		}
		fin.getline(buf,99);     
	}
	
  cout << fname << ": charges: " << ch.size() << ", net charge: " 
	<< sum << ", radius: " << rad << endl; 
}

/******************************************************************/
/******************************************************************//**
* writemdp: a function to write out a PQR file type.  
\param ifname a character string that contains the input file name
\param ofname a character name that contains the output file name
\param pnt a vector of charge positions
\param cen a position of the CG sphere center
\param sp a boolean that tells function whether to print out the CG
				sphere in the PQR file or not
\param chid a character of the chain ID of the molecule
******************************************************************/
void writemdp(const char * ifname, const char * ofname, 
							const vector<CPnt> & pnt, const CPnt& cen, bool sp,
							char chid)
{
  ifstream fin(ifname);
  ofstream fout(ofname, ios::app);
  if (!fin.is_open())
	{
		cout << "Could not open file " << ifname << endl;
		exit(0);
	}
	
  if (!fout.is_open())
	{
		cout << "Could not open file " << ofname << endl;
		exit(0);
	}
	
  char buf[100], temp[30];
  fin.getline(buf,99);
  int i = 0;
  while (!fin.eof())
	{
		if (strncmp(&(buf[0]),"ATOM",4) == 0)	  
		{
			buf[21] = chid;
			if (strncmp(&(buf[17]),"DUM",3) == 0)
	    {
	      if (sp)
				{
					CPnt P = cen;
					sprintf(temp,"%8.3f%8.3f%8.3f", P.x(), P.y(), P.z());
					strncpy(&(buf[30]), temp, 24);
					fout << buf << endl;
				}
	    }
			else
	    {
	      CPnt P = pnt[i] + cen;
	      sprintf(temp,"%8.3f%8.3f%8.3f", P.x(), P.y(), P.z());
	      strncpy(&(buf[30]), temp, 24);
	      fout << buf << endl;
	      i++;
	    }
		}
		
		fin.getline(buf,99);     
	}
}

/******************************************************************/
/******************************************************************//**
* readmdp: function that reads in an MDP file for a given molecule. 
The file should contain the PQR data of a molecule and a dummy molecule
that approximates the CG version of the molecule.
	\param fname a character pointer to the filename
	\param pnt a vector of xyz coordinates for the charges within the molecule
	\param ch a vector of floating point numbers for the charges
	\param rad a floating point number that stores the radius of the molecule
	\param cen a vector of xyz coordinates for the molecule's center of geom
******************************************************************/
void readmdp(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
						 REAL & rad, CPnt & cen)
{
  ifstream fin(fname);
  if (!fin.is_open())					// Check for file
	{
		cout << "Could not open file " << fname << endl;
		exit(0);
	}
	
  char buf[100], temp[10];
  REAL sum = 0.0;							// Total charge of the molecule
	
  pnt.clear();
  ch.clear();
  fin.getline(buf,99);
  
  while (!fin.eof())
	{
		double x,y,z,c,r;
		if (strncmp(&(buf[0]),"ATOM",4) == 0)
		{
			sscanf(&(buf[31]), "%lf %lf %lf", &x, &y, &z);
			strncpy(temp, &(buf[63]), 6);
			temp[6] = 0;
			if (strcmp(temp, "      ") == 0)
				r = 0.0;
			else
				sscanf(temp, "%lf", &r);
			sscanf(&(buf[55]), "%lf", &c);
			
			if (strncmp(&(buf[17]),"DUM",3) == 0)
	    {
	      rad = r;
	      cen = CPnt(x,y,z);
	    }
			else //if (c != 0.0)
	    {
	      pnt.push_back(CPnt(x,y,z));
	      ch.push_back(c);
	      sum += c;
	    }
		}
		
		fin.getline(buf,99);     
	}
	
  cout << fname << ": charges: " << ch.size() << ", net charge: " 
	<< sum << ", radius: " << rad << ", cen: " << cen << endl; 
}

/******************************************************************/
/******************************************************************//**
* buildUnitCell, function used to build a unit cell for use in 
infinite grid calculation.  
	\param ifname a character string describing the input file name
	\param ofname a character string describing the desired output file name
	\param mpe a vector of multipole expansions
	\param cen a vector of sphere centers
	\params bSave a boolean operator indicating whether or not to write 
						out grid
	\params stretch a floating point of how much to stretch the center by ??
******************************************************************/
void 
buildUnitCell(const char * ifname, const char * ofname, vector<CMPE*> & mpe, 
							vector<CPnt*> & cen, bool bSave, REAL stretch)
{
  REAL rad;
  vector<CPnt> V;
  vector<REAL> ch;
  CPnt cn, cn2;
	
  if (bSave)
	{
		ofstream fout(ofname);
		fout.close();
	}
  readmdp(ifname, V,ch,rad,cn);  // read in the molecule from PQR file, save atom point xyzs in V
  approxBBall(V, cn, rad);       // run through all the point charges and build a sphere large 
																 // enough to enclose all of them.  it has a radius of rad
	
  REAL lx = 79.1, ly = 79.1, lz = 37.9;
  REAL sx = 0.012642, sy = 0.012642, sz = 0.026385;
  CPnt t;
	
  for (int i = 0; i < V.size(); i++)			// for all point charges, move their position by a scaling factor of s
    V[i] = CPnt(V[i].x()*sx, V[i].y()*sy, V[i].z()*sz);
  cn = CPnt(cn.x()*sx, cn.y()*sy, cn.z()*sz); // and scale the center by that as well
  vector<CPnt> U(V.size());
	
	
	// create a grid of 8 points and write out each one
  // ASU #1
  t = CPnt(0.0, 0.0, 0.0);
  cn2 = CPnt(cn.x()*lx, cn.y()*ly, cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(V[i].x()*lx, V[i].y()*ly, V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'A');
  mpe.push_back(new CMPE(ch, U, rad, 0, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #2
  t = CPnt(lx, ly, 0.5*lz);
  cn2 = CPnt(-cn.x()*lx, -cn.y()*ly, cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(-V[i].x()*lx, -V[i].y()*ly, V[i].z()*lz) + (t - cn2);
  if (bSave) 
		writemdp(ifname, ofname, U, cn2, false, 'B');
  mpe.push_back(new CMPE(ch, U, rad, 1, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #3
  t = CPnt(0.5*lx, 0.5*ly, -0.25*lz);
  cn2 = CPnt(-cn.y()*lx, cn.x()*ly, cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(-V[i].y()*lx, V[i].x()*ly, V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'C');
  mpe.push_back(new CMPE(ch, U, rad, 2, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #4
  t = CPnt(0.5*lx, 0.5*ly, 0.25*lz);
  cn2 = CPnt(cn.y()*lx, -cn.x()*ly, cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(V[i].y()*lx, -V[i].x()*ly, V[i].z()*lz) + (t- cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'D');
  mpe.push_back(new CMPE(ch, U, rad, 3, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #5
  t = CPnt(0.5*lx, 0.5*ly, 0.75*lz);
  cn2 = CPnt(-cn.x()*lx, cn.y()*ly, -cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(-V[i].x()*lx, V[i].y()*ly, -V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'E');
  mpe.push_back(new CMPE(ch, U, rad, 4, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #6
  t = CPnt(0.5*lx, 0.5*ly, 1.25*lz);
  cn2 = CPnt(cn.x()*lx, -cn.y()*ly, -cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(V[i].x()*lx, -V[i].y()*ly, -V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'F');
  mpe.push_back(new CMPE(ch, U, rad, 5, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #7
  t = CPnt(0.0, 0.0, lz);
  cn2 = CPnt(cn.y()*lx, cn.x()*ly, -cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(V[i].y()*lx, V[i].x()*ly, -V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'G');
  mpe.push_back(new CMPE(ch, U, rad, 6, 15));
  cen.push_back(new CPnt(stretch*cn2));
	
  // ASU #8
  t = CPnt(lx, ly, 0.5*lz);
  cn2 = CPnt(-cn.y()*lx, -cn.x()*ly, -cn.z()*lz) + t;
  for (int i = 0; i < V.size(); i++)
    U[i] = CPnt(-V[i].y()*lx, -V[i].x()*ly, -V[i].z()*lz) + (t - cn2);
  if (bSave)
    writemdp(ifname, ofname, U, cn2, false, 'H');
  mpe.push_back(new CMPE(ch, U, rad, 7, 15));
  cen.push_back(new CPnt(stretch*cn2));
} //end buildUnitCell


/******************************************************************/
/******************************************************************//**
* buildSystem: a function that builds a system of num molecules 
								separated equidistantly in space
*			\param ifname a character pointer that has the name of an
							file of type MDP
			\param num an int that details the number of centers to 
							include for a system of multiple molecules
			\param dist a floating point number, describes the distance
							desired between molecules
			\param bSave a boolean that indicates whether or not the user
							desires to save the configuration
			\param mpe a vector of multipole expansion objects for use
							in mutual polarization
			\param cen a vector of positions to hold the position of each
							molecule in the system
			\return a floating point number that returns the radius of 
							the molecule input into the system
******************************************************************/
REAL buildSystem(const char * ifname, int num, REAL dist, bool bSave,
								 vector<CMPE*> & mpe, vector<CPnt*> & cen)
{
  REAL rad1;						// size of the molecule
  vector<CPnt> p1;			// vector of positions of point charges in the molecule
  vector<REAL> ch1;			// vector of point charges in the molecule
  CPnt cn;							// original position of center of geometry of molecule in space
	
  readmdp(ifname, p1,ch1,rad1,cn);
  REAL fact = rad1 + dist/2.0;
	
	if (num == 1)
	{
		cen.push_back(new CPnt(CPnt(0.0,0.0,0.0)));
	}
  else if (num == 2)
	{
		cen.push_back(new CPnt(fact*CPnt(0.0,0.0,1.0)));
		cen.push_back(new CPnt(fact*CPnt(0.0,0.0,-1.0)));
	}
  else if (num == 4)
	{
		cen.push_back(new CPnt(fact*CPnt(0.0,2/sqrt(3.0),0.0)));
		cen.push_back(new CPnt(fact*CPnt(-1.0,-1/sqrt(3.0),0.0)));
		cen.push_back(new CPnt(fact*CPnt(1.0,-1/sqrt(3.0),0.0)));
		cen.push_back(new CPnt(fact*CPnt(0.0,0.0,sqrt(8.0/3.0))));
	}
  else if (num == 6)
	{
		cen.push_back(new CPnt(fact*CPnt(sqrt(2.0),0.0,0.0)));
		cen.push_back(new CPnt(fact*CPnt(-sqrt(2.0),0.0,0.0)));
		cen.push_back(new CPnt(fact*CPnt(0.0,sqrt(2.0),0.0)));
		cen.push_back(new CPnt(fact*CPnt(0.0,-sqrt(2.0),0.0)));
		cen.push_back(new CPnt(fact*CPnt(0.0,0.0,sqrt(2.0))));
		cen.push_back(new CPnt(fact*CPnt(0.0,0.0,-sqrt(2.0))));
	}
	else if (num == 8)
	{
		cen.push_back(new CPnt(fact*CPnt(1.0,1.0,1.0)));
		cen.push_back(new CPnt(fact*CPnt(-1.0,1.0,1.0)));
		cen.push_back(new CPnt(fact*CPnt(1.0,-1.0,1.0)));
		cen.push_back(new CPnt(fact*CPnt(-1.0,-1.0,1.0)));
		cen.push_back(new CPnt(fact*CPnt(1.0,1.0,-1.0)));
		cen.push_back(new CPnt(fact*CPnt(-1.0,1.0,-1.0)));
		cen.push_back(new CPnt(fact*CPnt(1.0,-1.0,-1.0)));
		cen.push_back(new CPnt(fact*CPnt(-1.0,-1.0,-1.0)));
	}
	else
	{
		cout << "wrong number of molecules: " << num << endl;
		exit(0);
	}
	
	for (int i = 0; i < p1.size(); i++)			// moving the positions of each charge
    p1[i] -= cn;																// to be WRT the center of the mol
	
	char ofname[100], ofname_sp[100];
	if (bSave)																// print out molecules positions if desired
	{
		int l = strlen(ifname);
		
		strcpy(ofname,ifname);
		strcpy(ofname_sp,ifname);
		sprintf(&(ofname[l-7]), "_%dP.mdp",num);
		sprintf(&(ofname_sp[l-7]), "_%dP_sp.mdp",num);
		
		char rm[200];
		sprintf(rm,"rm -f %s", ofname);
		system(rm);
		sprintf(rm,"rm -f %s", ofname_sp);
		system(rm);
		
		cout << "saving files:" << endl;
		cout << "\t" << ofname << endl;
		cout << "\t" << ofname_sp << endl;
	}
	
	// Initialize parameters for the system, using hard-coded values,
	// salt conc = 0.05 M, T=298 K, eps_p = 4, eps_s = 78 etc.
	CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, num, rad1);
	for (int n = 0; n < num; n++)
	{
		CQuat Q = (num > 2 ? CQuat::chooseRandom() : CQuat());		// Chose random orientation if 3+ mol. in system
		vector<CPnt> pr(p1.size());
		for (int i = 0; i < p1.size(); i++)
		{
			pr[i] = Q*p1[i];																				// Rotate mol i by random orientation
			if (pr[i].norm() > rad1)
				cout << "ERROR" << endl;
		}
		
		CMPE * pm = new CMPE(ch1, pr, rad1, n, 15);								// Make MPE for molecule i with 15 poles
		mpe.push_back(pm);
		if (bSave)																								// If we wish to save configurations
		{																													// Write out a file without the CG sphere
			writemdp(ifname, ofname, pr, *(cen[n]), false, ' ');		// And one containing the CG sphere
			writemdp(ifname, ofname_sp, pr, *(cen[n]), true, ' ');
		}
	}
	
	cout << "Building system is complete" << endl;
	return rad1;
} // end buildSystem

/******************************************************************/
/******************************************************************//**
* Building a grid centered at ( 0, 0, 0 ) of bounded from -bd to +bd.
	
	\param ifname a character string of an input filename
	\param num an integer of the number of molecules in the system
	\param rad a floating point of the radius of the molecules
	\param mpe a vector of multipole expansions, one for each molecule
	\param cen a vector of molecule positions 
	\param fact a floating point of the conversion factor for the system
******************************************************************/
void buildGrid(const char * ifname, int num, REAL rad,
							 const vector<CMPE*> & mpe, const vector<CPnt*> & cen, REAL fact)
{
  cout << "building grid" << endl;
  REAL scale = 2.0;													// Scaling factor
//  int gsize = 301;													// Grid edge length
  int gsize = 151;													// Grid edge length
  int gcen = gsize/2 + 1;
  int bd = gcen-1;
  float * V = new float[gsize*gsize*gsize];	// Volume of the grid cell
	
  REAL iscale = 1.0/scale;									// 1/Scaling factor
  REAL r = rad*rad*scale*scale;							// Scaled radius, r = ract^2 / 2^2
	
  int gsizeq = gsize*gsize;									
  cout << "here" << endl;
  for (int i = -bd; i <= bd ; i++)
	{
		for (int j = -bd; j <= bd ; j++)
			for (int k = -bd; k <= bd ; k++)			// For all 3 Dimensions, i, j and k
			{
				bool cont = false;
				
				for (int n = 0; n < num; n++)				// for each particle, compute the distance^2 between
	      {																		// current pos. and 2*center of each particle in the sys
					REAL dis_sq = (CPnt(i,j,k) -(scale*(*cen[n]))).normsq();
					if (dis_sq < r)										// If the distance is less that scaled rad^2, then
					{																	// The point is within the molecule
						cont = true;
						break;
					}
	      }
				
				if (cont)														// If the point is within the molecule, then the potential there is zero
	      {
					V[(k+bd)*gsizeq+(j+bd)*gsize+(i+bd)] = 0.0;
					continue;
	      }
																						// otherwise, the potential can be computed and stored
				CPnt P(i,j,k);
				V[(k+bd)*gsizeq+(j+bd)*gsize+(i+bd)] = 
	      (float) CMPE::computePotAt(mpe, cen, iscale*P)*fact;
			}
		
		cout << i << "..";
		cout.flush();
	}
	
  cout << endl << "done" << endl;
  char ofname[100];
  strcpy(ofname,ifname);
  int l = strlen(ifname);
  sprintf(&(ofname[l-7]), "_%dP.gmpe",num);

	ofstream fd;
	fd.open (ofname);
	
	  for (int i = -bd; i <= bd ; i++)
		{
			for (int j = -bd; j <= bd ; j++)
			{
				for (int k = -bd; k <= bd ; k++)			// For all 3 Dimensions, i, j and k
				{
					fd << V[(k+bd)*gsizeq+(j+bd)*gsize+(i+bd)] << "  ";
				}
				fd << endl;
			}
			fd << "\n\n";
		}
		fd.close();
		return;
	
//  FILE * fd = fopen(ofname, "w+");
//  fwrite(V, sizeof(float), gsize*gsize*gsize, fd);
//	fclose( fd );
} // end buildGrid

/******************************************************************/
/******************************************************************//**
* Computing the perturbation
******************************************************************/
void perturb(int ct, int num, ofstream & fout, REAL Dtr, REAL Dr, REAL dt,
						 vector<CMPE*> & mpe, vector<CPnt*> & cen)
{
  CPnt per[num], rot;
  CQuat Q0[num];
  for (int j = 0; j < num; j++)
    Q0[j] = mpe[j]->getOrient();
	
  for (int k = 0; k < ct; k++)
	{
		cout << "------ " << k << " -------" << endl;
		for (int j = 0; j < num; j++)
		{
			per[j] = CBD::getRandVec(sqrt(2*dt*Dtr));
			*(cen[j]) += per[j];
			
			rot = CBD::getRandVec(sqrt(2*dt*Dr));
			CQuat Q(rot, rot.norm());
			mpe[j]->setOrient(Q*Q0[j]);
  	}
		
		for (int i = 0; i < cen.size(); i++)
			cout << *(cen[i]);
		cout << endl;
		
		CMPE::updateSolve(mpe, cen);

		cout << "Done in perturb of setting orient " << endl;
		for (int i = 0; i < num; i++)
			mpe[i]->saveUndo();

		
		//CMPE::computeForce(mpe, cen, force, torque);
		
	/*	fout << CMPE::npol << " " << CMPE::npol_t << " ";
		for (int i = 0; i < num; i++)
			fout << mpe[i]->ngpol << " " << mpe[i]->ngpol_t << " ";
		fout << endl;
	*/	
	
					
		CMPE::undoXForms();
		for (int i = 0; i < num; i++)
		{
			mpe[i]->undo();
			*(cen[i]) -= per[i];
		}
	}
} // end perturb

/******************************************************************/
/******************************************************************//**
* Main for BD simulation of Barnase/Barstar
* \param a floating point number of salt concentration
* \param fout an output file
* \param ind an integer for temporary file index
******************************************************************/
int main1(int argc, char ** argv)
{
	if ( argc != 7 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec sim [protein 1] [protein 2] [Salt conc] [outfile] [temp file #]" << endl;
		exit(0);
	}
  seedRand(-1);

  REAL kappa = sqrt(atof(argv[4]))/3.04;				// Inverse debye length
 // CBD bd("barnaseh.pdb", "barstarh.pdb",kappa);
  CBD bd(argv[2],argv[3],kappa);
  
  ofstream fout(argv[5]);
  sprintf(temp_file, "temp%d.out", atoi(argv[6]));

  int cdock = 0;
  for (int i = 0; i < 50000; i++)
    {
      CBD::STATUS status = bd.run();
      if (status == CBD::DOCKED)
				cdock++;

      fout << CBD::STATUS_STRING[status] << " " << (int)status 
	   << " " << scount << endl;
    }
  fout << cdock << " " << argv[1] << endl;
  return 0;
} // end main1

/******************************************************************/
/******************************************************************//**
* Main for computing energies, forces and torques on a replication
		of a single molecule places equidistantly in a system.  
* \param ifname an input file name
* \param num an integer number of molecules to include in the system. 
			Values are: 2, 4, 6 or 8						
* \param dist a floating point distance for molecule placement
******************************************************************/
int main2(int argc, char ** argv)
{
	if ( argc != 5 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec slv [PQR file] [2, 4, 6 or 8 molecules] [distance between molecules]" << endl;
		exit(0);
	}
	
  cout << "SOLVE" << endl;
  seedRand(-1);
  cout.precision(5);
 
  const char * ifname = argv[2];
  int num = atoi(argv[3]);
  int dist = atoi(argv[4]);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, mpe, cen);   // building a system of num equidistant mols

  CMPE::initXForms(mpe);
  CMPE::updateXForms(cen, mpe);
  CMPE::polarize(mpe, false); 
  CMPE::updateXForms(cen, mpe);
  CMPE::polarize(mpe, false); 
  
  vector<CPnt> force, torque;
  vector<REAL> pot(num);
  
  REAL fact = 1.0; //COUL_K/DIELECTRIC_WATER*IKbT;
  
  CMPE::computeForce(mpe, cen, pot, force, torque);
  for (int i = 0; i < num; i++)
    {
      cout << "MOLECULE #" << i+1 << endl;
      cout << "\tENERGY: " << fact*pot[i] << endl;
      force[i] *= fact;
      torque[i] *= fact;
      cout << "\tFORCE: " << force[i].norm() << " " << force[i] << endl;
      cout << "\tTORQUE: " << torque[i].norm() << " " << torque[i] << endl;
    }  

  return 0;
} //end main2

/******************************************************************/
/******************************************************************//**
* Main for perturbation run.
* \param ifname a character pointer holding input file name
* \param num an integer describing the number of iterations to 
						run force calculations
* \param dist an int describing the distance for molecule placement
* \param Dtr a floating point number containing the translational 
					diffusion coefficient of input molecule
* \param Dr a floating point number containing the rotational 
					diffusion coefficient of input molecule
* \param a string containing the perturbation file name 
******************************************************************/
int main3(int argc, char ** argv)
{
	if ( argc != 9 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec per [PQR file] [2, 4, 6 or 8 molecules] [distance between molecules]" <<
							"	[translational diff. coeff] [rotational diff coeff] [runname]" << endl;
		exit(0);
	}
  cout.precision(5);
  cout << "PERTURB" << endl;
  
  const char * ifname = argv[2];
  int num = atoi(argv[3]);
  int dist = atoi(argv[4]);
  REAL Dtr = atof(argv[5]);
  REAL Dr = atof(argv[6]);
	
  seedRand(-1);
	
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, mpe, cen); 
	
  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;
	
  char fn[100];
  sprintf(fn, "perturb_%s_%d_%d.txt", 
					argv[7], num, dist);
  cout << fn << endl;
  ofstream fout(fn);
  CMPE::initXForms(mpe);
  for (int i = 0; i < 200; i++)
	{
		for (int j = 0; j < num; j++)
		{
			CQuat Q = CQuat::chooseRandom();
			mpe[j]->reset(12, Q);
		}
		
		CMPE::solve(mpe,cen,false);
		CMPE::solve(mpe,cen,false);
		perturb(25, num, fout, Dtr, Dr, dist/2.0, mpe, cen);
		cout << i << "..";
		cout.flush();
	}
  cout << endl;
	
  return 0;;
} // end main3

/******************************************************************/
/******************************************************************//**
* Main for Computing the polarization forces, used for comparing the
*			effect of mutual polarization on force and torque.  The output is 
*			a file named polar_[force/torque]_name_nmol_dist.txt
*		The first line indicates the force or torque computed per molecule
*			in the absence of mutual polarization
*		The next line includes mutual polarization.
*   These two lines repeat for 1000 iterations.
* \param ifname is a character pointer for input file name
* \param num is an int of the number of molecules in the system
* \param dist is a floating point number indicating the 
							distance for molecule placement
*	\param a character string to describe the system
******************************************************************/
int main4(int argc, char ** argv)
{
	if ( argc != 6 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec pol [PQR file] [2, 4, 6 or 8 molecules] [distance between molecules] [runname]" << endl;
		exit(0);
	}
  cout.precision(5);
	cout << "POL" << endl;
	
  const char * ifname = argv[2];
  int num = atoi(argv[3]);
  int dist = atoi(argv[4]);
	
  seedRand(-1);
	
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, mpe, cen); 
	
  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;
	
  char fn1[100], fn2[100];
  sprintf(fn1, "polar_force_%s_%d_%d.txt", 
					argv[5], num, dist);
  sprintf(fn2, "polar_torque_%s_%d_%d.txt", 
					argv[5], num, dist);
  ofstream fout1(fn1);
  ofstream fout2(fn2);
	
  CMPE::initXForms(mpe);
  for (int i = 0; i < 1000; i++)
	{
		for (int j = 0; j < num; j++)
		{
			CQuat Q = CQuat::chooseRandom();
			mpe[j]->reset(12, Q);
		}
		
		vector<CPnt> force, torque;
		vector<REAL> pot;
		CMPE::updateXForms(cen, mpe);								// Compute forces and torques ignoring
		CMPE::reexpand(mpe);												// Mutual polarization effects
		CMPE::computeForce(mpe, cen, pot, force, torque);
		for (int j = 0; j < num; j++)
		{
			fout1 << mpe[j]->getOrder() << " " << force[j].x() << " " 
			<< force[j].y() << " " << force[j].z() << " ";
			fout2 << mpe[j]->getOrder() << " " << torque[j].x() << " " 
			<< torque[j].y() << " " << torque[j].z() << " ";
		}
		
		fout1 << endl;
		fout2 << endl;
		
		CMPE::polarize(mpe,false); 
		CMPE::updateXForms(cen, mpe);								// Compute forces and torques with
		CMPE::polarize(mpe,false); 									// Mutual polarization effects included
		CMPE::computeForce(mpe, cen, pot, force, torque);
		
		for (int j = 0; j < num; j++)
		{
			fout1 << mpe[j]->getOrder() << " " << force[j].x() << " " 
			<< force[j].y() << " " << force[j].z() << " ";
			fout2 << mpe[j]->getOrder() << " " << torque[j].x() << " " 
			<< torque[j].y() << " " << torque[j].z() << " ";
		}
		
		fout1 << endl;
		fout2 << endl;
		if (i % 20 == 0)
		{
			cout << ".." << i;
			cout.flush();
		}
	}
	
  cout << endl;
  return 0;;
} // end main4

/******************************************************************/
/******************************************************************//**
* Main for computing the potential on a grid, given a number of num 
* identical molecules placed dist apart in a salt solution.
* \param ifname a character string containing an input file name
* \param num an integer number of molecules to introduce into system
* \param dist a floating point number for distance for molecule placement
******************************************************************/
int main5(int argc, char ** argv)
{
	if ( argc != 5 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec dif [PQR file] [2, 4, 6 or 8 molecules] [distance between molecules]" << endl;
		exit(0);
	}
  cout.precision(5);
  cout << "GDIFF" << endl;
	
  const char * ifname = argv[2];
  int num = atoi(argv[3]);
  REAL dist = atof(argv[4]);
 
  seedRand(num+int(dist)*100);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  REAL rad = buildSystem(ifname, num, dist, true, mpe, cen);  // Build a system of num molecules

  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;

  CMPE::initXForms(mpe);
  if (num > 1)
    {
      CMPE::updateXForms(cen, mpe);
      CMPE::polarize(mpe,false); 
      CMPE::updateXForms(cen, mpe);
      CMPE::polarize(mpe, false); 
    }
 
  REAL fact = COUL_K/DIELECTRIC_WATER*IKbT;
  buildGrid(ifname, num, rad, mpe, cen, fact);

  return 0;
} // end main5

/******************************************************************/
/******************************************************************//**
* Main for computing self-polarization and potential grid for a 
* single coarse-grained molecule positioned at (0,0,0).  Similar to
the option dif, or main 5, but only for a single molecule and  
* \param ifname a character string describing the input file name
* \param fact a floating point scaling factor for the MPE radius
******************************************************************/
int main6(int argc, char ** argv)
{
	if ( argc != 4 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec rad [PQR file] [scaling factor]" << endl;
		exit(0);
	}
  cout.precision(5);
  cout << "RAD" << endl;

  const char * ifname = argv[2];
  REAL fact = atof(argv[3]);
 
  seedRand(-1);
  
  REAL rad1;
  vector<CPnt> p1;
  vector<REAL> ch1;
  CPnt cn;
  readmdp(ifname, p1,ch1,rad1,cn);
  for (int i = 0; i < p1.size(); i++)
    p1[i] -= cn;

  CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, 1, rad1);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  
  CMPE * pm = new CMPE(ch1, p1, fact*rad1, 0, 15);
  mpe.push_back(pm);
  cen.push_back(new CPnt(0.0,0.0,0.0));
  
  char ofname[100], ofname_sp[100];
  int l = strlen(ifname);
  
  strcpy(ofname,ifname);
  strcpy(ofname_sp,ifname);
  sprintf(&(ofname[l-7]), "_1P_%3.1f.mdp", fact);
  sprintf(&(ofname_sp[l-7]), "_1P_%3.1f_sp.mdp", fact);
  
  char rm[200];
  sprintf(rm,"rm -f %s", ofname);
  system(rm);
  sprintf(rm,"rm -f %s", ofname_sp);
  system(rm);
  
  cout << "saving files:" << endl;
  cout << "\t" << ofname << endl;
  cout << "\t" << ofname_sp << endl;
  
  writemdp(ifname, ofname, p1, *(cen[0]), false, ' ');		// print out a PQR file without the CG sphere
  writemdp(ifname, ofname_sp, p1, *(cen[0]), true, ' ');	// print out a PQR file including the CG sphere

  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;

  CMPE::initXForms(mpe);
  buildGrid(ifname, 1, rad1, mpe,cen, 1.0);								// Build a potential grid

  return 0;
} // end of main6

/******************************************************************/
/******************************************************************//**
* Main for computing the effect of a third molecule on the interaction
free energy of two molecules, as a function of separation distance.
The output file generated, [runname]_dist*10.txt has 900 lines, each line
representing the interaction free energy between two molecules.  
The first value in each line represents the system of just two molecules, the
next 30 values represents including a 3rd molecule in the system at 30 different
orientations.
	\param ifname a character string describing input file name
	\param dist a floating point of the distance between two molecules
	\param a string for output for runname
******************************************************************/
int main7(int argc, char ** argv)
{
	if ( argc != 5 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec cng [PQR file] [distance between molecles] [runname]" << endl;
		exit(0);
	}
  cout.precision(5);
  cout << "CHANGE" << endl;
	
  const char * ifname = argv[2];
  REAL dist = atof(argv[3]);
	
  seedRand(-1);
  
  REAL rad1;
  vector<CPnt> p1;
  vector<REAL> ch1;
  CPnt cn;
  readmdp(ifname, p1,ch1,rad1,cn);			// read in the PQR file
  for (int i = 0; i < p1.size(); i++)
    p1[i] -= cn;
	
  CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, 3, rad1);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  REAL fact = rad1 + dist/2.0;						// desired distance between two molecules
	
  for (int n = 0; n < 3; n++)
	{
		CMPE * pm = new CMPE(ch1, p1, rad1, n, 12);
		mpe.push_back(pm);
	}
	
  cen.push_back(new CPnt(fact*CPnt(-1.0,-1.0/sqrt(3.0),0.0)));
  cen.push_back(new CPnt(fact*CPnt(1.0,-1.0/sqrt(3.0),0.0)));
	
  CPnt c3;
  cen.push_back(&c3);										// create 3 molecules in an equilateral triangle
	
  CMPE::initXForms(mpe);
  REAL pi, pj;
	
  char fn[100];
  sprintf(fn,"%s_%dA.txt", argv[4], 
					(int)floor(1.1*dist)); 
  ofstream fout(fn);
  REAL c = 2.0/sqrt(3.0);
  int N = 30;
  REAL f[N], t[N];
  for (int i = 1; i <= N; i++)				// Generate orientations for the spheres to gather data over
	{
		REAL h = -1.0 + 2.0*(i-1.0)/(N-1.0);
		t[i-1] = acos(h);
		if (i == 1 or i == N)
			f[i-1] = 0;
		else 
			f[i-1] = (f[i-2] + 3.6/sqrt(N*(1.0 - h*h)));
		
		while(f[i-1] > 2*M_PI)
			f[i-1] -= 2*M_PI;
		
		while(f[i-1] < -2*M_PI)
			f[i-1] += 2*M_PI;
		cout << t[i-1]*180/M_PI << " " <<  f[i-1]*180/M_PI << endl; // print out generated rotations
	}
	
  for (int i = 0; i < N; i++)
	{
		CQuat Qt(CPnt(0,1,0),t[i]);			// for each orientation generated, rotate first sphere by it
		CQuat Qf(CPnt(0,0,1),f[i]);
		mpe[0]->reset(12, Qf*Qt);
		
		for (int j = 0; j <  N; j++)		// rotate second sphere by it as well 
		{
			Qt = CQuat(CPnt(0,1,0),t[j]);
			Qf = CQuat(CPnt(0,0,1),f[j]);
			mpe[1]->reset(12, Qf*Qt);
			
			mpe[2]->reset(12, CQuat());
			c3 = CPnt(0.0, 1000, 0.0);		// place third sphere very far away so as not to factor into computations
			CMPE::solve(mpe, cen, true);
			CMPE::computePairPot(mpe, 0, 1, pi, pj);  // compute a pairwise potential
			fout << pi  << " ";
			c3 = fact*CPnt(0, c, 0);		// now place it close and at random orientations like the other two
			for (int k = 0; k < N; k++)
	    {
	      Qt = CQuat(CPnt(0,1,0),t[k]);
	      Qf = CQuat(CPnt(0,0,1),f[k]);
	      mpe[2]->reset(12, Qf*Qt);
				
	      CMPE::solve(mpe, cen, true);
	      CMPE::computePairPot(mpe, 0, 1, pi, pj);
	      fout << pi  << " ";	 
	    }
			fout << endl;
		}
		cout << i <<  "..";
		cout.flush();
	}
  
  cout << endl;
  return 0;
} //main7

/******************************************************************/
/******************************************************************//**
* Main for making an infinte grid.  It prints out the potential, the force
and the torque for 8 molecules in a lattice with a given number of lattice
layers
* \param ifname a character string with input file name
* \param layer an integer describing the number of neighbors to 
			consider when performing energy, force and torque calculations
* \param stretch a floating point number describing the scaling of 
			the CG spheres

******************************************************************/
int main8(int argc, char ** argv)
{
	if ( argc != 5 )
	{
		cout << "Correct input format: " << endl;
		cout << " ./exec inf [PQR file] [# layers] [stretch factor]" << endl;
		exit(0);
	}
  cout.precision(5);
  cout << "INFINITE GRID" << endl;
	
  const char * ifname = argv[2];
  int layer = atoi(argv[3]);
  REAL stretch = atof(argv[4]);
	
  seedRand(-1);
  
  int N = 200;
  REAL f[N], t[N];
  for (int i = 1; i <= N; i++)     	// Generate orientations for the spheres to gather data over
	{
		REAL h = -1.0 + 2.0*(i-1.0)/(N-1.0);
		t[i-1] = acos(h);
		if (i == 1 or i == N)
			f[i-1] = 0;
		else 
			f[i-1] = (f[i-2] + 3.6/sqrt(N*(1.0 - h*h)));
		while(f[i-1] > 2*M_PI)
			f[i-1] -= 2*M_PI;
		while(f[i-1] < -2*M_PI)
			f[i-1] += 2*M_PI;
	}
	
  vector<REAL> fs, ts;
  for (int i = 0; i < N; i++)				// Save only orientations that are less than 180 deg
    if ((t[i] <= M_PI/2.0) && (f[i] <= M_PI/2.0))
		{
			fs.push_back(f[i]);
			ts.push_back(t[i]);
		}	
  cout << fs.size() << " orientations" << endl;  // printing them out
	
  char fon[100];
  strcpy(fon, ifname);
  sprintf(&(fon[strlen(fon)-4]), "_88.pdb");		// opening file to write out grid to
	
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  int m = 2*layer+1;
  CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, m*m*m*8, 26.0); // nmol=m^3*8
	
  buildUnitCell(ifname, fon, mpe, cen, true, stretch);
  
  REAL lx = stretch*79.1, ly = stretch*79.1, lz = stretch*37.9;
  
  for (int i = -layer; i <= layer; i++)     // make a grid of layer layers and place points at each location
    for (int j = -layer; j <= layer; j++)
      for (int k = -layer; k <= layer; k++)
			{
				CPnt c(lx*i, ly*j, lz*k);
				for (int n = 0; n < 8; n++)
				{
					if (i == 0 && j == 0 && k == 0)
						continue;
					
					cen.push_back(new CPnt(c + *(cen[n])));
					mpe.push_back(mpe[n]);
				}
			}
	
  int tot = cen.size();
  cout << "N_MOL = " << tot << endl;
	
  CMPE::m_bInfinite = true;
  CMPE::m_unit = 8; 
	
  CMPE::initXForms(mpe);
  REAL k = COUL_K/DIELECTRIC_WATER*IKbT;
  vector<CPnt> force, torque;
  vector<REAL> pot;
  CMPE::solve(mpe,cen,false);
  CMPE::computeForce(mpe, cen, pot, force, torque);
  for (int n = 0; n < 8; n++)
    cout << k*pot[n] << " " << k*force[n] << " " << k*torque[n] << endl;
	
	return 0;
} // end main8

/******************************************************************//**
* Main!
******************************************************************/
int main(int argc, char ** argv)
{
  if (strncmp(argv[1], "sim", 3) == 0)					// For running 2 molecule BD simulation
    return main1(argc,argv);

  else if (strncmp(argv[1], "slv", 3) == 0)			// For computing energies, torques and forces
    return main2(argc,argv);										// of many of the same molecules in solution

  else if (strncmp(argv[1], "per", 3) == 0)			// For computing the energy of many molecules 
    return main3(argc,argv);										// in solution as their rotations and locations are perturbed

  else if (strncmp(argv[1], "pol", 3) == 0)			// For computing the forces/torques of many molecules
    return main4(argc, argv);										// in solution with and without mutual polarization

  else if (strncmp(argv[1], "dif", 3) == 0)			// For computing a grid of potentials due to many molecules
    return main5(argc, argv);										// fixed in solution

  else if (strncmp(argv[1], "rad", 3) == 0)			// For computing a potential grid for a single molecule
    return main6(argc, argv);										// in solution

  else if (strncmp(argv[1], "cng", 3) == 0)			// For computing the effect of a 3rd molecule of the 
    return main7(argc, argv);										// total free energy of a system at distance dist (in A)

  else if (strncmp(argv[1], "inf", 3) == 0)			//
    return main8(argc, argv);										//
  else		// else, bad option
    {
      cout << "bad option!!! " << argv[1] << endl;
      exit(0);
    }
}

