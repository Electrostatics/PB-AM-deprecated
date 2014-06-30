#include <fstream>
#include <vector>

#include "BD.h"
#include "protein.h"
#include "mcoeff.h"

#define b_DIST 100.0
#define q_DIST 500.0
#define f_DIST 100.0

#define DIELECTRIC_WATER 78.0
#define DIELECTRIC_PROT 4.0

#define SALT_CONCENTRATION  0.0500
#define KAPPA (sqrt(SALT_CONCENTRATION)/3.04)

#define COULOMB_CONSTANT (8.988e9)
#define ELECTRON_CHARGE (1.60217733e-19)
#define AVOGADRO_NUM (6.02209e23)
#define KCAL 4186.0 
#define ANGSTROM (1e-10)
#define PICO_SEC (1e-12)
//#define COUL_K 332.058
#define COUL_K 332.141144

#define GAMMA (0.01*100/1000*ANGSTROM*PICO_SEC)
#define Kb (1.380658e-23)
#define TEMP 298.0
#define KbT (1.98718E-03*TEMP) //    (TEMP*Kb*AVOGADRO_NUM/KCAL)
#define IKbT (1.0/KbT)

#define TOL 2.5

#define PATCH_ANGLE 6.0
#define ROTATE_ANGLE 20.0
#define PATCH_SIZE (cos(PATCH_ANGLE*M_PI/180.0))
#define ROTATE_SIZE (cos(ROTATE_ANGLE*M_PI/180.0))


const char CBD::STATUS_STRING[3][10] = {"ESCAPED", "DOCKED", "RUNNING"};
int scount;
char temp_file[100];

CBD::CBD(const char * fname1, const char * fname2, REAL kappa) : m_fcnt(0)
{
  if (!CProtein::isInitialized())
    CProtein::initParameters(kappa, DIELECTRIC_PROT, DIELECTRIC_WATER, 2, 30);

  m_mol1 = new CProtein(fname1); 
  m_mol2 = new CProtein(fname2);
    //saveState();
  m_maxContDist = m_mol1->getRadius() + m_mol2->getRadius() + 1.0;
 
  m_Dtr = 0.030;
  m_Dr1 = 4e-5; m_Dr2 = 4.5e-5;

  computeInterfaceVectors(fname1, fname2);
  /*
  CPnt dir = (m_mol2->getPosition() - m_mol1->getPosition()).normalize();
  dir *= (m_maxContDist-0.3);
  cout << dir << " " << dir.norm() << endl;
  m_mol2->setPos(dir);
 m_mol1->setPos(CPnt());
 CQuat Q(CPnt(-1,1,-1), M_PI/33);
 m_orth2 = Q * m_orth2;
 m_patch2 = Q * m_patch2;

  if (isDocked())
    cout << "docked" << endl;
  else
    cout << "not docked" << endl;

  exit(0);
  */
}

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
  REAL fact = IKbT*COUL_K/DIELECTRIC_WATER;
  REAL srad = m_mol1->getRadius() + m_mol2->getRadius();
  STATUS status = RUNNING;
  scount = 0;
  CPnt dR2, dO1, dO2; 
  while (status == RUNNING)
    {
      REAL dt = compute_dt();     
      REAL dist = CProtein::computeDistance(*m_mol1, *m_mol2);
      if (dist < f_DIST)
	{
	  CProtein::computeForces(force, torque);
	  
	  dR2 = (m_Dtr*dt*fact)*force[1];
	  dO1 = (m_Dr1*dt*IKbT*COUL_K)*torque[0];
	  dO2 = (m_Dr2*dt*IKbT*COUL_K)*torque[1];
	}
      else
	{
	  dR2.zero();
	  dO1.zero();
	  dO2.zero();
	}

      if (scount % 1000 == 0)
	{
	  CPnt t =  m_mol2->getPosition() - m_mol1->getPosition();
	  fout << scount << ") " << m_mol1->getPosition() 
	       << m_mol2->getPosition()	<< "-> " << dist << endl;
	  fout << force[1] << " " << torque[1] << endl;
	  //saveState();
	}

      status = makeMove(dR2, dO1, dO2, dt);
      scount++;
    }

  fout << CBD::STATUS_STRING[status] << endl;
  return status;
}

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

void CBD::saveState()
{
  char fname[100];
  sprintf(fname, "out/output%0.3d.pdb", m_fcnt++);
  ofstream fout(fname);

  m_mol1->PDBoutput(fout, 1, "A", CQuat(), CPnt());

  CQuat rot;
  CPnt trans;
  m_mol1->computeTransformTo(*m_mol2, rot, trans);
  m_mol2->PDBoutput(fout, 1+m_mol1->getAtoms().size(), "B", rot, trans);  
}

void 
CBD::computeInterfaceVectors(const char * fname1, const char * fname2)
{
  CPnt mean;
  vector<CPnt> pnts;
  
  // CProtein::getInterfaceAtoms(*m_mol1, *m_mol2, pnts);

  //  for (int i = 0; i < pnts.size(); i++)
  //  mean += pnts[i];
  
  //mean *= 1.0/pnts.size();
  CPnt d = (m_mol2->getPosition() - m_mol1->getPosition()).normalize();
  m_patch1 = d; //(mean - m_mol1->getPosition()).normalize();
  m_patch2 = -d; //(mean - m_mol2->getPosition()).normalize();   
  m_orth1 = (CPnt(-d.z(), 0.0, d.x())).normalize(); // (cross(m_patch1, m_patch2)).normalize();
  m_orth2 = m_orth1;
  //cout << m_patch2 << m_orth1 << m_orth2 << endl;
}

REAL
CBD::computeRate(int numTraj, int nDocked)
{
  REAL cnst = 4.0*M_PI*m_Dtr; 
  REAL maxR = q_DIST;
  REAL dR = 10.0;
  
  /*
  REAL CONST = m_mol1->getSumCharge()*m_mol2->getSumCharge() * 
    COUL_K*IKbT/DIELECTRIC_WATER;
  
  REAL val = exp(CONST/maxR)/(maxR*maxR);
  while (val > 1e-10)
    {
      maxR += dR;
      val = exp(CONST/maxR)/(maxR*maxR);
    }
  REAL R = b_DIST;
  dR = (maxR - R)*0.000001;
  val = 0.0;
  for (; R <= maxR; R += dR)
    val += exp(CONST/R)*(dR/(R*R));      

  double k_b = cnst/val;

  R = q_DIST;
  dR = (maxR - R)*0.000001;
  val = 0.0;
  for (; R <= maxR; R += dR)
    val += exp(CONST/R)*(dR/(R*R));  
  
  REAL k_q = cnst/val;
  */

  REAL k_b = 4*M_PI*m_Dtr*b_DIST ;
  REAL k_q = 4*M_PI*m_Dtr*q_DIST;
  cout << k_b << " " << k_q << endl;

  REAL delta = ((REAL)nDocked)/numTraj;
  REAL gamma = k_b/k_q;

  REAL K = k_b*delta/(1 - (1 - delta)*gamma);  

  cout << "The relative rate is: " << K << endl;
  return K;
}

void 
approxBBall(const vector<CPnt> & V, CPnt & cen, REAL & rad)
{
  int n = V.size();
  REAL rad2; // radius and radius squared
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
}

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

void readmdp(const char * fname, vector<CPnt> & pnt, vector<REAL> & ch, 
	     REAL & rad, CPnt & cen)
{
  ifstream fin(fname);
  if (!fin.is_open())
    {
      cout << "Could not open file " << fname << endl;
      exit(0);
    }

  char buf[100], temp[10];
  REAL sum = 0.0;

  pnt.clear();
  ch.clear();
  fin.getline(buf,99);
  
  while (!fin.eof())
    {
      double x,y,z,c,r;
      if (strncmp(&(buf[0]),"ATOM",4) == 0)
	{
	  sscanf(&(buf[31]), "%lf %lf %lf", &x, &y, &z);
	  strncpy(temp, &(buf[55]), 6);
	  temp[6] = 0;
	  if (strcmp(temp, "      ") == 0)
	    r = 0.0;
	  else
	    sscanf(temp, "%lf", &r);
	  sscanf(&(buf[61]), "%lf", &c);
	  
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
  readmdp(ifname, V,ch,rad,cn);
  approxBBall(V, cn, rad);

  REAL lx = 79.1, ly = 79.1, lz = 37.9;
  REAL sx = 0.012642, sy = 0.012642, sz = 0.026385;
  CPnt t;

  for (int i = 0; i < V.size(); i++)
    V[i] = CPnt(V[i].x()*sx, V[i].y()*sy, V[i].z()*sz);
  cn = CPnt(cn.x()*sx, cn.y()*sy, cn.z()*sz);
  vector<CPnt> U(V.size());

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
}

REAL buildSystem(const char * ifname, int num, REAL dist, bool bSave, int ver, 
		 vector<CMPE*> & mpe, vector<CPnt*> & cen)
{
  REAL rad1;
  vector<CPnt> p1;
  vector<REAL> ch1;
  CPnt cn;

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

   for (int i = 0; i < p1.size(); i++)
    p1[i] -= cn;

   char ofname[100], ofname_sp[100];
   if (bSave)
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

   CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, num, rad1);
   for (int n = 0; n < num; n++)
    {
      CQuat Q = (num > 2 ? CQuat::chooseRandom() : CQuat());
      vector<CPnt> pr(p1.size());
      for (int i = 0; i < p1.size(); i++)
	{
	  pr[i] = Q*p1[i];
	  if (pr[i].norm() > rad1)
	    cout << "ERROR" << endl;
	}

      CMPE * pm = new CMPE(ch1, pr, rad1, n, 15);
      mpe.push_back(pm);
      if (bSave)
	{	  
	  writemdp(ifname, ofname, pr, *(cen[n]), false, ' ');
	  writemdp(ifname, ofname_sp, pr, *(cen[n]), true, ' ');
	}
    }

   cout << "Building system is complete" << endl;
   return rad1;
}

void buildGrid(const char * ifname, int num, REAL dist, int ver, REAL rad,
	       const vector<CMPE*> & mpe, const vector<CPnt*> & cen, REAL fact)
{
  cout << "building grid" << endl;
  REAL scale = 2.0;
  int gsize = 301;
  int gcen = gsize/2 + 1;
  int bd = gcen-1;
  float * V = new float[gsize*gsize*gsize];

  REAL iscale = 1.0/scale;
  REAL r = rad*rad*scale*scale;

  int gsizeq = gsize*gsize;
  cout << "here" << endl;
  for (int i = -bd; i <= bd ; i++)
    {
      for (int j = -bd; j <= bd ; j++)
	for (int k = -bd; k <= bd ; k++)
	  {
	    bool cont = false;

	    for (int n = 0; n < num; n++)
	      {
		REAL dis_sq = (CPnt(i,j,k) -(scale*(*cen[n]))).normsq();
		if (dis_sq < r)
		  {
		    cont = true;
		    break;
		  }
	      }
	    
	    if (cont)
	      {
		V[(k+bd)*gsizeq+(j+bd)*gsize+(i+bd)] = 0.0;
		continue;
	      }
	    
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
  FILE * fd = fopen(ofname, "w+");
  fwrite(V, sizeof(float), gsize*gsize*gsize, fd);
}

void perturb(int ct, int num, ofstream & fout, REAL Dtr, REAL Dr, REAL dt,
	     vector<CMPE*> & mpe, vector<CPnt*> & cen)
{
  CPnt per[num], rot;
  CQuat Q0[num];
  for (int j = 0; j < num; j++)
    Q0[j] = mpe[j]->getOrient();

  for (int k = 0; k < ct; k++)
    {
      // cout << "------ " << k << " -------" << endl;
      //for (int i = 0; i < cen.size(); i++)
      //	cout << *(cen[i]);
      //cout << endl;
      for (int j = 0; j < num; j++)
	{
	  per[j] = CBD::getRandVec(sqrt(2*dt*Dtr));
	  *(cen[j]) += per[j];

	  rot = CBD::getRandVec(sqrt(2*dt*Dr));
	  CQuat Q(rot, rot.norm());
	  mpe[j]->setOrient(Q*Q0[j]);
  	}

      //for (int i = 0; i < cen.size(); i++)
      //	cout << *(cen[i]);
      //cout << endl;

      CMPE::updateSolve(mpe, cen);
      //      CMPE::computeForce(mpe, cen, force, torque);

      fout << CMPE::npol << " " << CMPE::npol_t << " ";
      for (int i = 0; i < num; i++)
	fout << mpe[i]->ngpol << " " << mpe[i]->ngpol_t << " ";
      fout << endl;

      CMPE::undoXForms();
      for (int i = 0; i < num; i++)
	{
	  mpe[i]->undo();
	  *(cen[i]) -= per[i];
	}
    }
}

// param 2: salt concentration
// param 3: output file
// param 4: temporary file index
int main1(int argc, char ** argv)
{
  seedRand(-1);

  REAL kappa = sqrt(atof(argv[2]))/3.04;
  CBD bd("barnaseh.pdb", "barstarh.pdb",kappa);
  
  ofstream fout(argv[3]);
  sprintf(temp_file, "temp%d.out", atoi(argv[4]));

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
}

int main2(int argc, char ** argv)
{
  cout << "SOLVE" << endl;
  seedRand(-1);
  cout.precision(5);
 
  int num = atoi(argv[3]);
  const char * ifname = argv[2];
  int dist = atoi(argv[4]);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, 0, mpe, cen); 

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
}

int main3(int argc, char ** argv)
{
  cout.precision(5);
  cout << "PERTURB" << endl;

  int num = atoi(argv[3]);
  const char * ifname = argv[2];
  int dist = atoi(argv[4]);
  int ver = atoi(argv[5]);
  REAL Dtr = atof(argv[6]);
  REAL Dr = atof(argv[7]);

  //seedRand(dist+num*100+ver*10000);
  seedRand(-1);

  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, ver, mpe, cen); 

  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;

  char fn[100];
  sprintf(fn, "../../Documents/JCTC/res/perturb/perturb_%s_%d_%d.txt", 
	  argv[8], num, dist);
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

  //buildGrid(ifname, num, dist, ver, mpe[0]->getRad(), mpe, cen);

  return 0;;
}

int main4(int argc, char ** argv)
{
  cout.precision(5);

  int num = atoi(argv[3]);
  const char * ifname = argv[2];
  int dist = atoi(argv[4]);
  int ver = atoi(argv[5]);

  //seedRand(dist+num*100+ver*10000);
  seedRand(-1);

  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  buildSystem(ifname, num, dist, false, ver, mpe, cen); 

  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;

  char fn1[100], fn2[100];
  sprintf(fn1, "../../Documents/JCTC/res/polar/polar_force_%s_%d_%d.txt", 
	  argv[6], num, dist);
  sprintf(fn2, "../../Documents/JCTC/res/polar/polar_torque_%s_%d_%d.txt", 
	  argv[6], num, dist);
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
      CMPE::updateXForms(cen, mpe);
      CMPE::reexpand(mpe);
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
      CMPE::updateXForms(cen, mpe);
      CMPE::polarize(mpe,false); 
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
}

int main5(int argc, char ** argv)
{
  cout.precision(5);
  cout << "GDIFF" << endl;

  int num = atoi(argv[3]);
  const char * ifname = argv[2];
  REAL dist = atof(argv[4]);
  int ver = atoi(argv[5]);
 
  seedRand(num+ver*100);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  REAL rad = buildSystem(ifname, num, dist, true, ver, mpe, cen); 

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
  buildGrid(ifname, num, dist, ver, rad, mpe,cen,fact);

  return 0;
}

int main6(int argc, char ** argv)
{
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
  sprintf(&(ofname[l-7]), "_1P.mdp");
  sprintf(&(ofname_sp[l-7]), "_1P_sp.mdp");
  
  char rm[200];
  sprintf(rm,"rm -f %s", ofname);
  system(rm);
  sprintf(rm,"rm -f %s", ofname_sp);
  system(rm);
  
  cout << "saving files:" << endl;
  cout << "\t" << ofname << endl;
  cout << "\t" << ofname_sp << endl;
  
  writemdp(ifname, ofname, p1, *(cen[0]), false, ' ');
  writemdp(ifname, ofname_sp, p1, *(cen[0]), true, ' ');

  for (int i = 0; i < cen.size(); i++)
    cout << *(cen[i]);
  cout << endl;

  CMPE::initXForms(mpe);
  int f = (int) (fact*100);
  buildGrid(ifname, 1, f, 0, rad1, mpe,cen,1.0);

  return 0;
}

int main7(int argc, char ** argv)
{
  cout.precision(5);
  cout << "CHANGE" << endl;

  const char * ifname = argv[2];
  REAL dist = atof(argv[3]);
 
  seedRand(-1);
  
  REAL rad1;
  vector<CPnt> p1;
  vector<REAL> ch1;
  CPnt cn;
  readmdp(ifname, p1,ch1,rad1,cn);
  for (int i = 0; i < p1.size(); i++)
    p1[i] -= cn;

  CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, 3, rad1);
  
  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  REAL fact = rad1 + dist/2.0;

  for (int n = 0; n < 3; n++)
    {
      CMPE * pm = new CMPE(ch1, p1, rad1, n, 12);
      mpe.push_back(pm);
    }

  cen.push_back(new CPnt(fact*CPnt(-1.0,-1.0/sqrt(3.0),0.0)));
  cen.push_back(new CPnt(fact*CPnt(1.0,0.-1.0/sqrt(3.0),0.0)));
 
  CPnt c3;
  cen.push_back(&c3);

  CMPE::initXForms(mpe);
  REAL pi, pj;

  char fn[100];
  sprintf(fn,"../../Documents/JCTC/res/change/%s_%dA.txt", argv[4], 
	  (int)floor(10*dist)); 
  ofstream fout(fn);
  REAL c = 2.0/sqrt(3.0);
  int N = 30;
  REAL f[N], t[N];
  for (int i = 1; i <= N; i++)      
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
      cout << t[i-1]*180/M_PI << " " <<  f[i-1]*180/M_PI << endl; 
    }


  for (int i = 0; i < N; i++)
    {
      CQuat Qt(CPnt(0,1,0),t[i]);
      CQuat Qf(CPnt(0,0,1),f[i]);
      mpe[0]->reset(12, Qf*Qt);

      for (int j = 0; j <  N; j++)
	{
	  Qt = CQuat(CPnt(0,1,0),t[j]);
	  Qf = CQuat(CPnt(0,0,1),f[j]);
	  mpe[1]->reset(12, Qf*Qt);

	  mpe[2]->reset(12, CQuat());
	  c3 = CPnt(0.0, 1000, 0.0);
	  CMPE::solve(mpe, cen, true);
	  CMPE::computePairPot(mpe, 0, 1, pi, pj);
	  fout << pi  << " ";
	  c3 = fact*CPnt(0, c, 0);
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
}

int main8(int argc, char ** argv)
{
  cout.precision(5);
  cout << "INFINITE GRID" << endl;

  const char * ifname = argv[2];
  int dist = atoi(argv[3]);
  int layer = atoi(argv[4]);
  REAL stretch = atof(argv[5]);

  seedRand(-1);
  
  int N = 200;
  REAL f[N], t[N];
  for (int i = 1; i <= N; i++)      
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
  for (int i = 0; i < N; i++)
    if ((t[i] <= M_PI/2.0) && (f[i] <= M_PI/2.0))
      {
	fs.push_back(f[i]);
	ts.push_back(t[i]);
      }	
  cout << fs.size() << " orientations" << endl;

  char fon[100];
  strcpy(fon, ifname);
  sprintf(&(fon[strlen(fon)-4]), "_88.pdb");

  vector<CMPE*> mpe;
  vector<CPnt*> cen;
  int m = 2*layer+1;
  CProtein::initParameters(KAPPA, DIELECTRIC_PROT, DIELECTRIC_WATER, m*m*m*8, 26.0);

  buildUnitCell(ifname, fon, mpe, cen, true, stretch);
  
  REAL lx = stretch*79.1, ly = stretch*79.1, lz = stretch*37.9;
  
  for (int i = -layer; i <= layer; i++)
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

  /*
  char fn[100];
  sprintf(fn,"../../Documents/JCTC/res/infgrid/%s_%dA_%dL.txt", 
	  argv[5], dist, layer);
  ofstream fout(fn);
 
  
  for (int i = 0; i < fs.size(); i++)
    {
      CQuat Qt(CPnt(0,1,0),ts[i]);
      CQuat Qf(CPnt(0,0,1),fs[i]);
      pm->reset(12, Qf*Qt);

      CMPE::solve(mpe,cen,false);
      CMPE::computeForce(mpe, cen, pot, force, torque);

      fout << k*pot[0] << " " << k*force[0].x() << " " << k*force[0].y() << " "
	 << k*force[0].z() << " " << k*torque[0].x() << " " << k*torque[0].y() 
	 << " " << k*torque[0].z() << endl;
    }
  */
 return 0;
}

int main(int argc, char ** argv)
{
  if (strncmp(argv[1], "sim", 3) == 0)
    return main1(argc,argv);
  else if (strncmp(argv[1], "slv", 3) == 0)
    return main2(argc,argv);
  else if (strncmp(argv[1], "per", 3) == 0)
    return main3(argc,argv);
  else if (strncmp(argv[1], "pol", 3) == 0)
    return main4(argc, argv);
  else if (strncmp(argv[1], "dif", 3) == 0)
    return main5(argc, argv);
  else if (strncmp(argv[1], "rad", 3) == 0)
    return main6(argc, argv);
  else if (strncmp(argv[1], "cng", 3) == 0)
    return main7(argc, argv);
  else if (strncmp(argv[1], "inf", 3) == 0)
    return main8(argc, argv);
  else
    {
      cout << "bad option!!! " << argv[1] << endl;
      exit(0);
    }
}

