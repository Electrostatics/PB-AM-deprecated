#ifndef _UTIL_H_
#define _UTIL_H_

#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <memory>

using namespace std;

typedef double REAL;
typedef complex<REAL> Complex;

class CPnt;
class CQuat;

class CSpPnt
{
 public:
  CSpPnt() 
    { zero(); }
  CSpPnt(REAL rho_, REAL theta_, REAL phi_) :
    m_rho(rho_), m_theta(theta_), m_phi(phi_) {}
  const CSpPnt & operator=(const CSpPnt & s) 
    { m_rho = s.rho(); m_theta = s.theta(); m_phi = s.phi(); return *this; }

  void zero()
    { m_rho = 0.0; m_theta = 0.0; m_phi = 0.0; }

  REAL & rho() 
  { return m_rho; }
  REAL & theta() 
  { return m_theta; }
  REAL & phi() 
  {return m_phi; }

  const REAL rho() const
    { return m_rho; }
  const REAL theta() const
  { return m_theta; }
  const REAL phi() const
    {return m_phi; }

  friend ostream & operator<<(ostream & out, const CSpPnt & p);
  friend CPnt SphToCart(const CSpPnt & s);
  friend CSpPnt CartToSph(const CPnt & c);

 private:
  REAL m_rho, m_theta, m_phi;
};

REAL dot(const CPnt & p1, const CPnt & p2);

class CPnt
{
 public:
  CPnt() 
    { zero(); }
  CPnt(REAL x_, REAL y_, REAL z_) :
    m_x(x_), m_y(y_), m_z(z_) {}
  CPnt(const CPnt & p) : 
    m_x(p.x()), m_y(p.y()), m_z(p.z()) {}
  const CPnt & operator=(const CPnt & c)
    { m_x = c.x(); m_y = c.y(); m_z = c.z(); return *this; }
  void zero()
    { m_x = 0.0; m_y = 0.0; m_z = 0.0; }
  
  REAL normsq() const
    { return dot(*this,*this); }
  REAL norm() const
    { return sqrt(normsq()); }
  const CPnt & normalize()
    { (*this) *= (1.0/norm()); return *this; }
  CPnt normalize() const
  { CPnt P = *this; P *= (1.0/norm()); return P; }
  void operator+=(const CPnt & p)
  { m_x += p.x(); m_y += p.y(); m_z += p.z(); }
  void operator-=(const CPnt & p)
   { m_x -= p.x(); m_y -= p.y(); m_z -= p.z(); }
  void operator*=(REAL s)
    { m_x *= s; m_y *= s; m_z *= s; }
  void operator/=(REAL s)
  { REAL t = 1.0/s; *this *= t; }

  const REAL x() const
  { return m_x; }
  const REAL y() const
  { return m_y; }
  const REAL z() const
  { return m_z; }
  
  REAL & x() 
  { return m_x; }
  REAL & y() 
  { return m_y; }
  REAL & z() 
  { return m_z; }

  friend CPnt operator+(const CPnt & p1, const CPnt & p2);
  friend CPnt operator-(const CPnt & p1, const CPnt & p2);

  friend CPnt cross(const CPnt & p1, const CPnt & p2);
  friend  CPnt operator*(REAL c, const CPnt & p1);
  friend CPnt operator-(const CPnt & p);

  friend ostream & operator<<(ostream & out, const CPnt & p);
  friend CPnt SphToCart(const CSpPnt & s);
  friend CSpPnt CartToSph(const CPnt & c);
  friend REAL torsion(const CPnt & p1, const CPnt & p2,
		      const CPnt & p3, const CPnt & p4);
  friend REAL torsion(const CPnt & v1, const CPnt & v2,
		      const CPnt & v3);

  
 private:
  REAL m_x, m_y, m_z;
};

class CQuat
{
 public:
  static CQuat chooseRandom();
  CQuat()
    { identity(); }
  CQuat(REAL real, const CPnt & imag) :
    m_real(real), m_imag(imag) { normalize(); }
  CQuat(const CPnt & axis, REAL angle) :
    m_real(cos(0.5*angle)), m_imag(sin(0.5*angle)*(axis.normalize()))
    { normalize(); }
  CQuat(const CQuat & q)
    { *this = q; }
  CQuat & operator=(const CQuat & q)
    { m_imag = q.m_imag; m_real = q.m_real; return *this; }
  void identity()
    { m_imag.zero(); m_real = 1.0;}
  
  CQuat & operator*=(const CQuat & q);
  void normalize()
    { 
      REAL in_norm = 1.0/sqrt(m_imag.normsq() + m_real*m_real);
      m_imag *= in_norm; m_real *= in_norm;
    }  
  void conj()
    { m_imag = -m_imag; }

  REAL real() const
    { return m_real; }
  const CPnt & imag() const
    { return m_imag; }

  REAL x() const
  { return m_imag.x(); }
  REAL y() const
  { return m_imag.y(); }
  REAL z() const
  { return m_imag.z(); }
  REAL w() const
  { return m_real; }
  
  friend CQuat operator*(const CQuat & q1, const CQuat & q2);
  friend CPnt operator*(const CQuat & q1, const CPnt & p1);
  friend CQuat conj(const CQuat & q);
  friend ostream & operator<<(ostream & out, const CQuat & q);

 private:
  REAL m_real;
  CPnt m_imag;
};

CPnt randOrient();


inline CPnt 
operator-(const CPnt & c1, const CPnt & c2)
{
  CPnt c(c1.x()-c2.x(), c1.y()-c2.y(), c1.z()-c2.z());
  return c;
}

inline CPnt
operator+(const CPnt & c1, const CPnt & c2)
{
  CPnt c(c1.x()+c2.x(), c1.y()+c2.y(), c1.z()+c2.z());
  return c;
}

inline CPnt 
operator*(REAL c, const CPnt & p1)
{
  CPnt p(c * p1.x(), c * p1.y(), c * p1.z());
  return p;
}

inline REAL
dot(const CPnt & p1, const CPnt & p2)
{
  return (p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z());
}

inline CPnt 
cross(const CPnt & p1, const CPnt & p2)
{
  CPnt p(p1.y()*p2.z() - p1.z()*p2.y(),
	 p1.z()*p2.x() - p1.x()*p2.z(),
	 p1.x()*p2.y() - p1.y()*p2.x());
  
  return p;
}

inline CPnt 
operator-(const CPnt & p1)
{
  CPnt p(-p1.x(), -p1.y(), -p1.z()); 

  return p;
}

inline REAL
torsion(const CPnt & p1, const CPnt & p2, const CPnt & p3, const CPnt & p4)
{
  return torsion(p1 - p2, p3 - p2, p4 - p3);
}

inline REAL
torsion(const CPnt & v1, const CPnt & v2, const CPnt & v3)
{
  //      U     W      V
  //   o<----o----->o----->o
  //

  CPnt A = cross(v1,v2);
  CPnt B = cross(v3,v2);

  REAL F = dot(cross(A,B),v2);
  REAL T = v2.norm() * dot(A,B);

  return atan2(F,T);
}

inline ostream & 
operator<<(ostream & out, const CPnt & p)
{
  out << "[" << p.x() << "," << p.y() << "," << p.z() <<"]";
  return out;
}

inline ostream & 
operator<<(ostream & out, const CSpPnt & p)
{
  out << "[" << p.rho() << "," << p.theta() << "," << p.phi() << "]";
  return out;
}

inline CPnt 
SphToCart(const CSpPnt & s)
{
  REAL x = s.rho()*sin(s.theta())*cos(s.phi());
  REAL y = s.rho()*sin(s.theta())*sin(s.phi());
  REAL z = s.rho()*cos(s.theta());
  
  return CPnt(x,y,z);
}

inline CSpPnt
CartToSph(const CPnt & c)
{
  REAL theta, phi;
  REAL rho = sqrt(c.x()*c.x() + c.y()*c.y() + c.z()*c.z());

  if (rho < fabs(c.z())) 
    rho = fabs(c.z());

  if (rho == 0.0) 
    theta = 0.0;
  else 
    theta = acos(c.z()/rho);

  if ((c.x() == 0.0) && (c.y() == 0.0)) 
    phi = 0.0;
  else 
    phi = atan2(c.y(), c.x());

  return CSpPnt(rho, theta, phi);
}

inline CQuat &
CQuat::operator*=(const CQuat & q)
{
  REAL temp = q.m_real*m_real - dot(q.m_imag, m_imag);
  m_imag =  q.m_real*m_imag + m_real*q.m_imag +	cross(q.m_imag, m_imag);
  m_real = temp;

  return *this;
}

inline CQuat 
operator*(const CQuat & q1, const CQuat & q2)
{
  CQuat q = q2;
  q *= q1;
  
  return q;
}

inline CPnt 
operator*(const CQuat & q, const CPnt & r)
{
  CPnt qr = cross(q.m_imag, r);
  return (r + 2*(q.m_real*qr + cross(q.m_imag, qr)));
}

inline CQuat 
conj(const CQuat & q)
{
  return CQuat(q.m_real, -q.m_imag);
}


inline ostream & 
operator<<(ostream & out, const CQuat & q)
{
  out << "{" << q.m_real << " " << q.m_imag << "}";
  return out;
}

inline CQuat
CQuat::chooseRandom()
{
  REAL s = drand48();
  REAL sig1 = sqrt(1.0-s);
  REAL sig2 = sqrt(s);
  REAL t1 = 2.0*drand48()*M_PI;
  REAL t2 = 2.0*drand48()*M_PI;
  CQuat Q(cos(t2)*sig2, CPnt(sin(t1)*sig1, cos(t1)*sig1, sin(t2)*sig2));

  return Q;
}

inline CPnt
randOrient()
{
  REAL phi = drand48()*2*M_PI;
  REAL u = drand48()*2 - 1;

  return CPnt(sqrt(1 - u*u) * cos(phi), sqrt(1 - u*u) * sin(phi), u);
}

inline REAL
normRand()
{
  static double saved;
  static bool bSaved = false;

  double v1, v2, rsq;
  if (!bSaved)
    {
      do {
        v1 = 2.0 * drand48() - 1.0;
        v2 = 2.0 * drand48() - 1.0;
        rsq = v1*v1 + v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);

      double fac = sqrt(-2.0*log(rsq)/rsq);
      saved = v1*fac;
      bSaved = true;
      return (REAL) fac*v2;
    }
  else
    {
      bSaved = false;
      return (REAL) saved;
    }
}

inline void
seedRand(int s)
{
  if (s == -1)
    {
      struct timeval tvsd;
      gettimeofday(&tvsd, NULL);
      srand48(tvsd.tv_sec);
    }
  else
    srand48(s);
}


#endif
