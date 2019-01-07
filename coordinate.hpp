//
//    AttenuationVoronoi : A software used in a study of the attenuation of the Earths
//    inner core, see 
//
//      Pejic T., Hawkins R., Sambridge M. & Tkalcic H. "Trans-dimensional Bayesian attenuation tomography 
//      of the upper inner core", Journal of Geophysical Research, 2019, to appear.
//    
//    Copyright (C) 2014 - 2018 Rhys Hawkins
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#pragma once
#ifndef coordinate_hpp
#define coordinate_hpp

#include <cmath>

#include "attenuationutil.hpp"
#include "attenuationexception.hpp"

//
// Forward
//
template
<typename real>
class sphericalcoordinate;

template
<typename real>
class cartesiancoordinate;

template
<typename real>
class vector3;

template
<
  typename real,
  typename coord
>
class barycentrecoordinate {
public:

  static constexpr real BARYCENTRE_EPSILON = 1.0e-9;

  barycentrecoordinate() :
    alpha(0.0),
    beta(0.0),
    gamma(0.0)
  {
  }
  
  barycentrecoordinate(real _alpha, real _beta, real _gamma) :
    alpha(_alpha),
    beta(_beta),
    gamma(_gamma)
  {
  }

  bool inside() const {
    if (alpha < -BARYCENTRE_EPSILON || alpha > 1.0) {
      return false;
    }
    
    if (beta < -BARYCENTRE_EPSILON || beta > 1.0) {
      return false;
    }
    
    if (gamma < -BARYCENTRE_EPSILON || gamma > 1.0) {
      return false;
    }
    
    return true;
  }

  real interpolate(real a, real b, real c) const
  {
    return (alpha * a + beta * b + gamma * c);
  }
  
  real alpha, beta, gamma;
};

template
<
  typename real
>
class barycentrecoordinate<real, sphericalcoordinate<real>> {
public:
  
  static constexpr real BARYCENTRE_EPSILON = 1.0e-9;

  barycentrecoordinate() :
    alpha(0.0),
    beta(0.0),
    gamma(0.0)
  {
  }
  
  barycentrecoordinate(real _alpha, real _beta, real _gamma) :
    alpha(_alpha),
    beta(_beta),
    gamma(_gamma)
  {
  }

  bool inside() const {
    if (alpha < -BARYCENTRE_EPSILON || alpha > 1.0) {
      return false;
    }
    
    if (beta < -BARYCENTRE_EPSILON || beta > 1.0) {
      return false;
    }
    
    if (gamma < -BARYCENTRE_EPSILON || gamma > 1.0) {
      return false;
    }
    
    return true;
  }

  real interpolate(real a, real b, real c) const
  {
    return (alpha * a + beta * b + gamma * c)/(alpha + beta + gamma);
  }
  
  real alpha, beta, gamma;
};
  
template
<typename real>
class vector3 {
public:

  vector3() :
    x(0.0),
    y(0.0),
    z(0.0)
  {
  }

  vector3(real _x, real _y, real _z) :
    x(_x),
    y(_y),
    z(_z)
  {
  }

  real length() const
  {
    return sqrt(x*x + y*y + z*z);
  }

  real dot(const vector3 &v) const
  {
    return x*v.x + y*v.y + z*v.z;
  }

  vector3 &operator+=(const vector3 &v)
  {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
  }
  
  vector3 &operator-=(const vector3 &v)
  {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
  }

  vector3 &operator*=(const real &s)
  {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }

  vector3 &operator/=(const real &s)
  {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  friend vector3 operator+(const vector3 &a, const vector3 &b)
  {
    vector3 s(a);
    s += b;
    return s;
  }
  
  friend vector3 operator-(const vector3 &a, const vector3 &b)
  {
    vector3 s(a);
    s -= b;
    return s;
  }

  friend vector3 operator*(const vector3 &a, const real &s)
  {
    vector3 r(a);
    r *= s;
    return r;
  }

  friend vector3 operator*(const real &s, const vector3 &a)
  {
    vector3 r(a);
    r *= s;
    return r;
  }
  
  friend vector3 operator/(const vector3 &a, const real &s)
  {
    vector3 r(a);
    r /= s;
    return r;
  }

  friend vector3 cross(const vector3 &a, const vector3 &b)
  {
    return vector3(a.y*b.z - a.z*b.y,
		   a.z*b.x - a.x*b.z,
		   a.x*b.y - a.y*b.x);

  }

  friend bool computeaxisangle(const vector3 &a,
			       const vector3 &b,
			       vector3 &axis,
			       real &angle)
  {
    axis = cross(a, b);
    real l = axis.length();
    if (l == 0.0) {
      return false;
    }

    axis /= l;
    angle = acos(a.dot(b)/(a.length() * b.length()));

    return true;
  }

  friend vector3 rotate(const vector3 &a, const vector3 &axis, const real &angle)
  {
    return cos(angle) * a + sin(angle) * cross(axis, a) + (1.0 - cos(angle))*(axis.dot(a))*axis;
  }
			       
  static real determinant(const vector3 &a,
			  const vector3 &b,
			  const vector3 &c)
  {
    return
      (a.x * b.y * c.z) +
      (b.x * c.y * a.z) +
      (c.x * a.y * b.z) -
      (c.x * b.y * a.z) -
      (b.x * a.y * c.z) -
      (a.x * c.y * b.z);
  }

  real x, y, z;
};

template
<typename real>
class plane {
public:

  plane(const vector3<real> &p1,
	const vector3<real> &p2,
	const vector3<real> &p3)
  {
    normal = cross(p2 - p1, p3 - p1);
    real l = normal.length();
    if (l == 0.0) {
      throw ATTENUATIONEXCEPTION("Colinear points");
    }

    normal /= l;

    //
    // Hessian normal form
    //
    d = -normal.dot(p1);
  }

  real distance(const vector3<real> &p) const
  {
    return p.dot(normal) + d;
  }

  vector3<real> normal;
  real d;
};
  

template
<typename real>
class cartesiancoordinate {
public:

  cartesiancoordinate() :
    x(0.0),
    y(0.0)
  {
  }
  
  cartesiancoordinate(double _x, double _y) :
    x(_x),
    y(_y)
  {
  }

  //
  // Construct as midpoint
  //
  cartesiancoordinate(const cartesiancoordinate &a, const cartesiancoordinate &b) :
    x((a.x + b.x)/2.0),
    y((a.y + b.y)/2.0)
  {
  }

  //
  // Construct as linearly interpolated point alpha in -1 .. 1, alpha = 0 is equivalent to midpoint
  //
  cartesiancoordinate(const cartesiancoordinate &a, const cartesiancoordinate &b, real alpha) :
    x(((1.0 - alpha)*a.x + (1.0 + alpha)*b.x)/2.0),
    y(((1.0 - alpha)*a.y + (1.0 + alpha)*b.y)/2.0)
  {
  }

  cartesiancoordinate(const barycentrecoordinate<real, cartesiancoordinate> &bc,
		      const cartesiancoordinate &a,
		      const cartesiancoordinate &b,
		      const cartesiancoordinate &c) :
    x(a.x * bc.alpha + b.x * bc.beta + c.x * bc.gamma),
    y(a.y * bc.alpha + b.y * bc.beta + c.y * bc.gamma)
  {
  }
  
  barycentrecoordinate<real, cartesiancoordinate> get_barycentre_coordinate(const cartesiancoordinate &a,
									    const cartesiancoordinate &b,
									    const cartesiancoordinate &c) const
  {
    real detT =
      ((a.x - c.x) * (b.y - c.y)) -
      ((b.x - c.x) * (a.y - c.y));

    if (detT == 0.0) {
      throw ATTENUATIONEXCEPTION("Linearly dependent points\n");
    }

    real dx = x - c.x;
    real dy = y - c.y;

    barycentrecoordinate<real, cartesiancoordinate> bc;

    bc.alpha = ((b.y - c.y)*dx + (c.x - b.x)*dy)/detT;
    bc.beta = ((c.y - a.y)*dx + (a.x - c.x)*dy)/detT;
    bc.gamma = 1.0 - (bc.alpha + bc.beta);
    
    return bc;
  }

  real get_interpolation_alpha(const cartesiancoordinate &a,
			       const cartesiancoordinate &b) const
  {
    //
    // It is assumed that this point lies on the line from a to b
    //
    if (fabs(a.x - b.x) > fabs(a.y - b.y)) {

      return 2.0*(x - a.x)/(b.x - a.x) - 1.0;
      
    } else {

      return 2.0*(y - a.y)/(b.y - a.y) - 1.0;
      
    }
  }

  bool writetext(FILE *fp) const
  {
    double dx = scalartodouble<real>(x);
    double dy = scalartodouble<real>(y);

    fprintf(fp, "%.15g %.15g", dx, dy);
    return true;
  }

  bool readtext(FILE *fp)
  {
    double dx, dy;
    if (fscanf(fp, "%lf %lf", &dx, &dy) != 2) {
      return false;
    }

    x = doubletoscalar<real>(dx);
    y = doubletoscalar<real>(dy);
    return true;
  }
  bool writebinary(FILE *fp) const
  {
    return value_write<real>(fp, x) && value_write<real>(fp, y);
  }

  bool readbinary(FILE *fp)
  {
    return value_read<real>(fp, x) && value_read<real>(fp, y);
  }

  bool write(std::ostream &s)
  {
    s.write((char*)&x, sizeof(real));
    s.write((char*)&y, sizeof(real));

    return true;
  }
  
  bool read(std::istream &s)
  {
    s.read((char*)&x, sizeof(real));
    s.read((char*)&y, sizeof(real));

    return true;
  }

  friend bool operator==(const cartesiancoordinate &a,
			 const cartesiancoordinate &b)
  {
    return (a.x == b.x) && (a.y == b.y);
  }

  friend bool operator!=(const cartesiancoordinate &a,
			 const cartesiancoordinate &b)
  {
    return (a.x != b.x) || (a.y != b.y);
  }

  real x;
  real y;

};

template
<typename real>
class sphericalcoordinate {
public:

  sphericalcoordinate() :
    phi(0.0),
    theta(0.0)
  {
  }
  
  sphericalcoordinate(real _phi, real _theta) :
    phi(_phi),
    theta(_theta)
  {
  }

  sphericalcoordinate(real x, real y, real z)
  {
    real r;
    cartesiantospherical(vector3<real>(x, y, z),
			 phi,
			 theta,
			 r);
  }
  
  sphericalcoordinate(const vector3<real> &v)
  {
    real r;
    cartesiantospherical(v, phi, theta, r);
  }

  sphericalcoordinate(const sphericalcoordinate &a, const sphericalcoordinate &b)
  {
    vector3<real> v1, v2, vc;
    real r;

    //
    // This assumes that the two points are not antipodal
    //
    sphericaltocartesian(a.phi, a.theta, v1);
    sphericaltocartesian(b.phi, b.theta, v2);

    vector3<real> axis;
    plane<real> p(v1, vector3<real>(0.0, 0.0, 0.0), v2);
    
    real angle;

    if (!computeaxisangle(v1, v2, axis, angle)) {
      throw ATTENUATIONEXCEPTION("Coincident points\n");
    }

    vc = rotate(v1, axis, angle/2.0);

    //
    // Correction to ensure point is on correct plane 
    //
    // vc -= p.distance(vc) * p.normal;

    cartesiantospherical(vc,
			 phi, theta, r);
    
  }

  sphericalcoordinate(const barycentrecoordinate<real, sphericalcoordinate> &bc,
		      const sphericalcoordinate &a,
		      const sphericalcoordinate &b,
		      const sphericalcoordinate &c)
  {
    vector3<real> va, vb, vc;

    sphericaltocartesian(a.phi, a.theta, va);
    sphericaltocartesian(b.phi, b.theta, vb);
    sphericaltocartesian(c.phi, c.theta, vc);

    vector3<real> v = bc.alpha*va + bc.beta*vb + bc.gamma*vc;
    real r;
    cartesiantospherical(v, phi, theta, r);
  }

  sphericalcoordinate(const sphericalcoordinate &a,
		      const sphericalcoordinate &b,
		      real alpha)
  {
    vector3<real> v1, v2, vc;
    real r;

    //
    // This assumes that the two points are not antipodal
    //
    sphericaltocartesian(a.phi, a.theta, v1);
    sphericaltocartesian(b.phi, b.theta, v2);

    plane<real> p(v1, vector3<real>(0.0, 0.0, 0.0), v2);
    vector3<real> axis;
    real angle;

    if (!computeaxisangle(v1, v2, axis, angle)) {
      throw ATTENUATIONEXCEPTION("Coincident points\n");
    }

    vc = rotate(v1, axis, (alpha + 1.0)/2.0 * angle);

    //
    // Correction to ensure point is on correct plane 
    //
    // vc -= p.distance(vc) * p.normal;

    cartesiantospherical(vc,
			 phi, theta, r);
  }
  
    
  barycentrecoordinate<real, sphericalcoordinate> get_barycentre_coordinate(const sphericalcoordinate &a,
									    const sphericalcoordinate &b,
									    const sphericalcoordinate &c) const
  {
    //
    // This algorithm is from Alfeld et al, "Bernstein-Bezier polynomials on sphere and sphere-like surfaces",
    // Computer Aided Geometric Design, 1996
    //

    vector3<real> p;
    vector3<real> va, vb, vc;

    sphericaltocartesian(a.phi, a.theta, va);
    sphericaltocartesian(b.phi, b.theta, vb);
    sphericaltocartesian(c.phi, c.theta, vc);

    sphericaltocartesian(phi, theta, p);
    
    real detT = vector3<real>::determinant(va, vb, vc);

    if (detT == 0.0) {
      throw ATTENUATIONEXCEPTION("Linearly dependent points (%f %f : %f %f %f) (%f %f : %f %f %f) (%f %f : %f %f %f)\n",
				 a.phi, a.theta, va.x, va.y, va.z,
				 b.phi, b.theta, vb.x, vb.y, vb.z,
				 c.phi, c.theta, vc.x, vc.y, vc.z);
    }

    barycentrecoordinate<real, sphericalcoordinate> bc;

    bc.alpha = vector3<real>::determinant(p, vb, vc)/detT;
    bc.beta = vector3<real>::determinant(va, p, vc)/detT;
    bc.gamma = vector3<real>::determinant(va, vb, p)/detT;
    
    return bc;
  }

  real distance(const sphericalcoordinate &other) const
  {
    real hsinphi = sin((phi - other.phi)/2.0);
    real hsintheta = sin((theta - other.theta)/2.0);
      
    return 2.0 * asin(sqrt(hsinphi * hsinphi + sin(phi)*sin(other.phi)*(hsintheta*hsintheta)));
  }

  bool writetext(FILE *fp) const
  {
    double dphi = scalartodouble<real>(phi);
    double dtheta = scalartodouble<real>(theta);

    fprintf(fp, "%.15g %.15g", dphi, dtheta);
    return true;
  }

  bool readtext(FILE *fp)
  {
    double dphi, dtheta;
    if (fscanf(fp, "%lf %lf", &dphi, &dtheta) != 2) {
      return false;
    }

    phi = doubletoscalar<real>(dphi);
    theta = doubletoscalar<real>(dtheta);
    return true;
  }
  
  bool writebinary(FILE *fp) const
  {
    return value_write<real>(fp, phi) && value_write<real>(fp, theta);
  }

  bool readbinary(FILE *fp)
  {
    return value_read<real>(fp, phi) && value_read<real>(fp, theta);
  }

  bool write(std::ostream &s)
  {
    s.write((char*)&phi, sizeof(real));
    s.write((char*)&theta, sizeof(real));

    return true;
  }
  
  bool read(std::istream &s)
  {
    s.read((char*)&phi, sizeof(real));
    s.read((char*)&theta, sizeof(real));

    return true;
  }

  static void cartesiantospherical(const vector3<real> &v,
				   sphericalcoordinate &a)
  {
    real r = v.length();
    a.phi = acos(v.z/r);
    a.theta = atan2(v.y, v.x);
  }
  
  static void cartesiantospherical(const vector3<real> &v,
				   real &phi, real &theta, real &r)
  {
    r = v.length();
    phi = acos(v.z/r);
    theta = atan2(v.y, v.x);
  }

  static void sphericaltocartesian(const sphericalcoordinate &a,
				   vector3<real> &v)
  {
    v.x = cos(a.theta) * sin(a.phi);
    v.y = sin(a.theta) * sin(a.phi);
    v.z = cos(a.phi);
  }
  
  static void sphericaltocartesian(real phi, real theta,
				   vector3<real> &v)
  {
    v.x = cos(theta) * sin(phi);
    v.y = sin(theta) * sin(phi);
    v.z = cos(phi);
  }

  friend bool operator==(const sphericalcoordinate &a,
			 const sphericalcoordinate &b)
  {
    return (a.phi == b.phi) && (a.theta == b.theta);
  }

  friend bool operator!=(const sphericalcoordinate &a,
			 const sphericalcoordinate &b)
  {
    return (a.phi != b.phi) || (a.theta != b.theta);
  }

  real phi;   // Colatitude
  real theta; // Longitude

};

#endif // coordinate_hpp
