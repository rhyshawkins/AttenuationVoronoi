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
#ifndef attenuationdataS2_hpp
#define attenuationdataS2_hpp

#include <vector>

#include <stdio.h>

#include "coordinate.hpp"
#include "sphericalvoronoimodel.hpp"
#include "velocitymodel.hpp"

typedef double (*synthetic_model_f)(double phi, double theta);

template
< typename value >
class dataS2 {
public:
  typedef sphericalcoordinate<value> coord_t;

  //
  // Phi is colatitude in radians (0 (north pole) .. pi (south pole))
  // Theta is longitude in radians
  //
  dataS2(double _phi, double _theta, double _r, double _vp, double _vs) :
    phi(_phi),
    theta(_theta),
    r(_r),
    vp(_vp),
    vs(_vs),
    distance(0.0)
  {
  }

  double compute_distance(const dataS2<value> &rhs)
  {
    vector3<value> a;
    vector3<value> b;

    coord_t::sphericaltocartesian(phi, theta, a);
    a *= r;
    coord_t::sphericaltocartesian(rhs.phi, rhs.theta, b);
    b *= rhs.r;
    
    return (a - b).length();
  }

  double phi, theta, r;
  double vp, vs;
  double distance;
};

template
< typename value >
class pathS2 {
public:

  typedef sphericalcoordinate<value> coord_t;

  pathS2(double _tstar, double _noise) :
    tstar(_tstar),
    noise(_noise)
  {
  }

  void compute_distances()
  {
    for (size_t i = 1; i < points.size(); i ++) {

      double distance = points[i - 1].compute_distance(points[i]);

      points[i - 1].distance += distance/2.0;
      points[i].distance += distance/2.0;

    }
  }

  double compute_mean_Q()
  {
    double tt = 0.0;
    for (auto &p : points) {
      tt += p.distance/p.vp;
    }

    return tt/tstar;
  }

  value predicted_tstar_direct(const sphericalvoronoimodel<value> &model)
  {
    value tstar = 0.0;
    for (auto &d : points) {

      value Q = model.value_at_point(coord_t(d.phi, d.theta));

      tstar += d.distance/(Q * d.vp);
      
    }

    return tstar;
  }

  value predicted_tstar_synthetic(synthetic_model_f model)
  {
    value tstar = 0.0;
    for (auto &d : points) {

      value Q = model(d.phi, d.theta);

      tstar += d.distance/(Q * d.vp);
      
    }

    return tstar;
  }

  double tstar;
  double noise;
  std::vector<dataS2<value>> points;
};

template
< typename value >
class attenuationdataS2 {
public:

  attenuationdataS2(const char *filename) :
    phimin(1e99),
    phimax(-1e99),
    thetamin(1e99),
    thetamax(-1e99),
    Qmin(1e99),
    Qmax(-1e99),
    Qcount(0),
    Qmean(0.0)
  {

    FILE *fp;

    fp = fopen(filename, "r");
    if (fp == NULL) {
      throw ATTENUATIONEXCEPTION("Failed to open %s\n", filename);
    }

    while (!feof(fp)) {

      double tstar;
      double noise;
      int n;

      if (fscanf(fp, "%lf %lf %d\n", &tstar, &noise, &n) != 3) {
	if (feof(fp)) {
	  break;
	} else {
	  throw ATTENUATIONEXCEPTION("Failed to read data\n");
	}
      }

      pathS2<value> p(tstar, noise);
      
      for (int i = 0; i < n; i ++) {

	double lon, lat, r;
	if (fscanf(fp, "%lf %lf %lf\n", &lon, &lat, &r) != 3) {
	  throw ATTENUATIONEXCEPTION("Failed to read line\n");
	}
	
	double phi = (90.0 - lat) * M_PI/180.0;
	double theta = lon * M_PI/180.0;
	
	//
	// Roundabout way to ensure theta is between -pi .. pi
	//
	double x = cos(theta);
	double y = sin(theta);
	theta = atan2(y, x);
	
	if (phi < 0.0 || phi > M_PI) {
	  throw ATTENUATIONEXCEPTION("Latitude out of range: %f\n", lat);
	}
	
	phimin = std::min<double>(phimin, phi);
	phimax = std::max<double>(phimax, phi);
	
	thetamin = std::min<double>(thetamin, theta);
	thetamax = std::max<double>(thetamax, theta);

	//
	// Constants used temporarily for velocities
	//
	p.points.push_back(dataS2<value>(phi, theta, r, pwave_velocity<value>(r), 3.0));

      }

      p.compute_distances();
      double Q = p.compute_mean_Q();
      Qmin = std::min<double>(Qmin, Q);
      Qmax = std::max<double>(Qmax, Q);

      Qcount ++;
      double delta = Q - Qmean;
      Qmean += delta/(double)Qcount;
      
      data.push_back(p);
    }
    
    fclose(fp);
  }

  ~attenuationdataS2()
  {
  }

  value likelihood(const sphericalvoronoimodel<value> &model,
		   double lambda,
		   value *residuals)
  {
    value sum = 0.0;
    int i = 0;

    for (auto &d : data) {

      value pred = d.predicted_tstar_direct(model);
      value res = pred - d.tstar;
      double sigma = d.noise * lambda;

      residuals[i] = res;
      
      sum += res*res/(2.0 * sigma * sigma);
      i ++;
    }

    return sum;
  }

  value likelihood_partial(const sphericalvoronoimodel<value> &model,
			   double lambda,
			   int offset,
			   int size,
			   value *residuals)
  {
    value sum = 0.0;

    for (int i = 0; i < size; i ++) {

      auto &d = data[offset + i];

      value pred = d.predicted_tstar_direct(model);
      value res = pred - d.tstar;
      double sigma = d.noise * lambda;

      residuals[i] = res;
      
      sum += res*res/(2.0 * sigma * sigma);
      
    }

    return sum;
  }

  double phimin, phimax;
  double thetamin, thetamax;
  double Qmin, Qmax;

  int Qcount;
  double Qmean;
  
  std::vector<pathS2<value>> data;
};

#endif // attenuationdataS2_hpp
