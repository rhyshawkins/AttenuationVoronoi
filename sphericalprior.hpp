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
#ifndef sphericalprior_hpp
#define sphericalprior_hpp

#include <map>
#include <string>

#include "rng.hpp"

class SphericalPrior {
public:

  SphericalPrior();
  virtual ~SphericalPrior();

  virtual void sample(Rng &rng, double &phi, double &theta) = 0;

  virtual bool valid(double phi, double theta) = 0;

  virtual double pdf(double phi, double theta) = 0;

  virtual double logpdf(double phi, double theta) = 0;

  static SphericalPrior *load(FILE *fp);

  typedef SphericalPrior* (*spherical_prior_reader_f)(FILE *fp);

  static bool register_prior(const char *name, spherical_prior_reader_f reader);
  
private:

  static std::map<std::string, spherical_prior_reader_f> readers;
};
  
//
// Abstract proposal
//
class SphericalProposal {
public:

  SphericalProposal(SphericalPrior &prior);
  virtual ~SphericalProposal();

  virtual SphericalPrior *get_prior();
  
  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldphi,
		       double oldtheta,
		       double &newphi,
		       double &newtheta,
		       double &logpriorratio) = 0;

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldphi,
			      double oldtheta,
			      double newphi,
			      double newtheta) = 0;
  
  virtual double log_proposal_ratio(Rng &rng, double temperature,
				    double oldphi,
				    double oldtheta,
				    double newphi,
				    double newtheta) = 0;

  static SphericalProposal* load(FILE *fp, SphericalPrior &prior);

  typedef SphericalProposal* (*spherical_proposal_reader_f)(FILE *fp, SphericalPrior &prior);
  
  static bool register_spherical_proposal(const char *name, spherical_proposal_reader_f reader);

protected:

  SphericalPrior &prior;

private:

  static std::map<std::string, spherical_proposal_reader_f> readers;

};

class SphericalPriorProposal {
public:

  SphericalPriorProposal(SphericalPrior *prior, SphericalProposal *proposal);
  ~SphericalPriorProposal();
  
  void sample(Rng &rng, double &phi, double &theta);

  double pdf(double phi, double theta);

  double logpdf(double phi, double theta);

  bool propose(Rng &rng,
	       double temperature,
	       double oldphi, double oldtheta,
	       double &newphi, double &newtheta,
	       double &logpriorratio);
  
  double log_proposal_ratio(Rng &rng, double temperature,
			    double oldphi, double oldtheta,
			    double newphi, double newtheta);

  double log_proposal(Rng &rng, double temperature,
		      double oldphi, double oldtheta,
		      double newphi, double newtheta);

  SphericalPrior *get_prior()
  {
    return prior;
  }
  
  SphericalProposal *get_proposal()
  {
    return proposal;
  }

  static SphericalPriorProposal* load(const char *filename);

private:

  SphericalPrior *prior;
  SphericalProposal *proposal;

};

//
// Standard Priors
//

class UniformSphericalPrior : public SphericalPrior { 
public:

  UniformSphericalPrior();
  ~UniformSphericalPrior();

  virtual bool valid(double phi, double theta);
  
  virtual void sample(Rng &rng, double &phi, double &theta);

  virtual double pdf(double phi, double theta);

  virtual double logpdf(double phi, double theta);

  static SphericalPrior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
  double vmin;
  double vmax;
  
};

class CosineSphericalPrior : public SphericalPrior {
public:
  CosineSphericalPrior(double phi0, double theta0, double delta);
  ~CosineSphericalPrior();

  virtual bool valid(double phi, double theta);
  
  virtual void sample(Rng &rng, double &phi, double &theta);

  virtual double pdf(double phi, double theta);

  virtual double logpdf(double phi, double theta);

  static SphericalPrior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;

  double phi0;
  double theta0;
  double delta;
  
};
  
extern double vonMisesPhi(double u, double kappa);

class VonMisesSphericalPrior : public SphericalPrior {
public:

  VonMisesSphericalPrior(double phi0, double theta0, double kappa);
  ~VonMisesSphericalPrior();

  virtual bool valid(double phi0, double theta0);
  
  virtual void sample(Rng &rng, double &phi, double &theta);

  virtual double pdf(double phi, double theta);

  virtual double logpdf(double phi, double theta);

  static SphericalPrior* reader(FILE *fp);
  
private:

  static const bool REGISTRATION;
  
  double phi0;
  double theta0;
  double kappa;
  
};
  

//
// Standard Proposals
//

class PriorSampleSphericalProposal : public SphericalProposal {
public:
  
  PriorSampleSphericalProposal(SphericalPrior &prior);
  ~PriorSampleSphericalProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldphi,
		       double oldtheta,
		       double &newphi,
		       double &newtheta,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldphi,
			      double oldtheta,
			      double newphi,
			      double newtheta);
  
  virtual double log_proposal_ratio(Rng &rng,
				    double temperature,
				    double oldphi,
				    double oldtheta,
				    double newphi,
				    double newtheta);

  static SphericalProposal *read(FILE *fp, SphericalPrior &prior);

private:

  static const bool REGISTRATION;

};

class VonMisesSphericalProposal : public SphericalProposal {
public:

  VonMisesSphericalProposal(SphericalPrior &prior, double kappa);
  ~VonMisesSphericalProposal();

  virtual bool propose(Rng &rng,
		       double temperature,
		       double oldphi,
		       double oldtheta,
		       double &newphi,
		       double &newtheta,
		       double &logpriorratio);

  virtual double log_proposal(Rng &rng, double temperature,
			      double oldphi,
			      double oldtheta,
			      double newphi,
			      double newtheta);

  virtual double log_proposal_ratio(Rng &rng,
				    double temperature,
				    double oldphi,
				    double oldtheta,
				    double newphi,
				    double newtheta);

  static SphericalProposal *read(FILE *fp, SphericalPrior &prior);

private:

  static const bool REGISTRATION;
  
  double kappa;

};
  

#endif // prior_hpp
