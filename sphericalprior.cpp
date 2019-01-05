
#include <string.h>
#include <math.h>

#include "sphericalprior.hpp"

#include "attenuationexception.hpp"

#include "coordinate.hpp"

extern "C" {
  #include "slog.h"
};

typedef vector3<double> vec3_t;
typedef sphericalcoordinate<double> coord_t;

std::map<std::string, SphericalPrior::spherical_prior_reader_f> SphericalPrior::readers;

static double distance(double phi0, double theta0,
		       double phi1, double theta1);

static void rotatefrompole(double phipole,
			   double thetapole,
			   double phi0,
			   double theta0,
			   double &phi,
			   double &theta);

double vonMisesPhi(double u, double kappa);

static double
newtonsolve(double (*f)(void *user, double x),
	    double (*fprime)(void *user, double x),
	    double xmin,
	    double xmax,
	    double y,
	    double threshold,
	    int maxsteps,
	    void *user);

SphericalPrior::SphericalPrior()
{
}

SphericalPrior::~SphericalPrior()
{
}

SphericalPrior *
SphericalPrior::load(FILE *fp)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "sphericalprior") != 0) {
    fprintf(stderr, "SphericalPrior::load: exptected 'prior', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, SphericalPrior::spherical_prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "SphericalPrior::load: no prior name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp);

  }
}

bool
SphericalPrior::register_prior(const char *name, spherical_prior_reader_f reader)
{
  std::map<std::string, spherical_prior_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw ATTENUATIONEXCEPTION("Duplicate prior registered: %s\n", name);
}


std::map<std::string, SphericalProposal::spherical_proposal_reader_f> SphericalProposal::readers;

SphericalProposal::SphericalProposal(SphericalPrior &_prior) :
  prior(_prior)
{
}

SphericalProposal::~SphericalProposal()
{
}

SphericalPrior *
SphericalProposal::get_prior()
{
  return &prior;
}

SphericalProposal*
SphericalProposal::load(FILE *fp, SphericalPrior &prior)
{
  char type[256];
  char name[256];

  if (fscanf(fp, "%s %s\n", type, name) != 2) {
    return nullptr;
  }

  if (strcmp(type, "sphericalproposal") != 0) {
    fprintf(stderr, "SphericalProposal::load: exptected 'proposal', got '%s'\n", type);
    return nullptr;
  }

  std::map<std::string, spherical_proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    fprintf(stderr, "SphericalProposal::load: no proposal name '%s'\n", name);
    return nullptr;
  } else {

    return (i->second)(fp, prior);

  }
}

bool
SphericalProposal::register_spherical_proposal(const char *name, spherical_proposal_reader_f reader)
{
  std::map<std::string, spherical_proposal_reader_f>::iterator i = readers.find(name);
  if (i == readers.end()) {
    readers[name] = reader;
    return true;
  }

  throw ATTENUATIONEXCEPTION("Duplicate proposal registered: %s\n", name);
}


SphericalPriorProposal::SphericalPriorProposal(SphericalPrior *_prior, SphericalProposal *_proposal) :
  prior(_prior),
  proposal(_proposal)
{
}

SphericalPriorProposal::~SphericalPriorProposal()
{
  delete prior;
  delete proposal;
}
  
void
SphericalPriorProposal::sample(Rng &rng, double &phi, double &theta)
{
  
  prior->sample(rng, phi, theta);
}

double
SphericalPriorProposal::pdf(double phi, double theta)
{
  return prior->pdf(phi, theta);
}

double
SphericalPriorProposal::logpdf(double phi, double theta)
{
  return prior->logpdf(phi, theta);
}

bool
SphericalPriorProposal::propose(Rng &rng,
				double temperature,
				double oldphi,
				double oldtheta,
				double &newphi,
				double &newtheta,
				double &logpriorratio)
{
  return proposal->propose(rng, temperature, oldphi, oldtheta, newphi, newtheta, logpriorratio);
}
  
double
SphericalPriorProposal::log_proposal(Rng &rng,
					   double temperature,
					   double oldphi,
					   double oldtheta,
					   double newphi,
					   double newtheta)
{
  return proposal->log_proposal(rng, temperature, oldphi, oldtheta, newphi, newtheta);
}

double
SphericalPriorProposal::log_proposal_ratio(Rng &rng,
					   double temperature,
					   double oldphi,
					   double oldtheta,
					   double newphi,
					   double newtheta)
{
  return proposal->log_proposal_ratio(rng, temperature, oldphi, oldtheta, newphi, newtheta);
}

SphericalPriorProposal*
SphericalPriorProposal::load(const char *filename)
{
  FILE *fp;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "SphericalPriorProposal::load: failed to open file for reading\n");
    return nullptr;
  }

  SphericalPrior *prior = SphericalPrior::load(fp);
  if (prior == nullptr) {
    return nullptr;
  }

  SphericalProposal *proposal = SphericalProposal::load(fp, *prior);
  if (proposal == nullptr) {
    return nullptr;
  }

  fclose(fp);

  return new SphericalPriorProposal(prior, proposal);
}

//
// Uniform Prior
//

const bool UniformSphericalPrior::REGISTRATION =
  SphericalPrior::register_prior("UniformSpherical", UniformSphericalPrior::reader);

UniformSphericalPrior::UniformSphericalPrior()
{
}

UniformSphericalPrior::~UniformSphericalPrior()
{
}

bool
UniformSphericalPrior::valid(double phi, double theta)
{
  return true;
}
  
void
UniformSphericalPrior::sample(Rng &rng, double &phi, double &theta)
{
  theta = 2.0 * M_PI * rng.uniform();
  phi = acos(2.0 * rng.uniform() - 1.0);
}

double
UniformSphericalPrior::pdf(double phi, double theta)
{
  return 1.0/(4.0 * M_PI);
}
  

double
UniformSphericalPrior::logpdf(double phi, double theta)
{
  return -log(4.0 * M_PI);
}

SphericalPrior*
UniformSphericalPrior::reader(FILE *fp)
{
  return new UniformSphericalPrior();
}
  
const bool CosineSphericalPrior::REGISTRATION =
  SphericalPrior::register_prior("CosineSpherical", CosineSphericalPrior::reader);

CosineSphericalPrior::CosineSphericalPrior(double _phi0,
					   double _theta0,
					   double _delta) :
  phi0(_phi0),
  theta0(_theta0),
  delta(_delta)
{
}

CosineSphericalPrior::~CosineSphericalPrior()
{
}

bool
CosineSphericalPrior::valid(double phi, double theta)
{
  return distance(phi0, theta0, phi, theta) < delta;
}

void
CosineSphericalPrior::sample(Rng &rng, double &phi, double &theta)
{
  //
  // PDF is (cos(phi * M_PI/delta) + 1)/delta
  // CDF is (phi/delta + 1/M_PI*sin(phi * M_PI/delta)
  //
  // Can generate uniform t and then use newton method to find root of CDF - t
  //
  auto cdf = [](void *user, double x) -> double {
    double delta = ((CosineSphericalPrior*)user)->delta;
    return (sin(M_PI * x/delta)/M_PI + x/delta);
  };
  
  auto pdf = [](void *user, double x) -> double {
    double delta = ((CosineSphericalPrior*)user)->delta;
    return ((cos(M_PI * x/delta) + 1.0)/delta);
  };

  double u = rng.uniform();

  double phipole = newtonsolve(cdf,
			       pdf,
			       0.0,
			       delta,
			       u,
			       1.0e-9,
			       1000,
			       this);

  double thetapole = rng.uniform() * 2.0 * M_PI;

  rotatefrompole(phipole, thetapole, phi0, theta0, phi, theta);
}

double
CosineSphericalPrior::pdf(double phi, double theta)
{
  double d = distance(phi0, theta0, phi, theta);

  if (d < delta) {
    return (cos(d * M_PI/delta) + 1.0)/(delta * 2.0 * M_PI);
  } else {
    return 0.0;
  }
}

double
CosineSphericalPrior::logpdf(double phi, double theta)
{
  double p = pdf(phi, theta);
  if (p > 0.0) {
    return log(p);
  }
  return p;
}

SphericalPrior*
CosineSphericalPrior::reader(FILE *fp)
{
  double phi0, theta0, delta;
  if (fscanf(fp, "%lf %lf %lf", &phi0, &theta0, &delta) != 3) {
    ERROR("Failed to read parameters for prior\n");
    return nullptr;
  }

  return new CosineSphericalPrior(phi0, theta0, delta);
}
  
//
// von Mises Prior
//
const bool VonMisesSphericalPrior::REGISTRATION = SphericalPrior::register_prior("VonMisesSpherical", VonMisesSphericalPrior::reader);

VonMisesSphericalPrior::VonMisesSphericalPrior(double _phi0, double _theta0, double _kappa) :
  phi0(_phi0),
  theta0(_theta0),
  kappa(_kappa)
{
}

VonMisesSphericalPrior::~VonMisesSphericalPrior()
{
}

bool
VonMisesSphericalPrior::valid(double phi0, double theta0)
{
  return true;
}
  
void
VonMisesSphericalPrior::sample(Rng &rng, double &phi, double &theta)
{
  double phipole = vonMisesPhi(rng.uniform(), kappa);
  double thetapole = 2.0 * M_PI * rng.uniform();

  rotatefrompole(phipole, thetapole, phi0, theta0, phi, theta);
}

double
VonMisesSphericalPrior::pdf(double phi, double theta)
{
  double delta = distance(phi0, theta0, phi, theta);

  return kappa/(4.0 * M_PI * sinh(kappa)) * exp(cos(delta) * kappa) * sin(delta);
}

double
VonMisesSphericalPrior::logpdf(double phi, double theta)
{
  return log(pdf(phi, theta));
}

SphericalPrior*
VonMisesSphericalPrior::reader(FILE *fp)
{
  double phi0;
  double theta0;
  double kappa;

  if (fscanf(fp, "%lf %lf %lf", &phi0, &theta0, &kappa) != 3) {
    ERROR("Failed to read parameters\n");
    return nullptr;
  }

  return new VonMisesSphericalPrior(phi0, theta0, kappa);
}

//
//
//
const bool PriorSampleSphericalProposal::REGISTRATION = SphericalProposal::register_spherical_proposal("PriorSampleSpherical", PriorSampleSphericalProposal::read);

PriorSampleSphericalProposal::PriorSampleSphericalProposal(SphericalPrior &_prior) :
  SphericalProposal(_prior)
{
}

PriorSampleSphericalProposal::~PriorSampleSphericalProposal()
{
}

bool
PriorSampleSphericalProposal::propose(Rng &rng,
				      double temperature,
				      double oldphi,
				      double oldtheta,
				      double &newphi,
				      double &newtheta,
				      double &logpriorratio)
{
  prior.sample(rng, newphi, newtheta);
  logpriorratio = prior.logpdf(newphi, newtheta) - prior.logpdf(oldphi, oldtheta);
  
  return true;
}

double
PriorSampleSphericalProposal::log_proposal(Rng &rng,
					   double temperature,
					   double oldphi,
					   double oldtheta,
					   double newphi,
					   double newtheta)
{
  return prior.logpdf(newphi, newtheta);
}

double
PriorSampleSphericalProposal::log_proposal_ratio(Rng &rng,
						 double temperature,
						 double oldphi,
						 double oldtheta,
						 double newphi,
						 double newtheta)
{
  return prior.logpdf(newphi, newtheta) - prior.logpdf(oldphi, oldtheta);
}

SphericalProposal *
PriorSampleSphericalProposal::read(FILE *fp, SphericalPrior &prior)
{
  return new PriorSampleSphericalProposal(prior);
}

//
// von Mises Spherical Proposal
//

const bool VonMisesSphericalProposal::REGISTRATION = SphericalProposal::register_spherical_proposal("VonMisesSpherical", VonMisesSphericalProposal::read);

VonMisesSphericalProposal::VonMisesSphericalProposal(SphericalPrior &prior, double _kappa) :
  SphericalProposal(prior),
  kappa(_kappa)
{
}

VonMisesSphericalProposal::~VonMisesSphericalProposal()
{
}

bool
VonMisesSphericalProposal::propose(Rng &rng,
				   double temperature,
				   double oldphi,
				   double oldtheta,
				   double &newphi,
				   double &newtheta,
				   double &logpriorratio)
{
  double phipole = vonMisesPhi(rng.uniform(), kappa);
  double thetapole = 2.0 * M_PI * rng.uniform();

  rotatefrompole(phipole, thetapole, oldphi, oldtheta, newphi, newtheta);

  if (prior.valid(newphi, newtheta)) {
  
    logpriorratio = prior.logpdf(newphi, newtheta) - prior.logpdf(oldphi, oldtheta);
    return true;

  } else {

    return false;
  }
}

double
VonMisesSphericalProposal::log_proposal(Rng &rng, double temperature,
					double oldphi,
					double oldtheta,
					double newphi,
					double newtheta)
{
  double delta = distance(oldphi, oldtheta, newphi, newtheta);

  return log(kappa/(4.0 * M_PI * sinh(kappa)) * exp(cos(delta) * kappa) * sin(delta));
}

double VonMisesSphericalProposal::log_proposal_ratio(Rng &rng,
						     double temperature,
						     double oldphi,
						     double oldtheta,
						     double newphi,
						     double newtheta)
{
  return 0.0;
}

SphericalProposal *
VonMisesSphericalProposal::read(FILE *fp, SphericalPrior &prior)
{
  double kappa;

  if (fscanf(fp, "%lf", &kappa) != 1) {
    ERROR("Failed to read parameters\n");
    return nullptr;
  }

  return new VonMisesSphericalProposal(prior, kappa);
}

//
// Helper functions
//

static double distance(double phi0, double theta0,
		       double phi1, double theta1)
{
  double hsinphi = sin((phi0 - phi1)/2.0);
  double hsintheta = sin((theta0 - theta1)/2.0);
  
  return 2.0 * asin(sqrt(hsinphi * hsinphi + sin(phi0)*sin(phi1)*(hsintheta*hsintheta)));
}

static void rotatefrompole(double phipole,   // Coordinates of point with pole as origin
			   double thetapole,
			   double phi0,      // New origin
			   double theta0,
			   double &phi,
			   double &theta)
{
  //
  // Early bailout
  //
  if (phi0 == 0.0) {
    phi = phipole;
    theta = thetapole;
    return;
  }

  if (fabs(phi0 - M_PI) < 1.0e-9) {
    // Assumes rotational invariance.
    phi = M_PI - phipole;
    theta = thetapole;
    return;
  }

  //
  // Construct rotation from North pole to phi0, theta0
  //
  vec3_t pole(0.0, 0.0, 1.0);
  vec3_t newpole;
  coord_t::sphericaltocartesian(phi0, theta0, newpole);
  vec3_t axis;
  double angle;

  if (!computeaxisangle(pole, newpole, axis, angle)) {
    throw ATTENUATIONEXCEPTION("Failed to compute axis/angle");
  }

  vec3_t c;
  coord_t::sphericaltocartesian(phipole, thetapole, c);

  vec3_t cp = rotate(c, axis, angle);

  double r;
  coord_t::cartesiantospherical(cp, phi, theta, r);
}

double vonMisesPhi(double u, double kappa)
{
  if (kappa > 0.0) {
    return acos(1.0 + log(u + (1.0 - u)*exp(-2.0 * kappa))/kappa);
  } else if (kappa == 0.0) {
    return acos(2.0 * u - 1.0);
  } else {
    throw ATTENUATIONEXCEPTION("kappa must be greater or equal to zero.");
  }
}

static double
newtonsolve(double (*f)(void *user, double x),
	    double (*fprime)(void *user, double x),
	    double xmin,
	    double xmax,
	    double y,
	    double threshold,
	    int maxsteps,
	    void *user)
{
  double x = (xmin + xmax)/2.0;
  double yt = f(user, x);
  int i;

  i = 0;
  while ((i < maxsteps) && (fabs(yt - y) > threshold)) {

    double g = fprime(user, x);
    double dx;
    
    if (g == 0.0) {
      //
      // Since the cdf's should be monotonicly increasing, this should only happen at the boundaries
      //
      if (y < yt) {
	dx = 1.0e-6;
      } else {
	dx = -1.0e-6;
      }
    } else {
      dx = (y - yt)/g;
    }

    x += dx;
    if (x < xmin) {
      x = xmin + threshold;
    }
    if (x > xmax) {
      x = xmax - threshold;
    }

    yt = f(user, x);

    i ++;
  }

  return x;
}



