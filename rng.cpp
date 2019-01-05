
#include "rng.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Rng::impl {
public:

  impl(int seed) :
    rng(gsl_rng_alloc(gsl_rng_taus))
  {
    gsl_rng_set(rng, seed);
  }
  
  ~impl()
  {
    gsl_rng_free(rng);
  }

  gsl_rng *rng;
};

Rng::Rng(int seed) :
  pimpl(new impl(seed))
{
}

Rng::~Rng()
{
}

int
Rng::uniform(int n)
{
  return gsl_rng_uniform_int(pimpl->rng, n);
}

int
Rng::jeffreys(int n)
{
  double c;
  c = 0.0;
  for (int i = 1; i <= n; i ++) {
    c += 1.0/(double)i;
  }
  if (c == 0.0) {
    return -1;
  }

  c = 1.0/c;
  double u = uniform();

  int k = 1;
  double pk = c/(double)k;
  while (u > pk) {
    u -= pk;
    k ++;
    pk = c/(double)k;
  }

  return k;
}

int
Rng::select(int nweights, double *weights)
{
  double sum = 0.0;
  for (int i = 0; i < nweights; i ++) {
    sum += weights[i];
  }

  double u = uniform() * sum;

  for (int i = 0; i < nweights; i ++) {
    if (u < weights[i]) {
      return i;
    }

    u -= weights[i];
  }

  return nweights - 1;
}

void
Rng::shuffle(int nitems, int *items)
{
  gsl_ran_shuffle(pimpl->rng, items, nitems, sizeof(int));
}

double
Rng::uniform()
{
  return gsl_rng_uniform(pimpl->rng);
}

double
Rng::normal(double sigma)
{
  return gsl_ran_gaussian_ziggurat(pimpl->rng, sigma);
}

double
Rng::gamma(double a, double b)
{
  return gsl_ran_gamma(pimpl->rng, a, b);
}

double
Rng::pdf_normal(double x, double mean, double sigma)
{
  return gsl_ran_gaussian_pdf(x - mean, sigma);
}
