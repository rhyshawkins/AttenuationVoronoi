#pragma once
#ifndef rng_h
#define rng_h

#include <memory>

//
// A simple wrapper around gsl random number generator
//
class Rng {
public:

  Rng(int seed);
  ~Rng();

  //
  // Integer random
  //
  int uniform(int n);  // k in (0 .. n - 1) uniform
  int jeffreys(int n); // k in (1 .. n) proportional to 1/k

  int select(int nweights, double *weights);
  void shuffle(int nitems, int *items);
  
  //
  // Floating point random
  //
  double uniform();
  double normal(double sigma);
  double gamma(double a, double b);

  //
  // PDF
  //
  static double pdf_normal(double x, double mean, double sigma);

private:

  class impl;
  std::unique_ptr<impl> pimpl;

};

#endif // rng_h
  
