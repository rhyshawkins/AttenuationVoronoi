#pragma once
#ifndef perturbationS2Voronoi_hpp
#define perturbationS2Voronoi_hpp

#include <mpi.h>

#include "rng.hpp"
#include "prior.hpp"
#include "sphericalprior.hpp"
#include "hierarchical_model.hpp"

#include "sphericalvoronoimodel.hpp"

#include "chainhistoryVoronoi.hpp"

template
<
  typename coord,
  typename value
>
class deltaVoronoi;

template
<typename value>
class PerturbationS2Voronoi {
public:

  typedef sphericalcoordinate<value> coord_t;
  typedef deltaVoronoi<coord_t, value> delta_t;
  
  PerturbationS2Voronoi() :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1)
  {
  }
  
  virtual ~PerturbationS2Voronoi()
  {
  }

  void initialize_mpi(MPI_Comm _communicator)
  {
    MPI_Comm_dup(_communicator, &communicator);
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
  }
  
  virtual bool propose(int maxcells,
		       int nobs,
		       Rng &random,
		       PriorProposal &prior,
		       SphericalPriorProposal &position_prior,
		       sphericalvoronoimodel<value> &model,
		       PriorProposal &hierarchical_prior,
		       hierarchical_model &hierarchical,
		       double temperature,
		       double &log_prior_ratio,
		       delta_t *&perturbation) = 0;

  virtual double log_proposal_ratio(Rng &random,
				    PriorProposal &prior,
				    SphericalPriorProposal &position_prior,
				    sphericalvoronoimodel<value> &proposed_model,
				    PriorProposal &hierarchical_prior,
				    hierarchical_model &proposed_hierarchical,
				    double temperature) = 0;

  virtual void accept() = 0;

  virtual void reject(sphericalvoronoimodel<value> &model) = 0;

  virtual int proposal_count() const = 0;
  
  virtual int acceptance_count() const = 0;

  virtual const char *displayname() const = 0;

protected:

  bool primary()
  {
    return (communicator == MPI_COMM_NULL) || (rank == 0);
  }

  bool secondary()
  {
    return (communicator != MPI_COMM_NULL) && (rank != 0);
  }

  void communicate(bool &b)
  {
    if (communicator != MPI_COMM_NULL) {
      int i = (int)b;
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
      b = (bool)i;
    }
  }
  
  void communicate(int &i)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
    }
  }
  
  void communicate(value &d)
  {
    if (communicator != MPI_COMM_NULL) {
      double t = (double)d;
      MPI_Bcast(&t, 1, MPI_DOUBLE, 0, communicator);
      d = t;
    }
  }

  void communicate(coord_t &p)
  {
    if (communicator != MPI_COMM_NULL) {
      double t[2];
      t[0] = p.phi;
      t[1] = p.theta;
      MPI_Bcast(t, 2, MPI_DOUBLE, 0, communicator);
      p.phi = t[0];
      p.theta = t[1];
    }
  }

private:

  MPI_Comm communicator;
  int rank;
  int size;
  

};

#endif // perturbationS2Voronoi_hpp
