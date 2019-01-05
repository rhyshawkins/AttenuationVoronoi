#pragma once
#ifndef globalspherical_hpp
#define globalspherical_hpp

#include <mpi.h>

#include "sphericalvoronoimodel.hpp"

#include "attenuationdataS2.hpp"
#include "prior.hpp"
#include "sphericalprior.hpp"

#include "hierarchical_model.hpp"

extern "C" {
  #include "slog.h"
};

template
<typename value>
class globalS2Voronoi {
public:

  typedef sphericalcoordinate<value> coord_t;

  globalS2Voronoi(const char *input,
		  const char *initial_model,
		  const char *prior_file,
		  const char *hierarchical_prior_file,
		  const char *position_prior_file,
		  const char *birthdeath_proposal_file,
		  int _maxcells,
		  double _lambda,
		  double _temperature,
		  int seed,
		  bool posterior,
		  bool logspace) :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1),
    mpi_counts(nullptr),
    mpi_offsets(nullptr),
    data(nullptr),
    model(nullptr),
    prior(nullptr),
    positionprior(nullptr),
    hierarchicalprior(nullptr),
    birthdeathvalueproposal(nullptr),
    birthdeathpositionproposal(nullptr),
    hierarchical(new singlescaling_hierarchical_model(_lambda)),
    temperature(_temperature),
    residual_size(0),
    mean_residual_n(0),
    mean_residuals(nullptr),
    residuals(nullptr),
    last_valid_residuals(nullptr),
    maxcells(_maxcells),
    random(seed)
  {

    if (prior_file == nullptr) {
      Prior *p = new UniformPrior(-1.5, 1.5);
      Proposal *pp = new GaussianProposal(*p, 0.05);
      prior = new PriorProposal(p, pp);
    } else {
      prior = PriorProposal::load(prior_file);
    }

    if (hierarchical_prior_file == nullptr) {
      Prior *p = new UniformPrior(0.1, 5.0);
      Proposal *pp = new GaussianProposal(*p, 0.05);
      hierarchicalprior = new PriorProposal(p, pp);
    } else {
      hierarchicalprior = PriorProposal::load(hierarchical_prior_file);
    }

    if (position_prior_file == nullptr) {
      SphericalPrior *p = new UniformSphericalPrior();
      SphericalProposal *pp = new VonMisesSphericalProposal(*p, 1.0);
      positionprior = new SphericalPriorProposal(p, pp);
    } else {
      positionprior = SphericalPriorProposal::load(position_prior_file);
    }

    if (birthdeath_proposal_file == nullptr) {
      birthdeathvalueproposal = new PriorSampleProposal(*(prior->get_prior()));
      birthdeathpositionproposal = new PriorSampleSphericalProposal(*(positionprior->get_prior()));
    } else {
      FILE *fp = fopen(birthdeath_proposal_file, "r");
      if (fp == NULL) {
        throw ATTENUATIONEXCEPTION("Failed to open birth/death proposal file\n");
      }
      
      birthdeathvalueproposal = Proposal::load(fp, *(prior->get_prior()));
      birthdeathpositionproposal = SphericalProposal::load(fp, *(positionprior->get_prior()));
      
      fclose(fp);

      if (birthdeathvalueproposal == nullptr ||
          birthdeathpositionproposal == nullptr) {
        throw ATTENUATIONEXCEPTION("Failed to load birth/death proposal file\n");
      }
    }

    value Qmean = 1.0;
    
    if (!posterior) {

      data = new attenuationdataS2<value>(input);

      INFO(" Qmin: %10.6f", data->Qmin);
      INFO(" Qmax: %10.6f", data->Qmax);
      INFO("Qmean: %10.6f", data->Qmean);

      if (logspace) {
	Qmean = log(data->Qmean);
      } else {
	Qmean = data->Qmean;
      }

      residual_size = data->data.size();
      residuals = new value[residual_size];
      mean_residuals = new value[residual_size];
      last_valid_residuals = new value[residual_size];
      mean_residual_n = 0;
      for (int i = 0; i < residual_size; i ++) {
	residuals[i] = 0.0;
	mean_residuals[i] = 0.0;
	last_valid_residuals[i] = 0.0;
      }

    } else {

      data = nullptr;
      
    }
    
    model = new sphericalvoronoimodel<value>(logspace);

    if (initial_model == nullptr) {

      //
      // Initialize to mean Q computed from data with a single point at the North Pole
      //
      model->add_cell(coord_t(0.0, 0.0), Qmean);

    } else {

      if (!model->load(initial_model)) {
	throw ATTENUATIONEXCEPTION("Failed to load initial model from %s", initial_model);
      }

      INFO("Loaded model with %d cells", model->ncells());

    }
  }

  void initialize_mpi(MPI_Comm _communicator, double _temperature)
  {
    MPI_Comm_dup(_communicator, &communicator);

    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    temperature = _temperature;

    int observations = data->data.size();
    int processes = size;

    mpi_offsets = new int[size];
    mpi_counts = new int[size];
      
    for (int i = 0; i < size; i ++) {
      mpi_counts[i] = observations/processes;
      observations -= mpi_counts[i];
      processes --;
    }

    mpi_offsets[0] = 0;
    for (int i = 1; i < size; i ++) {
      mpi_offsets[i] = mpi_offsets[i - 1] + mpi_counts[i - 1];
    }

    if (mpi_offsets[size - 1] + mpi_counts[size - 1] != (int)data->data.size()) {
      throw ATTENUATIONEXCEPTION("Failed to distribute data points properly");
    }
  }

  value likelihood()
  {
    if (data) {
      if (communicator == MPI_COMM_NULL) {
	return data->likelihood(*model, hierarchical->get(0), residuals);
      } else {

	value plike = data->likelihood_partial(*model,
					       hierarchical->get(0),
					       mpi_offsets[rank],
					       mpi_counts[rank],
					       residuals + mpi_offsets[rank]);
	value sumlike;
	MPI_Reduce(&plike, &sumlike, 1, MPI_DOUBLE, MPI_SUM, 0, communicator);
	MPI_Bcast(&sumlike, 1, MPI_DOUBLE, 0, communicator);

	MPI_Allgatherv(residuals + mpi_offsets[rank],
		       mpi_counts[rank],
		       MPI_DOUBLE,
		       residuals,
		       mpi_counts,
		       mpi_offsets,
		       MPI_DOUBLE,
		       communicator);

	return sumlike;
	
      }
    } else {
      return 1.0;
    }
  }

  void accept()
  {
    for (int i = 0; i < residual_size; i ++) {
      last_valid_residuals[i] = residuals[i];
    }
    update_mean_residual();
  }

  void reject()
  {
    update_mean_residual();
  }

  void update_mean_residual()
  {
    mean_residual_n ++;

    for (int i = 0; i < residual_size; i ++) {
      value delta = last_valid_residuals[i] - mean_residuals[i];
      mean_residuals[i] += delta/(double)mean_residual_n;
    }
  }
  

  MPI_Comm communicator;
  int rank;
  int size;
  int *mpi_counts;
  int *mpi_offsets;

  attenuationdataS2<value> *data;
  sphericalvoronoimodel<value> *model;

  PriorProposal *prior;
  SphericalPriorProposal *positionprior;
  PriorProposal *hierarchicalprior;

  Proposal *birthdeathvalueproposal;
  SphericalProposal *birthdeathpositionproposal;

  hierarchical_model *hierarchical;
  double temperature;

  int residual_size;
  int mean_residual_n;
  value *mean_residuals;
  value *residuals;
  value *last_valid_residuals;

  int maxcells;
  
  Rng random;
  
};

#endif // globalspherical_hpp
