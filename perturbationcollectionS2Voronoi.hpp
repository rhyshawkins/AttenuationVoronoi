#pragma once
#ifndef perturbationcollectionS2Voronoi_hpp
#define perturbationcollectionS2_hpp

#include "globalS2Voronoi.hpp"
#include "chainhistoryVoronoi.hpp"
#include "perturbationS2Voronoi.hpp"

template
<typename value>
class PerturbationCollectionS2Voronoi {
public:

  typedef sphericalcoordinate<value> coord_t;
  typedef deltaVoronoi<coord_t, value> delta_t;
  
  
  PerturbationCollectionS2Voronoi() :
    communicator(MPI_COMM_NULL),
    rank(-1),
    size(-1),
    weight_sum(0.0),
    active(-1)
  {
  }
  
  ~PerturbationCollectionS2Voronoi()
  {
    for (auto &i : perturbations) {
      delete i.p;
    }
  }

  void initialize_mpi(MPI_Comm _communicator)
  {
    MPI_Comm_dup(_communicator, &communicator);
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);
  }
  
  void add(PerturbationS2Voronoi<value> *p, double weight)
  {
    if (weight <= 0.0) {
      throw ATTENUATIONEXCEPTION("Invalid weight: %f\n", weight);
    }
    
    weight_sum += weight;
    perturbations.push_back(WeightedPerturbation(p, weight));
  }

  bool propose(globalS2Voronoi<value> &g, double &log_prior_ratio, delta_t *&perturbation)
  {
    if (active >= 0) {
      throw ATTENUATIONEXCEPTION("Already have active perturbation\n");
    }
    
    if (primary()) {
      double u = g.random.uniform();
      active = 0;
      for (auto &wp: perturbations) {
	double w = wp.w/weight_sum;
	
	if (u < w) {
	  break;
	} else {
	  active ++;
	  u -= w;
	}
      }
    
      if (active >= (int)perturbations.size()) {
	throw ATTENUATIONEXCEPTION("Failed to determine active perturbation\n");
      }
    }

    communicate(active);
    
    WeightedPerturbation &wp = perturbations[active];
    
    int nobs = 0;
    if (g.data != nullptr) {
      nobs = g.data->data.size();
    }

    bool r = wp.p->propose(g.maxcells,
			   nobs,
			   g.random,
			   *g.prior,
			   *g.positionprior,
			   *g.model,
			   *g.hierarchicalprior,
			   *g.hierarchical,
			   g.temperature,
			   log_prior_ratio,
			   perturbation);
    if (!r) {
      active = -1;
    }
    
    return r;
  }

  double log_proposal_ratio(globalS2Voronoi<value> &g) 
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw ATTENUATIONEXCEPTION("Invalid active perturbation %d (%d)\n", active, (int)perturbations.size());
    }
    
    WeightedPerturbation &wp = perturbations[active];

    return wp.p->log_proposal_ratio(g.random,
				    *g.prior,
				    *g.positionprior,
				    *g.model,
				    *g.hierarchicalprior,
				    *g.hierarchical,
				    1.0);
  }

  void accept(globalS2Voronoi<value> &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw ATTENUATIONEXCEPTION("Invalid active perturbation %d (%d)\n", active, (int)perturbations.size());
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->accept();
    
    active = -1;
  }
  

  void reject(globalS2Voronoi<value> &g)
  {
    if (active < 0 || active >= (int)perturbations.size()) {
      throw ATTENUATIONEXCEPTION("Invalid active perturbation %d\n", active);
    }
    
    WeightedPerturbation &wp = perturbations[active];
    wp.p->reject(*g.model);
    
    active = -1;
  }
  

  void writeacceptancereport(FILE *fp)
  {
    fprintf(fp, "%s", generateacceptancereport().c_str());
  }

  std::string generateacceptancereport()
  {
    std::string s;
    char linebuffer[1024];
    
    for (auto &wp: perturbations) {
      int p = wp.p->proposal_count();
      int a = wp.p->acceptance_count();
      
      double f = 0.0;
      if (p > 0) {
	f = (double)a/(double)p * 100.0;
      }
      
      sprintf(linebuffer, "  %12s: %6d %6d : %6.2f\n", wp.p->displayname(), p, a, f);

      s += linebuffer;
    }

    return s;
  }
	       
private:

  bool primary()
  {
    return (communicator == MPI_COMM_NULL) || (rank == 0);
  }

  void communicate(int &i)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&i, 1, MPI_INT, 0, communicator);
    }
  }

  void communicate(double &d)
  {
    if (communicator != MPI_COMM_NULL) {
      MPI_Bcast(&d, 1, MPI_DOUBLE, 0, communicator);
    }
  }
  

  struct WeightedPerturbation {
    WeightedPerturbation(PerturbationS2Voronoi<value> *_p, double _w) :
      p(_p),
      w(_w)
    {
    }
    
    PerturbationS2Voronoi<value> *p;
    double w;
  };

  MPI_Comm communicator;
  int rank;
  int size;
  
  double weight_sum;
  std::vector<WeightedPerturbation> perturbations;
  int active;

};

#endif // perturbationcollectionS2Voronoi_hpp
