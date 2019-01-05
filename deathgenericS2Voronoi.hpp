#pragma once
#ifndef deathgenericS2Voronoi_hpp
#define deathgenericS2Voronoi_hpp

#include "globalS2Voronoi.hpp"
#include "perturbationS2Voronoi.hpp"

template
<typename value>
class DeathGenericS2Voronoi : public PerturbationS2Voronoi<value> {
public:

  typedef sphericalcoordinate<value> coord_t;
  typedef deltaVoronoi<coord_t, value> delta_t;
  typedef model_deltaVoronoi<coord_t, value> model_delta_t;
  typedef typename sphericalvoronoimodel<value>::cell_t cell_t;
  
  DeathGenericS2Voronoi(Proposal *_value_proposal,
			SphericalProposal *_position_proposal) :
    value_proposal(_value_proposal),
    position_proposal(_position_proposal),
    undo_index(-1),
    p(0),
    a(0)
  {
  }

  ~DeathGenericS2Voronoi()
  {
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
		       delta_t *&perturbation)
  {
    bool validproposal = false;
    int cell;

    if (this->primary()) {
      
      p ++;
    
      int k = model.ncells();
      if (k > 1) {
	
	cell = random.uniform(k);
	validproposal = true;
	perturbation = model_deltaVoronoi<coord_t, value>::mkdeath(cell);

      } else {
	perturbation = model_deltaVoronoi<coord_t, value>::mkdeath(-1);
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      this->communicate(cell);

      //
      // Store undo information
      //
      cell_t *c = model.get_cell_by_index(cell);

      undo_index = cell;
      undo_coord = c->c;
      undo_value = c->v;

      //
      // Compute ratios
      //
      last_log_proposal_ratio = 0.0;
      log_prior_ratio = -(prior.logpdf(undo_value) + position_prior.logpdf(undo_coord.phi,
									   undo_coord.theta));
      model.delete_cell(cell);
      
      value new_value = model.value_at_point(undo_coord);
      
      last_log_proposal_ratio =
	position_proposal->log_proposal(random,
					temperature,
					0.0,
					0.0,
					undo_coord.phi,
					undo_coord.theta) +
	value_proposal->log_proposal(random,
				     temperature,
				     new_value,
				     undo_value);
    }

    return validproposal;
  }

  virtual double log_proposal_ratio(Rng &random,
				    PriorProposal &prior,
				    SphericalPriorProposal &position_prior,
				    sphericalvoronoimodel<value> &proposed_model,
				    PriorProposal &hierarchical_prior,
				    hierarchical_model &proposed_hierarchical,
				    double temperature)
  {
    return last_log_proposal_ratio;
  }
  
  void accept()
  {
    a ++;
    
    if (undo_index < 0) {
      throw ATTENUATIONEXCEPTION("No undo information\n");
    }
    
    undo_index = -1;
  }

  void reject(sphericalvoronoimodel<value> &model)
  {
    if (undo_index < 0) {
      throw ATTENUATIONEXCEPTION("No undo information\n");
    }

    //
    // Re-insert the deleted node
    //
    model.insert_cell(undo_index, undo_coord, undo_value);

    undo_index = -1;
  }

  virtual int proposal_count() const
  {
    return p;
  }
  
  virtual int acceptance_count() const
  {
    return a;
  }

  virtual const char *displayname() const
  {
    return "Death";
  }

private:

  Proposal *value_proposal;
  SphericalProposal *position_proposal;
  
  int undo_index;
  value undo_value;
  coord_t undo_coord;

  double last_log_proposal_ratio;

  int p;
  int a;

};

#endif // deathgenericS2Voronoi_hpp
