#pragma once
#ifndef moveS2Voronoi_hpp
#define moveS2Voronoi_hpp

#include "globalS2Voronoi.hpp"
#include "perturbationS2Voronoi.hpp"

template
<typename value>
class MoveS2Voronoi : public PerturbationS2Voronoi<value> {
public:

  typedef sphericalcoordinate<value> coord_t;
  typedef typename sphericalvoronoimodel<value>::cell_t cell_t;
  typedef deltaVoronoi<coord_t, value> delta_t;
  typedef model_deltaVoronoi<coord_t, value> model_delta_t;

  MoveS2Voronoi() :
    undo_cell(nullptr),
    last_log_proposal_ratio(0.0),
    p(0),
    a(0)
  {
  }
  
  ~MoveS2Voronoi()
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
    coord_t newposition;


    if (this->primary()) {
      
      p ++;

      cell = random.uniform(model.ncells());
    
      //
      // Next get the active node
      //
      cell_t *c = model.get_cell_by_index(cell);
      
      value newphi;
      value newtheta;
      if (position_prior.propose(random,
				 temperature,
				 c->c.phi,
				 c->c.theta,
				 newphi,
				 newtheta,
				 log_prior_ratio)) {

	validproposal = true;
	newposition = coord_t(newphi, newtheta);
	perturbation = model_deltaVoronoi<coord_t, value>::mkmove(cell, c->c, newposition);
	
      } else {
	
	perturbation = model_deltaVoronoi<coord_t, value>::mkmove(cell, coord_t(), coord_t());

      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      this->communicate(cell);
      this->communicate(newposition);

      cell_t *c = model.get_cell_by_index(cell);
      
      undo_cell = c;
      undo_coord = c->c;
      c->c = newposition;

      last_log_proposal_ratio =
	position_prior.log_proposal_ratio(random,
					  temperature,
					  undo_coord.phi,
					  undo_coord.theta,
					  newposition.phi,
					  newposition.theta);
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
    
    if (undo_cell == nullptr) {
      throw ATTENUATIONEXCEPTION("No undo information\n");
    }
    
    undo_cell = nullptr;
  }
  
  void reject(sphericalvoronoimodel<value> &model)
  {
    if (undo_cell == nullptr) {
      throw ATTENUATIONEXCEPTION("No undo information\n");
    }

    undo_cell->c = undo_coord;
    
    undo_cell = nullptr;
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
    return "Move";
  }
  
private:
  
  cell_t *undo_cell;
  coord_t undo_coord;

  double last_log_proposal_ratio;
  
  int p;
  int a;
  
};

#endif // moveS2Voronoi_hpp
