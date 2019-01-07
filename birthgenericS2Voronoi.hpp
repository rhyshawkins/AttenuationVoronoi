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
#ifndef birthgenerics2Voronoi_hpp
#define birthgenerics2Voronoi_hpp

#include "globalS2Voronoi.hpp"
#include "perturbationS2Voronoi.hpp"

template
<typename value>
class BirthGenericS2Voronoi : public PerturbationS2Voronoi<value> {
public:

  typedef sphericalcoordinate<value> coord_t;
  typedef deltaVoronoi<coord_t, value> delta_t;
  typedef model_deltaVoronoi<coord_t, value> model_delta_t;

  BirthGenericS2Voronoi(Proposal *_value_proposal,
			SphericalProposal *_position_proposal) :
    value_proposal(_value_proposal),
    position_proposal(_position_proposal),
    undo_available(false),
    p(0),
    a(0),
    last_log_proposal_ratio(0.0)
  {
  }

  ~BirthGenericS2Voronoi()
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
    value oldvalue;
    value newvalue;
    coord_t newposition;

    if (undo_available) {
      throw ATTENUATIONEXCEPTION("Proposal in progress\n");
    }
    
    if (this->primary()) {
      p ++;
      
      int k = model.ncells();
    
      if (k == maxcells) {
	perturbation = model_deltaVoronoi<coord_t, value>::mkbirth(coord_t(),
								   0.0);
	validproposal = false;
      } else {
	  

	//
	// Generate new point
	//
	value newphi;
	value newtheta;
	value logpriorratio;
	if (!position_proposal->propose(random,
					temperature,
					0.0,
					0.0,
					newphi,
					newtheta,
					logpriorratio)) {
	  perturbation = model_deltaVoronoi<coord_t, value>::mkbirth(coord_t(),
								     0.0);
	  validproposal = false;
	} else {
	  

	  //
	  // Get existing value at point
	  //
	  newposition = coord_t(newphi, newtheta);
	  oldvalue = model.value_at_point(newposition);
	  
	  //
	  // Propose new value
	  //
	  if (!value_proposal->propose(random,
				       temperature,
				       oldvalue,
				       newvalue,
				       logpriorratio)) {
	    perturbation = model_deltaVoronoi<coord_t, value>::mkbirth(newposition,
								       0.0);
	    validproposal = false;
	  } else {

	    

	    perturbation = model_deltaVoronoi<coord_t, value>::mkbirth(newposition,
								       newvalue);
	    
	    validproposal = true;

	  }
	}
      }
    }

    this->communicate(validproposal);

    if (validproposal) {

      this->communicate(newposition);
      this->communicate(oldvalue);
      this->communicate(newvalue);

      log_prior_ratio =
	position_prior.logpdf(newposition.phi, newposition.theta) +
	prior.logpdf(newvalue);

      last_log_proposal_ratio =
	-position_proposal->log_proposal(random,
					 temperature,
					 0.0,
					 0.0,
					 newposition.phi,
					 newposition.theta)
	-value_proposal->log_proposal(random,
				      temperature,
				      oldvalue,
				      newvalue);      
      
      //
      // Proposal valid
      //
      model.add_cell(newposition, newvalue);

      undo_available = true;
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
    if (!undo_available) {
      throw ATTENUATIONEXCEPTION("No proposal in progress\n");
    }
    
    a ++;
    undo_available = false;
  }

  void reject(sphericalvoronoimodel<value> &model)
  {
    if (!undo_available) {
      throw ATTENUATIONEXCEPTION("No proposal in progress\n");
    }

    model.pop();
    
    undo_available = false;
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
    return "Birth";
  }

private:

  Proposal *value_proposal;
  SphericalProposal *position_proposal;

  bool undo_available;
  
  int p;
  int a;

  double last_log_proposal_ratio;
  
};

#endif // birthgenerics2_hpp
