#pragma once
#ifndef hierarchical_model_hpp
#define hierarchical_model_hpp

#include <ostream>
#include <istream>

class hierarchical_model {
public:

  hierarchical_model(int nhierarchical);
  virtual ~hierarchical_model();

  virtual int get_nhierarchical() const;
  virtual void set(int i, double v);
  virtual double get(int i) const;

  virtual bool write(std::ostream &s) const;
  virtual bool read(std::istream &s);

protected:

  int nhierarchical;
  double *hierarchical;
  
};

class singlescaling_hierarchical_model : public hierarchical_model {
public:

  singlescaling_hierarchical_model(double lambda = 1.0);
  ~singlescaling_hierarchical_model();

};

#endif // hierarchical_model_hpp
  
