
#include "hierarchical_model.hpp"

#include "attenuationexception.hpp"

hierarchical_model::hierarchical_model(int _nhierarchical) :
  nhierarchical(_nhierarchical),
  hierarchical(new double[_nhierarchical])
{
}

hierarchical_model::~hierarchical_model()
{
  delete [] hierarchical;
}

int
hierarchical_model::get_nhierarchical() const
{
  return nhierarchical;
}

void
hierarchical_model::set(int i, double v)
{
  if (i < 0 || i >= nhierarchical) {
    throw ATTENUATIONEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  hierarchical[i] = v;
}

double
hierarchical_model::get(int i) const
{
  if (i < 0 || i >= nhierarchical) {
    throw ATTENUATIONEXCEPTION("Index out of range %d (%d)\n", i, nhierarchical);
  }

  return hierarchical[i];
}

bool
hierarchical_model::write(std::ostream &s) const
{
  s.write((char*)&nhierarchical, sizeof(int));
  
  for (int i = 0; i < nhierarchical; i ++) {
    s.write((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

bool
hierarchical_model::read(std::istream &s)
{
  int n;

  s.read((char*)&n, sizeof(int));
  if (n != nhierarchical) {
    return false;
  }

  for (int i = 0; i < nhierarchical; i ++) {
    s.read((char*)(&hierarchical[i]), sizeof(double));
  }

  return true;
}

singlescaling_hierarchical_model::singlescaling_hierarchical_model(double lambda) :
  hierarchical_model(1)
{
  set(0, lambda);
}

singlescaling_hierarchical_model::~singlescaling_hierarchical_model()
{
}
