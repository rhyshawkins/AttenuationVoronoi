#pragma once
#ifndef chainhistoryVoronoi_hpp
#define chainhistoryVoronoi_hpp

#include <vector>

#include "coordinate.hpp"
#include "hierarchical_model.hpp"
#include "attenuationexception.hpp"
#include "util.hpp"

#include "sphericalvoronoimodel.hpp"

extern "C" {
  #include "slog.h"
};

template
<
  typename coord,
  typename value
>
class deltaVoronoi;

template
<
  typename coord,
  typename value
>
class model_initializationVoronoi;
  
template
<
  typename coord,
  typename value
>
class model_deltaVoronoi;

template
<
  typename coord,
  typename value
>
class hierarchical_deltaVoronoi;

template
<
  typename coord,
  typename value
>
class deltaVoronoi {
public:

  enum {
    DELTA_INITIALIZATION = 0,
    DELTA_DELTA = 1,
    DELTA_HIERARCHICAL = 2
  };

  deltaVoronoi(int _id) :
    id(_id),
    proposed_like(0.0),
    accepted(false)
  {
  }
  
  virtual ~deltaVoronoi()
  {
  }

  virtual int write_header(FILE *fp)
  {
    if (fwrite(&id, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    if (fwrite(&proposed_like, sizeof(double), 1, fp) != 1) {
      return -1;
    }
    
    int a = (int)accepted;
    if (fwrite(&a, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    return 0;
  }

  static int read_header(FILE *fp, double &like, bool &accepted)
  {
    if (fread(&like, sizeof(double), 1, fp) != 1) {
      return -1;
    }
    
    int a;
    if (fread(&a, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    accepted = (bool)a;
    
    return 0;
  }

  virtual int write(FILE *fp) = 0;

  virtual int apply(sphericalvoronoimodel<value> &model, hierarchical_model &hierarchical) = 0;

  virtual void accept()
  {
    accepted = true;
  }
    

  virtual void reject()
  {
    // Nothing to do: assumed rejected on creation
  }
    

  virtual void set_proposed_likelihood(double like)
  {
    proposed_like = like;
  }

  virtual bool isaccepted() const
  {
    return accepted;
  }
  
  virtual double get_proposed_likelihood() const
  {
    return proposed_like;
  }

  static deltaVoronoi *read(FILE *fp)
  {
    int id;
    
    if (fread(&id, sizeof(int), 1, fp) != 1) {
      return nullptr;
    }

    switch (id) {
    case DELTA_INITIALIZATION:
      return model_initializationVoronoi<coord, value>::read(fp);

    case DELTA_DELTA:
      return model_deltaVoronoi<coord, value>::read(fp);

    case DELTA_HIERARCHICAL:
      return hierarchical_deltaVoronoi<coord, value>::read(fp);

    default:
      fprintf(stderr, "deltaVoronoi::read: invalid unique index: %d (%d)\n", id, (int)readers.size());
      return nullptr;
    }
  }

private:

  int id;
  double proposed_like;
  bool accepted;
  
  static std::vector<deltaVoronoi<coord, value>* (*)(FILE *fp)> readers;
  
};

template
<
  typename coord,
  typename value
>
std::vector<deltaVoronoi<coord, value>* (*)(FILE *fp)> deltaVoronoi<coord, value>::readers;


template
<
  typename coord,
  typename value
>
class model_initializationVoronoi : public deltaVoronoi<coord, value> {
public:
  typedef typename sphericalvoronoimodel<value>::cell_t cell_t;
  
  model_initializationVoronoi(sphericalvoronoimodel<value> &_model,
			      hierarchical_model &_hierarchical,
			      double _likelihood) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_INITIALIZATION)
  {
    // Initialization always accepted
    deltaVoronoi<coord, value>::accept();
    
    deltaVoronoi<coord, value>::set_proposed_likelihood(_likelihood);

    int ncells = _model.ncells();
    for (int i = 0; i < ncells; i ++) {
      cell_t *c = _model.get_cell_by_index(i);
      initial_cells.push_back(cellinitialization(c->c, c->v));
    }
    
    int nh = _hierarchical.get_nhierarchical();
    for (int i = 0; i < nh; i ++) {
      hierarchical.push_back(_hierarchical.get(i));
    }
  }
    
  virtual ~model_initializationVoronoi()
  {
  }

  virtual int write(FILE *fp)
  {
    if (deltaVoronoi<coord, value>::write_header(fp) < 0) {
      return -1;
    }
    
    int ncells = initial_cells.size();

    if (fwrite(&ncells, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    for (auto &n : initial_cells) {
      
      if (!n.c.writebinary(fp)) {
	return -1;
      }
      
      if (fwrite(&(n.v), sizeof(double), 1, fp) != 1) {
	return -1;
      }
      
    }

    int nh = (int)hierarchical.size();
    
    if (fwrite(&nh, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    for (auto &h : hierarchical) {
      if (fwrite(&h, sizeof(double), 1, fp) != 1) {
	return -1;
      }
    }
    
    return 0;
  }

  virtual int apply(sphericalvoronoimodel<value> &_model, hierarchical_model &_hierarchical)
  {
    _model.reset();
    
    for (auto &c : initial_cells) {

      _model.add_cell(c.c, c.v);
      
    }

    int i = 0;
    for (auto &h : hierarchical) {
      _hierarchical.set(i, h);
      i ++;
    }
    
    return 0;
  }

  static deltaVoronoi<coord, value> *read(FILE *fp)
  {
    double like;
    bool accepted;
    
    if (deltaVoronoi<coord, value>::read_header(fp, like, accepted) < 0) {
      fprintf(stderr, "model_initialization::read: failed to read header\n");
      return nullptr;
    }
    
    model_initializationVoronoi *r = new model_initializationVoronoi();
    
    r->accept();
    r->set_proposed_likelihood(like);
    
    int ncells;
    
    if (fread(&ncells, sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "model_initialization::read: failed to read no. cells\n");
      return nullptr;
    }
    
    for (int i = 0; i < ncells; i ++) {
      
      coord cd;
      if (!cd.readbinary(fp)) {
	fprintf(stderr, "model_initialization::read: failed to read cell coord\n");
	return nullptr;
      }
      
      double cellvalue;
      if (fread(&cellvalue, sizeof(double), 1, fp) != 1) {
	fprintf(stderr, "model_initialization::read: failed to read cell value\n");
	return nullptr;
      }
      
      r->initial_cells.push_back(cellinitialization(cd, cellvalue));
    }

    int nh;

    if (fread(&nh, sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "model_initialization::read: failed to no. hierarchical\n");
      return nullptr;
    }
    
    for (int i = 0; i < nh; i ++) {
      double h;
      if (fread(&h, sizeof(double), 1, fp) != 1) {
	fprintf(stderr, "model_initialization::read: failed to hierarchical %d\n", i);
	return nullptr;
      }
      
      r->hierarchical.push_back(h);
    }
    
    return r;
  }
  
private:

  model_initializationVoronoi() :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_INITIALIZATION)
  {
  }
  
  struct cellinitialization {

    cellinitialization(const coord &_c, const value &_v) :
      c(_c),
      v(_v)
    {
    }
    
    coord c;
    value v;
  };

  std::vector<cellinitialization> initial_cells;
  std::vector<double> hierarchical;
};

template
<
  typename coord,
  typename value
> 
class model_deltaVoronoi : public deltaVoronoi<coord, value> {
public:

  //
  // Change value
  //
  static model_deltaVoronoi *mkvalue(int cellindex,
				     value oldvalue,
				     value newvalue)
  {
    return new model_deltaVoronoi(VALUE,
				  cellindex,
				  oldvalue,
				  newvalue);
  }

  //
  // Move
  //
  static model_deltaVoronoi *mkmove(int cellindex,
				    coord oldposition,
				    coord newposition)
  {
    return new model_deltaVoronoi(MOVE,
				  cellindex,
				  oldposition,
				  newposition);
  }

  //
  // Birth
  //
  static model_deltaVoronoi *mkbirth(coord newposition,
				     value newvalue)
  {
    return new model_deltaVoronoi(BIRTH,
				  newposition,
				  newvalue);
  }

  static model_deltaVoronoi *mkdeath(int cellindex)
  {
    return new model_deltaVoronoi(DEATH,
				  cellindex);
  }
			      
  ~model_deltaVoronoi()
  {
  }
  
  virtual int write(FILE *fp)
  {
    if (deltaVoronoi<coord, value>::write_header(fp) < 0) {
      return -1;
    }

    int itype = (int)type;
    if (fwrite(&itype, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to write type");
      return -1;
    }

    if (fwrite(&cellindex, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to write cell index\n");
      return -1;
    }

    double v = oldvalue;
    if (fwrite(&v, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old value\n");
      return -1;
    }

    v = newvalue;
    if (fwrite(&v, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to write old value\n");
      return -1;
    }

    if (!oldposition.writebinary(fp)) {
      ERROR("Failed to write old coordinate\n");
      return -1;
    }

    if (!newposition.writebinary(fp)) {
      ERROR("Failed to write new coordinate\n");
      return -1;
    }

    return 0;
  }
  
    
  virtual int apply(sphericalvoronoimodel<value> &model, hierarchical_model &hierarchical)
  {
    if (deltaVoronoi<coord, value>::isaccepted()) {

      switch (type) {
      case BIRTH:
	{
	  model.add_cell(newposition, newvalue);
	}
	break;
	
      case DEATH:
	{
	  model.delete_cell(cellindex);
	}
	break;

      case MOVE:
	{
	  model[cellindex].c = newposition;
	}
	break;

      case VALUE:
	{
	  model[cellindex].v = newvalue;
	}
	break;
	  
      default:
	throw ATTENUATIONEXCEPTION("Unhandled attenuation delta type\n");
      }
    }
    return 0;
  }
  
  static deltaVoronoi<coord, value> *read(FILE *fp)
  {
    double like;
    bool accepted;
    if (deltaVoronoi<coord, value>::read_header(fp, like, accepted) < 0) {
      fprintf(stderr, "model_deltaVoronoi::read: failed to read header\n");
      return nullptr;
    }
    
    model_deltaVoronoi *r = new model_deltaVoronoi();
    
    if (accepted) {
      r->accept();
    }
    r->set_proposed_likelihood(like);

    int itype;
    if (fread(&itype, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to read type");
      return nullptr;
    }
    if (itype < BIRTH || itype > MOVE) {
      ERROR("Invalid type: %d\n", itype);
    }
    r->type = (delta_t)itype;

    if (fread(&r->cellindex, sizeof(int), 1, fp) != 1) {
      ERROR("Failed to read cell index\n");
      return nullptr;
    }

    double v;
    if (fread(&v, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old value\n");
      return nullptr;
    }
    r->oldvalue = v;

    if (fread(&v, sizeof(double), 1, fp) != 1) {
      ERROR("Failed to read old value\n");
      return nullptr;
    }
    r->newvalue = v;

    if (!r->oldposition.readbinary(fp)) {
      ERROR("Failed to read old coordinate\n");
      return nullptr;
    }

    if (!r->newposition.readbinary(fp)) {
      ERROR("Failed to read new coordinate\n");
      return nullptr;
    }
    
    return r;
  }

private:

  typedef enum {
    UNKNOWN = -1,
    BIRTH = 0,
    DEATH = 1,
    VALUE = 2,
    MOVE = 3
  } delta_t;

  model_deltaVoronoi() :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_DELTA)
  {
  }
  
  //
  // Birth Cell
  //
  model_deltaVoronoi(delta_t _type,
		     const coord &_newposition,
		     const value &_newvalue) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_DELTA),
    type(_type),
    cellindex(-1),
    oldvalue(0.0),
    newvalue(_newvalue),
    newposition(_newposition)
  {
  }

  //
  // Death Cell
  model_deltaVoronoi(delta_t _type,
		     int _cellindex) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_DELTA),
    type(_type),
    cellindex(_cellindex),
    oldvalue(0.0),
    newvalue(0.0)
  {
  }
    

  //
  // Value
  //
  model_deltaVoronoi(delta_t _type,
		     int _cellindex,
		     value _oldvalue,
		     value _newvalue) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_DELTA),
    type(_type),
    cellindex(_cellindex),
    oldvalue(_oldvalue),
    newvalue(_newvalue)
  {
  }

  //
  // Move
  //
  model_deltaVoronoi(delta_t _type,
		     int _cellindex,
		     const coord &_oldposition,
		     const coord &_newposition) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_DELTA),
    type(_type),
    cellindex(_cellindex),
    oldvalue(0.0),
    newvalue(0.0),
    oldposition(_oldposition),
    newposition(_newposition)
  {
  }

  delta_t type;

  int cellindex;
  value oldvalue;
  value newvalue;
  coord oldposition;
  coord newposition;
  
};

template
<
  typename coord,
  typename value
>
class hierarchical_deltaVoronoi : public deltaVoronoi<coord, value> {
public:

  static const int DELTAID = deltaVoronoi<coord, value>::register_reader(hierarchical_deltaVoronoi::read);

  hierarchical_deltaVoronoi(int nhierarchical, int *indices, double *old_value, double *new_value) :
    deltaVoronoi<coord, value>(deltaVoronoi<coord, value>::DELTA_HIERARCHICAL)
  {
    for (int i = 0; i < nhierarchical; i ++) {
      
      hierarchical_term_delta d;
      
      d.index = indices[i];
      d.old_value = old_value[i];
      d.new_value = new_value[i];
      
      hierarchical.push_back(d);
    }
  }
  ~hierarchical_deltaVoronoi()
  {
  }
  
  virtual int write(FILE *fp)
  {
    if (deltaVoronoi<coord, value>::write_header(fp) < 0) {
      return -1;
    }
    
    int nh = (int)hierarchical.size();
    if (fwrite(&nh, sizeof(int), 1, fp) != 1) {
      return -1;
    }
    
    for (auto &dh : hierarchical) {
      if (fwrite(&(dh.index), sizeof(int), 1, fp) != 1) {
	return -1;
      }
      if (fwrite(&(dh.old_value), sizeof(double), 1, fp) != 1) {
	return -1;
      }
      if (fwrite(&(dh.new_value), sizeof(double), 1, fp) != 1) {
	return -1;
      }
    }
    
    return 0;
  }

  virtual int apply(sphericalvoronoimodel<value> &model, hierarchical_model &_hierarchical)
  {
    if (deltaVoronoi<coord, value>::isaccepted()) {
      for (auto &dh : hierarchical) {
	
	if (dh.index < 0 || dh.index >= _hierarchical.get_nhierarchical()) {
	  throw ATTENUATIONEXCEPTION("Hierarchical index out of range: %d (%d)\n",
				  dh.index, _hierarchical.get_nhierarchical());
	}
	
	_hierarchical.set(dh.index, dh.new_value);
      }
    }
  
    return 0;
  }
  
  static deltaVoronoi<coord, value> *read(FILE *fp)
  {
    double like;
    bool accepted;
    if (deltaVoronoi<coord, value>::read_header(fp, like, accepted) < 0) {
      return nullptr;
    }
    
    hierarchical_deltaVoronoi *r = new hierarchical_deltaVoronoi(0, nullptr, nullptr, nullptr);
    
    if (accepted) {
      r->accept();
    }
    r->set_proposed_likelihood(like);
    
    int nh;
    if (fread(&nh, sizeof(int), 1, fp) != 1) {
      return nullptr;
    }
    
    for (int i = 0; i < nh; i ++) {
      
      hierarchical_term_delta dh;
      
      if (fread(&(dh.index), sizeof(int), 1, fp) != 1) {
	return nullptr;
      }

      if (fread(&(dh.old_value), sizeof(double), 1, fp) != 1) {
	return nullptr;
      }
      if (fread(&(dh.new_value), sizeof(double), 1, fp) != 1) {
	return nullptr;
      }
      
      r->hierarchical.push_back(dh);
    }
    
    return r;
  }
  
private:
  
  struct hierarchical_term_delta {
    int index;
    double old_value;
    double new_value;
  };
  std::vector<hierarchical_term_delta> hierarchical;
};

template
<
  typename coord,
  typename value
>
class chainhistorywriterVoronoi {
public:

  chainhistorywriterVoronoi(const char *_filename,
			    sphericalvoronoimodel<value> &_initial_model,
			    hierarchical_model &_hierarchical,
			    double _likelihood) :
    fp(fopen(_filename, "w"))
  {
    if (fp == NULL) {
      throw ATTENUATIONEXCEPTION("Failed to create chain history file: %s\n", _filename);
    }
    
    steps.push_back(new model_initializationVoronoi<coord, value>(_initial_model, _hierarchical, _likelihood));
    
    flush();
  }

  ~chainhistorywriterVoronoi()
  {
    flush();
    fclose(fp);
  }


  void add(deltaVoronoi<coord, value> *d)
  {
    steps.push_back(d);
  }

  void flush()
  {
    for (auto &s : steps) {
      s->write(fp);
      delete s;
    }
    
    steps.clear();
  }

private:

  FILE *fp;
  std::vector<deltaVoronoi<coord, value>*> steps;
};

template
<
  typename coord,
  typename value
>

class chainhistoryreaderVoronoi {
public:
  
  chainhistoryreaderVoronoi(const char *filename) :
    fp(fopen(filename, "r")),
    current_likelihood(-1.0)
  {
    if (fp == NULL) {
      throw ATTENUATIONEXCEPTION("Failed to open file for reading: %s\n", filename);
    }
  }
  ~chainhistoryreaderVoronoi()
  {
    fclose(fp);
  }

  int step(sphericalvoronoimodel<value> &model, hierarchical_model &hierarchical, double &likelihood)
  {
    deltaVoronoi<coord, value> *d = deltaVoronoi<coord, value>::read(fp);

    if (d == nullptr) {
      if (feof(fp)) {
	return 0;
      } else {
	fprintf(stderr, "chainhistoryreaderVoronoi::step: failed to read next step\n");
	return -1;
      }
    }
    
    if (d->apply(model, hierarchical) < 0) {
      fprintf(stderr, "chainhistoryreaderVoronoi::step: failed to apply step to model/hierarchical/likelihood\n");
      return -1;
    }
    
    if (d->isaccepted()) {
      current_likelihood = d->get_proposed_likelihood();
    }
    
    likelihood = current_likelihood;
    
    delete d;
    
    return 1;
  }
  
private:

  FILE *fp;
  double current_likelihood;
  
};

#endif // chainhistoryVoronoi_hpp
  
