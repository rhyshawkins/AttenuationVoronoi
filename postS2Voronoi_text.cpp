
#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "coordinate.hpp"
#include "sphericalvoronoimodel.hpp"
#include "chainhistoryVoronoi.hpp"

typedef sphericalcoordinate<double> coord_t;
typedef chainhistoryreaderVoronoi<coord_t, double> chainhistoryreader_t;

static char short_options[] = "i:o:t:s:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  
  {"thin", required_argument, 0, 't'},
  {"skip", required_argument, 0, 's'},

  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  //
  // Chain processing
  //
  int skip;
  int thin;

  //
  // Input files
  //
  char *input;

  //
  // Output Files
  //
  char *output;
  
  //
  // Defaults
  //

  input = nullptr;
  output = nullptr;
  
  skip = 0;
  thin = 0;
  
  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      input = optarg;
      break;

    case 'o':
      output = optarg;
      break;

    default:
      fprintf(stderr, "error: invalid option '%c'\n", c);
      
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (input == nullptr) {
    fprintf(stderr, "error: required input file parameter missing\n");
    return -1;
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  chainhistoryreader_t reader(input);
    
  sphericalvoronoimodel<double> model(false);
  singlescaling_hierarchical_model hierarchical_model;
  double likelihood;
    
  int status = reader.step(model, hierarchical_model, likelihood);
  int step = 0;

  FILE *fp = fopen(output, "w");
  if (fp == NULL) {
    fprintf(stderr, "error: failed to create text file\n");
    return -1;
  }

  while (status > 0) {
      
    if (step >= skip) {
	
      if (thin <= 1 || (step - skip) % thin == 0) {
	
	fprintf(fp, "%6d %15.9f %2d",
		step, likelihood, hierarchical_model.get_nhierarchical());

	for (int i = 0; i < hierarchical_model.get_nhierarchical(); i ++) {
	  fprintf(fp, "%15.9f ", hierarchical_model.get(i));
	}

	fprintf(fp, "%2d ", model.ncells());

	for (int i = 0; i < model.ncells(); i ++) {
	  typename sphericalvoronoimodel<double>::cell_t &c = model[i];
	
	  fprintf(fp, "%15.9f %15.9f %15.9f ",
		  c.c.phi, c.c.theta, c.v);
	}

	fprintf(fp, "\n");
      }
    }

    status = reader.step(model, hierarchical_model, likelihood);
    step ++;
    
    if ((step % 100000) == 0) {
      printf("%d\n", step);
    }
  }
  
  if (status < 0) {
    fprintf(stderr, "error: failed to step through chain history\n");
    return -1;
  }

  fclose(fp);

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>        Input observations file (required)\n"
	  " -o|--output <filename>       Output text file (required)\n"
	  "\n"
	  " -t|--thin <int>             Only use every nth model\n"
	  " -s|--skip <int>             Skip first n models\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}
