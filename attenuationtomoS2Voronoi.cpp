#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "attenuationdataS2.hpp"

#include "globalS2Voronoi.hpp"

#include "perturbationcollectionS2Voronoi.hpp"
#include "valueS2Voronoi.hpp"
#include "birthgenericS2Voronoi.hpp"
#include "deathgenericS2Voronoi.hpp"
#include "moveS2Voronoi.hpp"
#include "hierarchicalS2Voronoi.hpp"

#include "pathutil.hpp"

typedef chainhistorywriterVoronoi<sphericalcoordinate<double>, double> chainhistorywriter_t;

static char short_options[] = "i:I:o:P:H:M:B:D:T:S:t:l:v:b:pLh";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"initial", required_argument, 0, 'I'},
  {"output", required_argument, 0, 'o'},

  {"prior", required_argument, 0, 'P'},
  {"hierarchical-prior", required_argument, 0, 'H'},
  {"move-prior", required_argument, 0, 'M'},
  {"birth-death-prior", required_argument, 0, 'B'},

  {"max-cells", required_argument, 0, 'T'},

  {"total", required_argument, 0, 't'},
  {"lambda", required_argument, 0, 'l'},
  
  {"verbosity", required_argument, 0, 'v'},

  {"birth-probability", required_argument, 0, 'b'},
  {"posterior", no_argument, 0, 'p'},
  {"logspace", no_argument, 0, 'L'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;

  //
  // Configuration
  //
  char *input;
  char *initial;
  char *output;

  char *prior;
  char *hierarchicalprior;
  char *positionprior;
  char *birthdeathprior;

  int maxcells;

  double lambda;

  int seed;

  bool posterior;

  int total;
  int verbosity;

  double Pb;

  bool logspace;
  
  //
  // State
  //
  globalS2Voronoi<double> *global;
  
  double current_likelihood;

  //
  // Misc
  //
  char filename[1024];

  //
  // Initialize defaults
  //
  input = nullptr;
  initial = nullptr;
  output = nullptr;

  prior = nullptr;
  hierarchicalprior = nullptr;
  positionprior = nullptr;
  birthdeathprior = nullptr;

  maxcells = 1000;

  lambda = 1.0;

  seed = 983;

  posterior = false;

  total = 1000;
  verbosity = 1000;

  Pb = 0.05;

  logspace = false;
  
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

    case 'I':
      initial = optarg;
      break;

    case 'o':
      output = optarg;
      break;

    case 'P':
      prior = optarg;
      break;

    case 'H':
      hierarchicalprior = optarg;
      break;

    case 'M':
      positionprior = optarg;
      break;

    case 'B':
      birthdeathprior = optarg;
      break;

    case 'T':
      maxcells = atoi(optarg);
      if (maxcells < 1) {
	fprintf(stderr, "error: max cells must be 1 or greater\n");
	return -1;
      }
      break;

    case 't':
      total = atoi(optarg);
      if (total < 1) {
	fprintf(stderr, "error: total must be greater than 0\n");
	return -1;
      }
      break;

    case 'l':
      lambda = atof(optarg);
      if (lambda <= 0.0) {
	fprintf(stderr, "error: lambda must be greater than 0\n");
	return -1;
      }
      break;

    case 'v':
      verbosity = atoi(optarg);
      if (verbosity < 0) {
	fprintf(stderr, "error: verbosity must be greater than 0\n");
	return -1;
      }
      break;

    case 'b':
      Pb = atof(optarg);
      if (Pb < 0.0 || Pb >= 0.5) {
	fprintf(stderr, "error: Pb must be between 0 and 0.5\n");
	return -1;
      }
      break;
      
    case 'p':
      posterior = true;
      break;

    case 'L':
      logspace = true;
      break;
	
      
    case 'h':
    default:
      usage(argv[0]);
      return -1;

    }
  }
  
  if (input == nullptr) {
    fprintf(stderr, "error: required parameter input not set\n");
    return -1;
  }

  global = new globalS2Voronoi<double>(input,
				       initial,
				       prior,
				       hierarchicalprior,
				       positionprior,
				       birthdeathprior,
				       maxcells,
				       lambda,
				       1.0, // temperature
				       seed,
				       posterior,
				       logspace);
  
  current_likelihood = global->likelihood();
  printf("Initial likelihood: %10.6f\n", current_likelihood);
  global->accept();

  int *khistogram = new int[global->maxcells + 1];
  for (int i = 0; i <= global->maxcells; i ++) {
    khistogram[i] = 0;
  }

  PerturbationCollectionS2Voronoi<double> pc;
  mkpath(output, "ch.dat", filename);
  chainhistorywriter_t *history = new chainhistorywriter_t(filename,
							   *(global->model),
							   *(global->hierarchical),
							   current_likelihood);

  if (posterior) {
    pc.add(new ValueS2Voronoi<double>(), 0.1);

    pc.add(new BirthGenericS2Voronoi<double>(global->birthdeathvalueproposal,
					     global->birthdeathpositionproposal), 1.0);
    pc.add(new DeathGenericS2Voronoi<double>(global->birthdeathvalueproposal,
					     global->birthdeathpositionproposal), 1.0);
  } else {
    pc.add(new ValueS2Voronoi<double>(), 1.0);

    if (Pb > 0.0) {
      pc.add(new MoveS2Voronoi<double>(), 0.5);
      
      pc.add(new BirthGenericS2Voronoi<double>(global->birthdeathvalueproposal,
					       global->birthdeathpositionproposal), Pb);
      pc.add(new DeathGenericS2Voronoi<double>(global->birthdeathvalueproposal,
					       global->birthdeathpositionproposal), Pb);
    }

    if (hierarchicalprior != nullptr) {
      pc.add(new HierarchicalS2Voronoi<double>(), 0.5);
    }
  }

  for (int i = 0; i < total; i ++) {
    
    double log_prior_ratio;
    double log_proposal_ratio;
    PerturbationS2Voronoi<double>::delta_t *perturbation = nullptr;
    
    if (pc.propose(*global, log_prior_ratio, perturbation)) {
      
      if (perturbation == nullptr) {
	throw ATTENUATIONEXCEPTION("Valid proposal has null perturbation\n");
      }

      double u = log(global->random.uniform());

      double proposed_likelihood = global->likelihood();
      perturbation->set_proposed_likelihood(proposed_likelihood);

      log_proposal_ratio = pc.log_proposal_ratio(*global);

      if (u < (current_likelihood - proposed_likelihood + log_prior_ratio + log_proposal_ratio)) {

	pc.accept(*global);
	perturbation->accept();
	global->accept();
	
	current_likelihood = proposed_likelihood;
      } else {
	pc.reject(*global);
	perturbation->reject(); 
	global->reject();
      }
    }

    if (verbosity > 0 && (i + 1) % verbosity == 0) {

      printf("%5d: Cells %d Likelihood %10.6f Lambda %10.6f\n",
             i + 1,
             (int)global->model->ncells(),
             current_likelihood,
             global->hierarchical->get(0));

      pc.writeacceptancereport(stdout);
    }
      
    int k = global->model->ncells();
    if (k < 1 || k > maxcells) {
      throw ATTENUATIONEXCEPTION("k out of range: %d (%d)\n", k, maxcells);
    }
    khistogram[k] ++;
      
    history->add(perturbation);
  }

  //
  // Save khistogram
  //
  mkpath(output, "khistogram.txt", filename);
  FILE *fp = fopen(filename, "w");
  for (int i = 0; i <= global->maxcells; i ++) {
    fprintf(fp, "%d %d\n", i, khistogram[i]);
  }
  fclose(fp);
  delete [] khistogram;

  //
  // Save residuals
  //
  mkpath(output, "residuals.txt", filename);
  fp = fopen(filename, "w");
  for (int i = 0; i < global->residual_size; i ++) {
    fprintf(fp, "%.9g\n", global->mean_residuals[i]);
  }
  fclose(fp);
  
  history->flush();

  delete history;

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -i|--input <filename>                   Observations input file\n"
	  " -I|--initial <filename>                 Initial model input file\n"
	  " -o|--output <path>                      Path/prefix for output files\n"
	  "\n"
	  " -P|--prior <filename>                   Prior input file\n"
	  " -H|--hierarchical-prior <filename>      Hierarchical prior input file\n"
	  " -M|--move-prior <filename>              Move prior input file\n"
	  " -B|--birth-death-prior <filename>       Birth/Death proposal file\n"
	  "\n"
	  " -t|--total <int>                        Total number of iterations\n"
	  " -v|--verbosity <int>                    Number of iterations between status updates (0 = none)\n"
	  "\n"
	  " -l|--lambda <float>                     Initial/fixed lambda parameter\n"
	  "\n"
	  " -L|--logspace                           Model is in log(Q)\n"
	  " -T|--max-cells <int>                    Max no. Voronoi cells\n"
	  " -b|--birth-probability <float>          Relative probability of birth\n"
	  " -p|--posterior                          Posterior test\n"
	  "\n"
	  " -h|--help                               Usage information\n"
	  "\n",
	  pname);

}
