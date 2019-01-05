#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include <mpi.h>

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

static char short_options[] = "i:I:o:P:H:M:B:D:T:S:t:l:v:b:pLc:K:h";
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

  {"chains", required_argument, 0, 'c'},
  {"temperatures", required_argument, 0, 'K'},
  
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

  int seed_base;
  int seed_mult;

  bool posterior;

  int total;
  int verbosity;

  double Pb;
  double Pm;
  bool logspace;

  int chains;
  int temperatures;
  double max_temperature;
  
  //
  // State
  //
  globalS2Voronoi<double> *global;
  
  double current_likelihood;

  //
  // Misc
  //
  char filename[1024];

  int mpi_size;
  int mpi_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

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

  seed_base = 983;
  seed_mult = 101;

  posterior = false;

  total = 1000;
  verbosity = 1000;

  Pb = 0.05;
  Pm = 0.25;

  logspace = false;

  chains = 1;
  temperatures = 1;
  max_temperature = 1000.0;
  
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

    case 'c':
      chains = atoi(optarg);
      if (chains <= 0) {
	fprintf(stderr, "error: no. chains must be greater than 0\n");
	return -1;
      }
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

  if (mpi_size % chains != 0) {
    fprintf(stderr, "error: no. chains (%d) must be a divisor of MPI processes (%d)\n",
	    chains,
	    mpi_size);
    return -1;
  }
  
  mkrankpath(mpi_rank, output, "log.txt", filename);
  if (slog_set_output_file(filename,
                           SLOG_FLAGS_CLEAR) < 0) {
    fprintf(stderr, "error: failed to redirect log file\n");
    return -1;
  }

  MPI_Comm chain_communicator;
  double temperature;
  int chain_rank;
  int chain_id;
  
  if (chains == 1) {
    chain_id = 0;
    temperature = 1.0;
    chain_communicator = MPI_COMM_WORLD;
  } else {
    int processesperchain = mpi_size/chains;
    chain_id = mpi_rank/processesperchain;

    if (temperatures > 1) {
      int temperature_id = chain_id % temperatures;
      temperature = pow(10.0, log10(max_temperature) * (double)temperature_id/(double)(temperatures - 1));
    } else {
      temperature = 1.0;
    }
    MPI_Comm_split(MPI_COMM_WORLD, chain_id, mpi_rank, &chain_communicator);
  }

  MPI_Comm_rank(chain_communicator, &chain_rank);

  const char *initial_model_ptr = nullptr;
  char initial_model_filename[1024];
  if (initial != nullptr) {
    mkrankpath(chain_id, initial, "finalmodel.txt", initial_model_filename);
    initial_model_ptr = initial_model_filename;
  }
					    
  global = new globalS2Voronoi<double>(input,
				       initial_model_ptr,
				       prior,
				       hierarchicalprior,
				       positionprior,
				       birthdeathprior,
				       maxcells,
				       lambda,
				       temperature,
				       seed_base + seed_mult * mpi_rank,
				       posterior,
				       logspace);

  ValueS2Voronoi<double> *value = new ValueS2Voronoi<double>();
  MoveS2Voronoi<double> *move = new MoveS2Voronoi<double>();
  BirthGenericS2Voronoi<double> *birth = new BirthGenericS2Voronoi<double>(global->birthdeathvalueproposal,
									   global->birthdeathpositionproposal);
  DeathGenericS2Voronoi<double> *death = new DeathGenericS2Voronoi<double>(global->birthdeathvalueproposal,
									   global->birthdeathpositionproposal);
  
  global->initialize_mpi(chain_communicator, temperature);
  value->initialize_mpi(chain_communicator);
  move->initialize_mpi(chain_communicator);
  birth->initialize_mpi(chain_communicator);
  death->initialize_mpi(chain_communicator);

  current_likelihood = global->likelihood();
  if (chain_rank == 0) {
    INFO("Chain %03d: Initial likelihood: %10.6f\n", chain_id, current_likelihood);
  }
  global->accept();

  int *khistogram = nullptr;
  chainhistorywriter_t *history = nullptr;
  
  if (chain_rank == 0) {
    khistogram = new int[global->maxcells + 1];
    for (int i = 0; i <= global->maxcells; i ++) {
      khistogram[i] = 0;
    }

    mkrankpath(chain_id, output, "ch.dat", filename);

    history = new chainhistorywriter_t(filename,
				       *(global->model),
				       *(global->hierarchical),
				       current_likelihood);
  }

  //
  //
  //
  PerturbationCollectionS2Voronoi<double> pc;
  
  if (posterior) {
    Pb = 0.45;
  }
  
  pc.add(value, (1.0 - (2.0 * Pb)) * (1.0 - Pm));
  pc.add(move, (1.0 - (2.0 * Pb)) * Pm);
  
  pc.add(birth, Pb);
  pc.add(death, Pb);

  pc.initialize_mpi(chain_communicator);

  if (hierarchicalprior != nullptr) {
    //
    // Add Hierarchical Sampling if required
    //
    HierarchicalS2Voronoi<double> *hierarchical = new HierarchicalS2Voronoi<double>();

    hierarchical->initialize_mpi(chain_communicator);
    
    pc.add(hierarchical, 0.5);
  }

  if (temperatures > 1) {
    //
    // Add PT Exchanges if required
    //
  }
  
  for (int i = 0; i < total; i ++) {
    
    double log_prior_ratio;
    double log_proposal_ratio;
    PerturbationS2Voronoi<double>::delta_t *perturbation = nullptr;
    bool accepted;
    
    if (pc.propose(*global, log_prior_ratio, perturbation)) {

      double proposed_likelihood = global->likelihood();
      
      if (chain_rank == 0) {
	
	if (perturbation == nullptr) {
	  throw ATTENUATIONEXCEPTION("Valid proposal has null perturbation\n");
	}
	
	double u = log(global->random.uniform());
	
	perturbation->set_proposed_likelihood(proposed_likelihood);
	
	log_proposal_ratio = pc.log_proposal_ratio(*global);
	
	accepted = u < (current_likelihood - proposed_likelihood + log_prior_ratio + log_proposal_ratio);

	if (accepted) {
	  perturbation->accept();
	} else {
	  perturbation->reject();
	}
      }

      int t = accepted;
      MPI_Bcast(&t, 1, MPI_INT, 0, chain_communicator);
      accepted = t;

      if (accepted) {
	pc.accept(*global);
	current_likelihood = proposed_likelihood;
	global->accept();
      } else {
	pc.reject(*global);
	global->reject();
      }
    }

    if (chain_rank == 0) {
      if (verbosity > 0 && (i + 1) % verbosity == 0) {
	
	INFO("%5d: Cells %d Likelihood %10.6f Lambda %10.6f\n",
	     i + 1,
	     global->model->ncells(),
	     current_likelihood,
	     global->hierarchical->get(0));

	std::string report = pc.generateacceptancereport();
	INFO("%s", report.c_str());
      }
      
      int k = global->model->ncells();
      if (k < 1 || k > maxcells) {
	throw ATTENUATIONEXCEPTION("k out of range: %d (%d)\n", k, maxcells);
      }
      khistogram[k] ++;
      
      history->add(perturbation);
    }
  }

  if (chain_rank == 0) {
    //
    // Save khistogram
    //
    mkrankpath(chain_id, output, "khistogram.txt", filename);
    FILE *fp = fopen(filename, "w");
    for (int i = 0; i <= global->maxcells; i ++) {
      fprintf(fp, "%d %d\n", i, khistogram[i]);
    }
    fclose(fp);
    delete [] khistogram;
    
    history->flush();
    
    delete history;

    mkrankpath(chain_id, output, "finalmodel.txt", filename);
    if (!global->model->save(filename)) {
      throw ATTENUATIONEXCEPTION("Failed to save final model\n");
    }

    //
    // Save residuals
    //
    mkrankpath(chain_id, output, "residuals.txt", filename);
    fp = fopen(filename, "w");
    for (int i = 0; i < global->residual_size; i ++) {
      fprintf(fp, "%.9g\n", global->mean_residuals[i]);
    }
    fclose(fp);
      
  }

  MPI_Finalize();

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
	  " -T|--max-cells <int>                    Max no. Voronoi cells\n"
	  " -b|--birth-probability <float>          Relative probability of birth\n"
	  " -p|--posterior                          Posterior test\n"
	  "\n"
	  " -c|--chains <int>                       No. of chains to run\n"
	  " -K|--temperatures <int>                 No. of temperatures to run\n"
	  "\n"
	  " -h|--help                               Usage information\n"
	  "\n",
	  pname);

}
