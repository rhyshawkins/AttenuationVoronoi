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

#include <math.h>

#include <map>
#include <string>

#include <getopt.h>

#include "attenuationdataS2.hpp"
#include "coordinate.hpp"
#include "rng.hpp"

double synthetic_constant(double phi, double theta)
{
  return 250.0;
}

double synthetic_eastwest(double phi, double theta)
{
  if (theta < 0.0) {
    return 100.0;
  } else {
    return 400.0;
  }
}

double synthetic_northsouth(double phi, double theta)
{
  if (phi < M_PI/2.0) {
    return 100.0;
  } else {
    return 400.0;
  }
}

double synthetic_cubedsphere(double phi, double theta)
{
  vector3<double> v;

  sphericalcoordinate<double>::sphericaltocartesian(phi, theta, v);

  if (fabs(v.x) > fabs(v.y)) {
    if (fabs(v.x) > fabs(v.z)) {
      if (v.x < 0.0) {
	return 200.0;
      } else {
	return 500.0;
      }
    } else {
      if (v.z < 0.0) {
	return 400.0;
      } else {
	return 300.0;
      }
    }
  } else if (fabs(v.y) > fabs(v.z)) {
    if (v.y < 0.0) {
      return 100.0;
    } else {
      return 600.0;
    }
  } else {

    if (v.z < 0.0) {
      return 400.0;
    } else {
      return 300.0;
    }
  }
}

std::map<std::string, synthetic_model_f> models = { {"Constant", synthetic_constant},
						    {"EastWest", synthetic_eastwest},
						    {"NorthSouth", synthetic_northsouth},
						    {"CubedSphere", synthetic_cubedsphere} };

static char short_options[] = "i:o:O:I:m:ln:S:W:H:h";
static struct option long_options[] = {
  {"input", required_argument, 0, 'i'},
  {"output", required_argument, 0, 'o'},
  {"output-true", required_argument, 0, 'O'},

  {"model", required_argument, 0, 'm'},
  {"list-models", no_argument, 0, 'l'},
  
  {"noise", required_argument, 0, 'n'},
  {"seed", required_argument, 0, 'S'},

  {"image-output", required_argument, 0, 'I'},
  {"image-width", required_argument, 0, 'W'},
  {"image-height", required_argument, 0, 'H'},
  
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  char *source_paths;
  char *output_file;
  char *output_true;
  int seed;
  double noise_sigma;
  std::string model_name;

  char *image_output;
  int image_width;
  int image_height;
    
  //
  // Defaults
  //
  model_name = "Constant";
  source_paths = nullptr;
  output_file = nullptr;
  output_true = nullptr;
  noise_sigma = 0.1;
  seed = 983;

  image_output = nullptr;
  image_width = 128;
  image_height = 64;

  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'i':
      source_paths = optarg;
      break;

    case 'o':
      output_file = optarg;
      break;

    case 'O':
      output_true = optarg;
      break;
      
    case 'I':
      image_output = optarg;
      break;

    case 'm':
      model_name = optarg;
      break;

    case 'l':
      {
	fprintf(stderr, "Models:\n");
	for (auto &mn : models) {
	  fprintf(stderr, "  `%s'\n", mn.first.c_str());
	}
	return -1;
      }
      break;

    case 'n':
      noise_sigma = atof(optarg);
      if (noise_sigma <= 0.0) {
	fprintf(stderr, "error: noise must be greater than 0\n");
	return -1;
      }
      break;

    case 'S':
      seed = atoi(optarg);
      break;

    case 'W':
      image_width = atoi(optarg);
      if (image_width < 16) {
	fprintf(stderr, "error: image width must be 16 or greater\n");
	return -1;
      }
      break;

    case 'H':
      image_height = atoi(optarg);
      if (image_height < 16) {
	fprintf(stderr, "error: image_height must be 16 or greater\n");
	return -1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return -1;

    }
  }
  
  if (source_paths == nullptr) {
    fprintf(stderr, "error: required parameter input not set\n");
    return -1;
  }

  if (output_file == nullptr) {
    fprintf(stderr, "error: required parameter output not set\n");
    return -1;
  }
      
  attenuationdataS2<double> data(source_paths);

  std::map<std::string, synthetic_model_f>::const_iterator i = models.find(model_name);

  if (i == models.end()) {
    fprintf(stderr, "error: invalid model name: %s\n", model_name.c_str());
    return -1;
  }

  synthetic_model_f model = i->second;

  FILE *fp_out = fopen(output_file, "w");
  if (fp_out == NULL) {
    fprintf(stderr, "error: failed to create output file\n");
    return -1;
  }

  FILE *fp_true = nullptr;
  if (output_true != NULL) {
    fp_true = fopen(output_true, "w");
    if (fp_true == NULL) {
      fprintf(stderr, "error: faield to create output true file\n");
      return -1;
    }
  }

  Rng random(seed);
  
  for (auto &p : data.data) {

    double tstar = p.predicted_tstar_synthetic(model);

    fprintf(fp_out, "%15.9f %6.3f %d\n", tstar + random.normal(noise_sigma), noise_sigma, (int)p.points.size());

    for (auto &pp : p.points) {

      double lon = pp.theta * 180.0/M_PI;
      double lat = 90.0 - pp.phi * 180.0/M_PI;

      fprintf(fp_out, "%15.9f %15.9f %15.9f\n", lon, lat, pp.r);
    }

    if (fp_true != NULL) {
      fprintf(fp_true, "%15.9f %6.3f %d\n", tstar, noise_sigma, (int)p.points.size());
      for (auto &pp : p.points) {
	
	double lon = pp.theta * 180.0/M_PI;
	double lat = 90.0 - pp.phi * 180.0/M_PI;
	
	fprintf(fp_true, "%15.9f %15.9f %15.9f\n", lon, lat, pp.r);
      }
    }
    
  }

  fclose(fp_out);

  if (image_output != nullptr) {

    FILE *fp_image = fopen(image_output, "w");
    if (fp_image == NULL) {
      fprintf(stderr, "error: failed to create image output file\n");
      return -1;
    }

    for (int j = 0; j < image_height; j ++) {
      
      // North Pole to South Pole
      double imagephi = ((double)j + 0.5)/(double)image_height * M_PI;
      
      for (int i = 0; i < image_width; i ++) {
	
	// -180 .. 180
	double imagetheta = ((double)i + 0.5)/(double)image_width * 2.0 * M_PI - M_PI;

	fprintf(fp_image, "%15.9f ", model(imagephi, imagetheta));

      }

      fprintf(fp_image, "\n");
    }

    fclose(fp_image);
  }
  
  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of\n"
	  "\n"
	  " -i | --input <filename>           Input data to base recompute with synthetic model\n"
	  " -o | --output <filename>          Output file to write synthetic noisy observations to\n"
	  " -O | --output-true <filename>     Output file to write synthetic true observations to\n"
	  "\n"
	  " -m | --model <name>               Synthetic model to use\n"
	  " -l | --list-models                List available synthetic models and exit\n"
	  "\n"
	  " -n | --noise <float>              Std dev of gaussian noise to add to observations\n"
	  " -S | --seed <int>                 Random seed\n"
	  "\n"
	  " -I | --image-output <filename>    Write synthetic model image\n"
	  " -W | --image-width <int>          Image width\n"
	  " -H | --image-height <int>         Image height\n"
	  "\n"
	  " -h | --help                       Usage\n"
	  "\n",
	  pname);
}
