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

#include <stdio.h>
#include <stdlib.h>

#include <getopt.h>

#include "coordinate.hpp"
#include "sphericalvoronoimodel.hpp"
#include "rng.hpp"
#include "sphericalprior.hpp"

typedef sphericalcoordinate<double> coord_t;

static char short_options[] = "i:o:W:H:N:v:V:h";

static struct option long_options[] = {
  {"output", required_argument, 0, 'o'},
  {"lonsamples", required_argument, 0, 'W'},
  {"latsamples", required_argument, 0, 'H'},
  {"points", required_argument, 0, 'N'},
  {"vmin", required_argument, 0, 'v'},
  {"vmax", required_argument, 0, 'V'},
  {"help", no_argument, 0, 'h'},
  {0, 0, 0, 0}
 
};

static void usage(const char *pname);

static int saveimage(const char *filename, double *image, int width, int height);

int main(int argc, char *argv[])
{
  int c;
  int option_index;
  
  int lonsamples;
  int latsamples;
  
  //
  // Output Files
  //
  char *output;

  int npoints;
  double vmin;
  double vmax;
  
  //
  // State
  //
  int imagesize;
  double *image;

  //
  // Defaults
  //

  lonsamples = 16;
  latsamples = 16;

  output = nullptr;

  npoints = 10;
  vmin = 0.0;
  vmax = 1.0;
  
  option_index = 0;
  while (1) {

    c = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {

    case 'o':
      output = optarg;
      break;

    case 'W':
      lonsamples = atoi(optarg);
      if (lonsamples < 1) {
	fprintf(stderr, "error: xsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 'H':
      latsamples = atoi(optarg);
      if (latsamples < 1) {
	fprintf(stderr, "error: latsamples must be 1 or greater\n");
	return -1;
      }
      break;

    case 'N':
      npoints = atoi(optarg);
      if (npoints <= 0) {
	fprintf(stderr, "error: npoints must be 1 or more\n");
	return -1;
      }
      break;

    case 'v':
      vmin = atof(optarg);
      break;

    case 'V':
      vmax = atof(optarg);
      break;
      
    default:
      fprintf(stderr, "error: invalid option '%c'\n", c);
      
    case 'h':
      usage(argv[0]);
      return -1;
    }
  }

  if (output == nullptr) {
    fprintf(stderr, "error: required output file parameter missing\n");
    return -1;
  }

  imagesize = lonsamples * latsamples;
  image = new double[imagesize];

  sphericalvoronoimodel<double> model(false);

  Rng rng(983);
  UniformSphericalPrior prior;

  for (int i = 0; i < npoints; i ++) {

    double v = ((double)i + 0.5)/(double)npoints * (vmax - vmin) + vmin;

    double phi, theta;
    prior.sample(rng, phi, theta);

    model.add_cell(coord_t(phi, theta), v);
  }
  
  
  for (int j = 0; j < latsamples; j ++) {
    
    // North Pole to South Pole
    double imagephi = ((double)j + 0.5)/(double)latsamples * M_PI;
    for (int i = 0; i < lonsamples; i ++) {
      
      // -180 .. 180
      double imagetheta = ((double)i + 0.5)/(double)lonsamples * 2.0 * M_PI - M_PI;
      
      image[j * lonsamples + i] = model.value_at_point(coord_t(imagephi, imagetheta));
      
    }
  }
  
  //
  // Always save the mean
  //
  if (saveimage(output, image, lonsamples, latsamples) < 0) {
    fprintf(stderr, "error: failed to save mean\n");
    return -1;
  }

  char filename[1024];
  sprintf(filename, "%s.points", output);
  if (!model.save(filename)) {
    fprintf(stderr, "error: failed to save model\n");
    return -1;
  }
  

  return 0;
}

static void usage(const char *pname)
{
  fprintf(stderr,
	  "usage: %s [options]\n"
	  "where options is one or more of:\n"
	  "\n"
	  " -o|--output <filename>      Output likelihood file (required)\n"
	  "\n"
	  " -W|--lonsamples <int>       No. samples in longitude direction\n"
	  " -H|--latsamples <int>       No. samples in latitude direction\n"
	  "\n"
	  "\n"
	  " -h|--help                   Usage\n"
	  "\n",
	  pname);
}


static int saveimage(const char *filename, double *image, int width, int height)
{
  FILE *fp;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    return -1;
  }

  for (int j = 0; j < height; j ++) {
    for (int i = 0; i < width; i ++) {

      fprintf(fp, "%10.6f ", image[width * j + i]);

    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}
