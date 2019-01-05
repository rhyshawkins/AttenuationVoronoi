
#include <stdio.h>

#include "pathutil.hpp"

void mkpath(const char *prefix, const char *filename, char *path)
{
  if (prefix != nullptr) {
    sprintf(path, "%s%s", prefix, filename);
  } else {
    sprintf(path, "%s", filename);
  }
}

void mkrankpath(int rank, const char *prefix, const char *filename, char *path)
{
  if (prefix != nullptr) {
    sprintf(path, "%s%s-%03d", prefix, filename, rank);
  } else {
    sprintf(path, "%s-%03d", filename, rank);
  }
}
