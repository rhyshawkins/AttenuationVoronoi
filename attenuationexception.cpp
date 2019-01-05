
#include <stdio.h>
#include <stdarg.h>

#include "attenuationexception.hpp"

attenuationexception::attenuationexception(const char *srcfile,
				     const char *function,
				     int lineno,
				     const char *fmt, ...)
{
  va_list ap;
  
  fprintf(stderr, "Attenuation Exception: %s: %s: %d:", srcfile, function, lineno);

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  fprintf(stderr, "\n");
}

attenuationexception::~attenuationexception()
{
}
