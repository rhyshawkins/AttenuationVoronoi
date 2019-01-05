#pragma once
#ifndef attenuationutil_hpp
#define attenuationutil_hpp

#include <ostream>
#include <istream>

#include <stdio.h>

//
// Templated IO of values that can be over-ridden
//
template <typename value>
bool value_read(std::istream &s, value &v) {
  s.read((char*)&v, sizeof(v));
    
  return true;
}

template <typename value>
bool value_read(FILE *fp, value &v) {
  if (fread(&v, sizeof(v), 1, fp) != 1) {
    return false;
  }

  return true;
}

template <typename value>
bool value_write(std::ostream &s, const value &v)
{
  s.write((char*)&v, sizeof(v));
  
  return true;
}

template <typename value>
bool value_write(FILE *fp, const value &v) {
  if (fwrite(&v, sizeof(v), 1, fp) != 1) {
    return false;
  }

  return true;
}

template <typename value>
double scalartodouble(const value &v)
{
  return (double)v;
}

template <typename value>
value doubletoscalar(const double &d)
{
  return (value)d;
}

#endif // attenuationutil_hpp
