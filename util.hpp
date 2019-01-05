#pragma once
#ifndef util_hpp
#define util_hpp

template
<typename T>
void encode(std::ostream &s, const T &v)
{
  s.write((char*)&v, sizeof(v));
}

template
<typename T>
void decode(std::istream &s, T &v)
{
  s.read((char*)&v, sizeof(v));
}

template
<typename value>
double todouble(const value &v) {
  return (double)v;
}

template
<typename value>
void fromdouble(double d, value &v) {
  v = (value)d;
}

#endif // util_hpp
