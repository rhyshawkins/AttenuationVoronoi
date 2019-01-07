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
