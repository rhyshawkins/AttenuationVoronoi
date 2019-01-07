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
#ifndef attenuationexception_hpp
#define attenuationexception_hpp

#include <exception>

#define ATTENUATIONEXCEPTION(fmt, ...) attenuationexception(__FILE__, __FUNCTION__, __LINE__, fmt, ##__VA_ARGS__)

class attenuationexception : public std::exception {
public:

  
  attenuationexception(const char *srcfile,
		       const char *function,
		       int lineno,
		       const char *fmt, ...);
  ~attenuationexception();
  
};

#endif // attenuationexception_hpp
