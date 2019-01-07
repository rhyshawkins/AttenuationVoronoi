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
#ifndef velocitymodel_hpp
#define velocitymodel_hpp

template
<typename value>
value pwave_velocity(value r)
{
  //
  // Example polynomial (put in your one and uncomment)
  //

  return 1.0e-6 *r*r*r + 1.0e-5*r*r + 1.0e-3*r*r + 11.0;


  // Constant model
  //return 11.0;
  
}


#endif // velocitymodel_hpp
