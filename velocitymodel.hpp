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
