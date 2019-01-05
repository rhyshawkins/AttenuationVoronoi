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
