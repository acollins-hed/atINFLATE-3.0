#ifndef CONVERSIONS_HPP
#define CONVERSIONS_HPP

#include<iostream>
#include<sstream>


struct BadConversion : std::runtime_error {
  BadConversion(const std::string& s);
};
/*
inline double convertToDouble(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}
*/
#endif
