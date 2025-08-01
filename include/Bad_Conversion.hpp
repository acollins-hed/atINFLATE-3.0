#ifndef CONVERSIONS_HPP
#define CONVERSIONS_HPP

#include<iostream>
#include<sstream>


struct Bad_Conversion : std::runtime_error {
  Bad_Conversion(const std::string& s);
};
/*
inline double Convert_To_Double(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw Bad_Conversion("Convert_To_Double(\"" + s + "\")");
  return x;
}
*/
#endif
