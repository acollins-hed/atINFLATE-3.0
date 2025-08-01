#include<sstream>
#include "Bad_Conversion.hpp"

double Convert_To_Double(const std::string& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    throw Bad_Conversion("Convert_To_Double(\"" + s + "\")");
  return x;
}
