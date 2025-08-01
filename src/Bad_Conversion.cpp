#include<iostream>
#include<sstream>

struct Bad_Conversion : std::runtime_error {
  Bad_Conversion(const std::string& s);
};

Bad_Conversion::Bad_Conversion(const std::string& s)
    : std::runtime_error(s)
    { }

