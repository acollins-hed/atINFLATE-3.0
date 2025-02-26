#include<iostream>
#include<sstream>

struct BadConversion : std::runtime_error {
  BadConversion(const std::string& s);
};

BadConversion::BadConversion(const std::string& s)
    : std::runtime_error(s)
    { }

