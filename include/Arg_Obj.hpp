#ifndef ARG_OBJ_HPP
#define ARG_OBJ_HPP

#include<iostream>
#include<cctype>
#include<typeinfo>
#include<vector>
#include<climits>
#include<sstream>
#include "Convert_To_Double.hpp"
#include "Bad_Conversion.hpp"
#include "Common_Variables.hpp"

template <typename T> struct Arg_Obj{

  const std::string arg_name;
  std::string arg_name2;
  std::string arg_description;
  std::string arg_error_message;
  bool limit_error;
  T value;
  T upper_lim;
  T lower_lim;
  int count;
  Arg_Obj(std::string name, T v, T ll, T ul, std::string desc, std::string em, std::string& mp);
  Arg_Obj(std::string name, std::string name2, T v, T ll, T ul, std::string desc, std::string em, std::string& mp);
  void check_limits();
  void add_to_manpage(std::string& mp);
};

#endif
