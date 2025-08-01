#ifndef AARS_HPP
#define AARS_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<map>
#include "Common_Variables.hpp"

struct aaRSs
{
  Eigen::VectorXd aas;
  Eigen::MatrixXi iis;
  int N_int_interface;
  aaRSs();
  aaRSs(int traj, struct Common_Variables * common_variables);
  void  print();
  std::string get_aars_int(int i, int j);
};

#endif
