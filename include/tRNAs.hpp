#ifndef TRNA_HPP
#define TRNA_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"

struct tRNAs
{
  Eigen::VectorXi anticodons;
  Eigen::MatrixXi iis;
  int N_int_interface;
  
  tRNAs();
  
  tRNAs(int traj, struct Common_Variables * common_variables);

  void  print();

  std::string get_trna_int(int i, int j);
};

#endif
