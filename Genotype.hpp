#ifndef GENOTYPE_HPP
#define GENOTYPE_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"
#include "tRNAs.hpp"
#include "aaRSs.hpp"
#include "fitness_functions.hpp"

struct Genotype
{
  tRNAs trnas;
  aaRSs aarss;
  Eigen::MatrixXd code;
  Eigen::MatrixXd kd;
  Eigen::MatrixXd matches;
  void Get_Code(struct Common_Variables * common_variables);

  Genotype();

  Genotype(tRNAs t, aaRSs a,struct Common_Variables * common_variables);

  Eigen::MatrixXd Get_kd(struct Common_Variables * common_variables);

};

#endif
