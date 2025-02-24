#ifndef POPULATION_HPP
#define POPULATION_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"
#include "Genotype.hpp"
#include "Mutant.hpp"

struct Population
{
  Genotype genotype;
  double fitness;
  Eigen::MatrixXd codon_frequency;
  Population(struct Common_Variables * common_variables);

  Population(struct Genotype g, struct Common_Variables * common_variables);

  void Get_Codon_Freq(struct Common_Variables * common_variables);

  void Sum_to_one(Eigen::MatrixXd * emat);

  Eigen::MatrixXd Diagonal(struct Common_Variables * common_variables);

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals);
};

#endif
