#ifndef POPULATION_HPP
#define POPULATION_HPP

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"
#include "Genotype.hpp"
#include "Mutant.hpp"

struct Population
{
  Genotype genotype; 
  double fitness;
  double mu_per_codon;
  double phi;
  int N_site_type;
  std::map<double, int> aa_to_st;
  Eigen::VectorXd site_types;
  Eigen::MatrixXd mutation_mat;
  Eigen::MatrixXd selection_mat;
  Eigen::MatrixXd codon_frequency;
  Population();

  Population(struct Genotype g, struct Common_Variables * common_variables);

  void Get_Codon_Freq();

  void Sum_to_one(Eigen::MatrixXd * emat);

  Eigen::MatrixXd Diagonal();

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals);
};

#endif
