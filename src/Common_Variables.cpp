//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<map>
#include "Mutation.hpp"
#include "Mutant.hpp"

struct Common_Variables{
  int N_aaRS;
  int N_tRNA;
  int N_int_interface;
  int N_trajectory;
  int halting_fixation;
  int end_fixation;
  int N_site_type;
  int N_population;
  int N_base;
  int N_total_mutants;
  int N_single_mutants;
  int N_threads;
  unsigned long int seed;
  double kmin;
  double kmax;
  double epsilon;
  double mu_id_feat;
  double mu_per_codon;
  double phi;
  double binom_p;
  double halting_fitness;
  double transition_bias;
  double rate_exponent;
  double rate_constant;
  //long double prob_of_1_mutation;
  //long double prob_of_2_mutations;
  std::string output_filename;
  std::string mutation_type;
  std::string checkpoint_filename;
  std::string params_filename;
  Eigen::VectorXi site_type_freqs;
  Eigen::VectorXd amino_acids;
  Eigen::VectorXd site_types;
  Eigen::MatrixXd selection_mat;
  Eigen::MatrixXd mutation_mat;
  Eigen::MatrixXi tRNA_State_bits;
  Eigen::MatrixXi tRNA_Mask_bits;
  Eigen::MatrixXi aaRS_State_bits;
  Eigen::MatrixXi aaRS_Mask_bits;
  std::random_device rand_dev;
  std::map<double, int> aa_to_st;
  std::mt19937_64 mersene_twister;
  std::uniform_real_distribution<long double> gillespie_rand{0,1};
  //std::uniform_real_distribution<double> unit_unif{0,1};
  bool rate;
  bool mask;
  bool codon_ring_space;
  bool codonspace2;
  bool codonspace4;
  bool proofreading;
  bool bl_rate_const;
  bool bl_checkpoint_filename;
  bool bl_params_filename;
  bool bl_double_mutants;
  struct Mutant fixed_mutant;
};

//Common_Variables.cpp is part of atinflate_3
//atinflate_3 simulates the aaRS-tRNA Interaction Fitness LAndscape
//Topographer Express (atINFLATE) model, a dynamic extension of the
//atINFLAT model published in Collins-Hed and Ardell 2019
//https://www.sciencedirect.com/science/article/pii/S0040580918300789
//Copyright (C) 2026 Andrea Collins-Hed
//
//This file is part of atinflate_3
//
//atinflate_3 is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//atinflate_3 is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program. If not, see <https://www.gnu.org/licenses/>. 
