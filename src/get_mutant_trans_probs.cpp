//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<algorithm>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"
#include "Mutant.hpp"
#include "fitness_functions.hpp"
#include "Transition_Probability.hpp"
#include "Population.hpp"


void get_trans_prob_rd(std::vector<Mutant*> mut_vec, struct Population * population, struct Common_Variables * common_variables, long double * sum){

  long double sum1=0;
  int hamming_dist = 1;
  
  for(unsigned int i=0;i<mut_vec.size();i++){
    struct Mutant mutant = *(mut_vec[i]);
    struct Population mutant_pop = *population;
    if(mutant.mutation_1.kind_trans_machinery == 't'){
      mutant_pop.genotype.trnas.iis(mutant.mutation_1.mask,mutant.mutation_1.which_trans_machinery) ^= 1<<mutant.mutation_1.mutation_position;
    } else{
      mutant_pop.genotype.aarss.iis(mutant.mutation_1.mask,mutant.mutation_1.which_trans_machinery) ^= 1<<mutant.mutation_1.mutation_position;
    }
    if(mutant.mutation_2.kind_trans_machinery != '0'){
      hamming_dist = 2;
      if(mutant.mutation_2.kind_trans_machinery == 't'){
	mutant_pop.genotype.trnas.iis(mutant.mutation_2.mask,mutant.mutation_2.which_trans_machinery) ^= 1<<mutant.mutation_2.mutation_position;
      } else{
	mutant_pop.genotype.aarss.iis(mutant.mutation_2.mask,mutant.mutation_2.which_trans_machinery) ^= 1<<mutant.mutation_2.mutation_position;
      }
    }
    mutant_pop.genotype.Get_Code();
    mutant_pop.fitness = Fitness_rate_dep(&(population->codon_frequency),&mutant_pop.genotype.code,&mutant_pop.genotype.kd,common_variables);
    
    long double trans_prob = Transition_Probability(population->fitness,mutant_pop.fitness,hamming_dist,common_variables);
    mut_vec[i]->trans_prob = trans_prob;
    sum1 += trans_prob;

  }
  *sum = sum1;
}

void get_trans_prob_ri(std::vector<Mutant*> mut_vec, struct Population * population, struct Common_Variables * common_variables, long double * sum){

  long double sum1=0;
  int hamming_dist = 1;
  
  for(unsigned int i=0;i<mut_vec.size();i++){
    struct Mutant mutant = *(mut_vec[i]);
    struct Population mutant_pop = *population;
    if(mutant.mutation_1.kind_trans_machinery == 't'){
      mutant_pop.genotype.trnas.iis(mutant.mutation_1.mask,mutant.mutation_1.which_trans_machinery) ^= 1<<mutant.mutation_1.mutation_position;
    } else{
      mutant_pop.genotype.aarss.iis(mutant.mutation_1.mask,mutant.mutation_1.which_trans_machinery) ^= 1<<mutant.mutation_1.mutation_position;
    }
    if(mutant.mutation_2.kind_trans_machinery != '0'){
      hamming_dist = 2;
      if(mutant.mutation_2.kind_trans_machinery == 't'){
	mutant_pop.genotype.trnas.iis(mutant.mutation_2.mask,mutant.mutation_2.which_trans_machinery) ^= 1<<mutant.mutation_2.mutation_position;
      } else{
	mutant_pop.genotype.aarss.iis(mutant.mutation_2.mask,mutant.mutation_2.which_trans_machinery) ^= 1<<mutant.mutation_2.mutation_position;
      }
    }
    mutant_pop.genotype.Get_Code();
    mutant_pop.fitness = Fitness_rate_indep(&(population->codon_frequency),&mutant_pop.genotype.code,common_variables);
    
    long double trans_prob = Transition_Probability(population->fitness,mutant_pop.fitness,hamming_dist,common_variables);
    mut_vec[i]->trans_prob = trans_prob;
    sum1 += trans_prob;

  }
  *sum = sum1;
}

//get_mutant_trans_probs.cpp is part of atinflate_3
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
