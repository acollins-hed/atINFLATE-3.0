//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"
//#include "Stats.hpp"
#include "inverted_hsum.hpp"
#include "Transition_Probability.hpp"

double A_fit(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, int stype, struct Common_Variables * common_variables){
  double f=0;

  for(int aa=0;aa<eff_mutant_code->cols();aa++){
    for(int codon=0;codon<eff_mutant_code->rows();codon++){
      f += common_variables->selection_mat(stype,common_variables->aa_to_st[common_variables->amino_acids(aa)])*eff_mutant_code->coeff(codon,aa)*codon_freq->coeff(stype,codon);
    }
  }
  return f;
}

double Fitness_rate_indep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, struct Common_Variables * common_variables){
  double f=1;
  
  for(int stype=0;stype<common_variables->N_site_type;stype++){
    f *= pow(A_fit(codon_freq,eff_mutant_code,stype,common_variables),common_variables->site_type_freqs(stype));
  }
  
  return f;
}

double Fitness_rate_dep(Eigen::MatrixXd* codon_freq, Eigen::MatrixXd* eff_mutant_code, Eigen::MatrixXd* kd, struct Common_Variables * common_variables){
  
  double kdbar=0;
  
  //As described in Collins-Hed, Ardell 2019
  for(int i=0; i<eff_mutant_code->rows();i++){
    double sum = 0;
    for(int j=0; j<eff_mutant_code->cols(); j++){
      sum += kd->coeff(i,j);
    }
    kdbar+=1/sum;
  }
  kdbar = eff_mutant_code->rows()/kdbar;

  //////////////////
  
  return Fitness_rate_indep(codon_freq,eff_mutant_code,common_variables)*(1 - exp(-common_variables->rate_constant*kdbar));  
}

//fitness_functions.cpp is part of atinflate_3
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
