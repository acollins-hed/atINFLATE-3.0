//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include "Common_Variables.hpp"

long double Transition_Probability(double current_fitness, double mutant_fitness,int hamming, struct Common_Variables * common_variables){
  long double Fixation_Probability;
  if(current_fitness == mutant_fitness)
    Fixation_Probability = 1/((long double) common_variables->N_population);
  else
    Fixation_Probability = (long double) (1 - (current_fitness)/(mutant_fitness))/(1 - pow((current_fitness)/(mutant_fitness),common_variables->N_population));
  return ((long double) common_variables->N_population*pow(common_variables->mu_id_feat,hamming)*Fixation_Probability); 
}

//Transition_Porbability.cpp is part of atinflate_3
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
