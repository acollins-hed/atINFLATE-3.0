//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include "Common_Variables.hpp"

int Get_Mutant_Vector_Size(struct Common_Variables * common_variables){
  int single_mutants;
  int double_mutants;

//Explanation for the calculation of total mutants:
  /*
  the number of possible single mutants is 2*(N_tRNA+N_aaRS)*N_int_interface 
  is the total number of mutants 1 Hamming distance away. The factor of 2
  comes from the fact that there are state bits and mask bits.
  
  the number of possible double mutants is single mutants choose 2
  */

  if(common_variables->mask)
    single_mutants = 2*(common_variables->N_tRNA+common_variables->N_aaRS)*common_variables->N_int_interface;
  else{
    single_mutants = (common_variables->N_tRNA+common_variables->N_aaRS)*common_variables->N_int_interface;
  }
  double_mutants = (single_mutants*(single_mutants-1))/2;
  
  if(common_variables->bl_double_mutants)
    return (double_mutants + single_mutants) ;
  else
    return single_mutants;
}

//Get_Mutant_Vector_Size.cpp is part of atinflate_3
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
