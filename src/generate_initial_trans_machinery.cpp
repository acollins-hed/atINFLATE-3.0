//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include "Common_Variables.hpp"
#include "Genotype.hpp"

void gen_trans_machry(struct Common_Variables * common_variables){
  std::bernoulli_distribution dist_bern(common_variables->binom_p);
  if(common_variables->mask){
    if(common_variables->binom_p == 0){
      common_variables->tRNA_State_bits.setZero();
      common_variables->tRNA_Mask_bits.setZero();
      common_variables->aaRS_State_bits.setZero();
      common_variables->aaRS_Mask_bits.setZero();
    } else if(common_variables->binom_p == 1){
      common_variables->tRNA_State_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->tRNA_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_State_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);										       } else{
      
      common_variables->tRNA_State_bits.setZero();
      common_variables->tRNA_Mask_bits.setZero();
      common_variables->aaRS_State_bits.setZero();
      common_variables->aaRS_Mask_bits.setZero();
      for(int i=0;i<common_variables->N_tRNA;i++){
	for(int j=0;j<common_variables->N_trajectory;j++){
	  for(int m=0;m<common_variables->N_int_interface;m++){
	    common_variables->tRNA_State_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	    common_variables->tRNA_Mask_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	  }
	}
      }
      
      for(int i=0;i<common_variables->N_aaRS;i++){
	for(int j=0;j<common_variables->N_trajectory;j++){
	  for(int m=0;m<common_variables->N_int_interface;m++){
	    common_variables->aaRS_State_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	    common_variables->aaRS_Mask_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	  }
	}
      }
    }
  }
  else{
    if(common_variables->binom_p == 0){
      common_variables->tRNA_State_bits.setZero();
      common_variables->tRNA_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_State_bits.setZero();
      common_variables->aaRS_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);
    } else if(common_variables->binom_p == 1){
      common_variables->tRNA_State_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->tRNA_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_State_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);										       } else{
      
      common_variables->tRNA_State_bits.setZero();
      common_variables->tRNA_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);
      common_variables->aaRS_State_bits.setZero();
      common_variables->aaRS_Mask_bits.setConstant((1<<(common_variables->N_int_interface))-1);	
      for(int i=0;i<common_variables->N_tRNA;i++){
	for(int j=0;j<common_variables->N_trajectory;j++){
	  for(int m=0;m<common_variables->N_int_interface;m++){
	    common_variables->tRNA_State_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	  }
	}
      }
      
      for(int i=0;i<common_variables->N_aaRS;i++){
	for(int j=0;j<common_variables->N_trajectory;j++){
	  for(int m=0;m<common_variables->N_int_interface;m++){
	    common_variables->aaRS_State_bits(j,i) += (dist_bern(common_variables->mersene_twister)<<m);
	  }
	}
      }
    }
  }
}

//generate_initial_trans_machinery.cpp is part of atinflate_3
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
