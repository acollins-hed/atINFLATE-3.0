#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include "common_variables.hpp"
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
