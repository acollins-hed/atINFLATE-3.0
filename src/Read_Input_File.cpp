//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Bad_Conversion.hpp"
#include "Convert_To_Double.hpp"
#include "Common_Variables.hpp"

int Read_Checkpoint(struct Common_Variables * common_variables, bool bl_halting_fixation){
  std::ifstream icheckpoint_file;
  //icheckpoint_file.open(common_variables->checkpoint_filename);
  icheckpoint_file.open(common_variables->checkpoint_filename+"_checkpoint.log");
  if(icheckpoint_file){
    std::string cp;
    icheckpoint_file>>cp;
    common_variables->output_filename = cp;
    icheckpoint_file>>cp;
    common_variables->N_int_interface = stoi(cp);
    //icheckpoint_file>>cp;
    //k = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->N_tRNA = stoi(cp);
    //icheckpoint_file>>cp;
    //Tcap = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->N_aaRS = stoi(cp);
    //icheckpoint_file>>cp;
    //Acap = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->N_site_type = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->phi = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->transition_bias = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->N_population = stoi(cp);
    common_variables->site_type_freqs.resize(common_variables->N_site_type);
    for(int i=0;i<common_variables->N_site_type;i++){
      icheckpoint_file>>cp;
      common_variables->site_type_freqs(i) = stoi(cp);
    }
    icheckpoint_file>>cp;
    common_variables->rate = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->mask = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->proofreading = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->rate_constant = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->mu_id_feat = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->mu_per_codon = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->kmax = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->kmin = Convert_To_Double(cp);
    if(common_variables->proofreading){
      common_variables->kmax *= common_variables->kmax;
      common_variables->kmin *= common_variables->kmin;
    }
    common_variables->site_types.resize(common_variables->N_site_type);
    for(int i = 0;i<common_variables->N_site_type;i++){
      icheckpoint_file>>cp;
      common_variables->site_types(i) = Convert_To_Double(cp);
    }
    icheckpoint_file>>cp;
    common_variables->codon_ring_space = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->N_base = stoi(cp);
    if(!(common_variables->codon_ring_space) && common_variables->N_base == 2)
      common_variables->codonspace2=true;
    else{
      if(common_variables->N_base == 4)
	common_variables->codonspace4=true;
    }
    icheckpoint_file>>cp;
    common_variables->binom_p = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->N_trajectory = stoi(cp);
    common_variables->tRNA_State_bits.resize(common_variables->N_trajectory,common_variables->N_tRNA);
    common_variables->aaRS_State_bits.resize(common_variables->N_trajectory,common_variables->N_aaRS);
    common_variables->tRNA_Mask_bits.resize(common_variables->N_trajectory,common_variables->N_tRNA);
    common_variables->aaRS_Mask_bits.resize(common_variables->N_trajectory,common_variables->N_aaRS);
    icheckpoint_file>>cp;
    common_variables->end_fixation = stoi(cp);
    icheckpoint_file>>cp;
    common_variables->halting_fitness = Convert_To_Double(cp);
    icheckpoint_file>>cp;
    common_variables->N_threads = stoi(cp);
    if(!bl_halting_fixation){
      std::cout<<"You need to use --halting-fixation with -i or --ifile. See --help or -h for more.\n\n";
      return 0;
    }
    else{
      if(common_variables->halting_fixation < common_variables->end_fixation){
	std::cout<<"The new halting fixation, --halting-fixation, must be larger than the ending fixation in the checkpoint.log file. See --help or -h for more.\n\n";
	return 0;
      }
    }
    common_variables->amino_acids.resize(common_variables->N_aaRS);
    for(int i=0;i<common_variables->N_aaRS;i++){
      common_variables->amino_acids(i) = common_variables->site_types(i);
      common_variables->aa_to_st[common_variables->amino_acids(i)] = i;
    }
    
    int index=0;
    for(int i = 0;i<4*common_variables->N_trajectory;i++){
      if(i != 0 && i%4 == 0)
	index++;
      if(i%4 == 0){
	for(int j = 0;j<common_variables->N_tRNA;j++){
	  icheckpoint_file>>cp;
	  common_variables->tRNA_State_bits(index,j) = stoi(cp);
	}
      }
      else{
	if(i%4 == 1){
	  for(int j = 0;j<common_variables->N_tRNA;j++){
	    icheckpoint_file>>cp;
	    common_variables->tRNA_Mask_bits(index,j) = stoi(cp);
	  }
	}
	else{
	  if(i%4 == 2){
	    for(int j = 0;j<common_variables->N_aaRS;j++){
	      icheckpoint_file>>cp;
	      common_variables->aaRS_State_bits(index,j) = stoi(cp);
	    }
	  }
	  else{
	    for(int j = 0;j<common_variables->N_aaRS;j++){
	      icheckpoint_file>>cp;
	      common_variables->aaRS_Mask_bits(index,j) = stoi(cp);
	    }
	  }
	}
      }
    }
  } else{
    std::cout<<"The file "<<common_variables->checkpoint_filename+"_checkpoint.log"<<" does not exist\n";
    return 0;
  }
  
  icheckpoint_file.close();
  return 1;
}

//Read_Input_File.cpp is part of atinflate_3
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
