//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"

struct tRNAs
{
  Eigen::VectorXi anticodons;
  //The iis of the tRNAs are 2xT matrices.
  Eigen::MatrixXi iis;
  int N_int_interface;
  
  tRNAs();
  
  tRNAs(int traj, struct Common_Variables * common_variables);

  void  print();

  std::string get_trna_int(int i, int j);

};

tRNAs::tRNAs(){
  Eigen::MatrixXi t(2,4);
  Eigen::VectorXi a(4);
  t.setZero();
  iis = t;
  N_int_interface = 4;
  for(int i=0;i<4;i++)
    a(i) = i;
  anticodons=a;
}

tRNAs::tRNAs(int traj, struct Common_Variables * common_variables){
  iis.resize(2,common_variables->N_tRNA);
  iis.row(0) = common_variables->tRNA_State_bits.row(traj);
  iis.row(1) = common_variables->tRNA_Mask_bits.row(traj);
  N_int_interface = common_variables->N_int_interface;
  Eigen::VectorXi z(iis.cols());
  for(int i=0;i<iis.cols();i++)
    z(i)=i;
  anticodons = z;
  
}

void  tRNAs::print(){
  for(int i = 0;i<2;i++){
    for(int j=0;j<iis.cols();j++){
      for(int ii=N_int_interface-1;ii>=0;ii--){
	if(iis(i,j)&(1<<ii))
	  std::cout<<"1";
	else
	  std::cout<<"0";
      }
      std::cout<<" ";
      }
    std::cout<<std::endl;
  }
}

std::string tRNAs::get_trna_int(int i, int j){
  std::string s = "";
  for(int ii=N_int_interface-1;ii>=0;ii--){
    if(iis(i,j)&(1<<ii))
      s += "1";
    else
      s += "0";
  }
  return s;
}

//tRNAs.cpp is part of atinflate_3
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
