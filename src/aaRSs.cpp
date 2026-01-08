//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<map>
#include "Common_Variables.hpp"

struct aaRSs
{
  Eigen::VectorXd aas;
  //The iis of the aaRSs are 2xA matrices.
  Eigen::MatrixXi iis;
  int N_int_interface;
  aaRSs();
  aaRSs(int traj, struct Common_Variables * common_variables);
  
  void  print();

  std::string get_aars_int(int i, int j);

};

aaRSs::aaRSs(){
  Eigen::MatrixXi a(2,4);
  Eigen::VectorXd aa(4);
  a.setZero();
  N_int_interface = 4;
  iis = a;
  for(int i=0;i<4;i++)
    aa(i) = ((float)i)/(4-1);
  aas=aa;
}
aaRSs::aaRSs(int traj, struct Common_Variables * common_variables){
  iis.resize(2,common_variables->N_aaRS);
  iis.row(0) = common_variables->aaRS_State_bits.row(traj);
  iis.row(1) = common_variables->aaRS_Mask_bits.row(traj);
  aas = common_variables->amino_acids;
  N_int_interface = common_variables->N_int_interface;
}
void  aaRSs::print(){
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

std::string aaRSs::get_aars_int(int i, int j){
  //This function returns either the state bits (row 0)
  //or the mask bits (row 1) of a specfic aaRS.
  std::string s = "";
  for(int ii=N_int_interface-1;ii>=0;ii--){
    if(iis(i,j)&(1<<ii))
      s += "1";
    else
      s += "0";
  }
  return s;
}

//aaRSs.cpp is part of atinflate_3
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
