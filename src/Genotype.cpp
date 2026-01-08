//Copyright (C) 2026 Andrea Collins-Hed
//See end of file for extended copyright information.

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Common_Variables.hpp"
#include "tRNAs.hpp"
#include "aaRSs.hpp"
#include "fitness_functions.hpp"

struct Genotype
{
  double kmax;
  double epsilon;
  tRNAs trnas;
  aaRSs aarss;
  Eigen::MatrixXd code;
  Eigen::MatrixXd kd;
  Eigen::MatrixXd matches;

  Genotype();
  
  Genotype(tRNAs t, aaRSs a,struct Common_Variables * common_variables);

  void Get_Code();
  
  Eigen::MatrixXd Get_kd();
};

Genotype::Genotype(){
  tRNAs t;
  aaRSs a;
  trnas=t;
  aarss=a;
  matches.resize(trnas.iis.cols(),aarss.iis.cols());
  kmax = 10000;
  epsilon = 0.954178206405955;
  Get_Code();
}

Genotype::Genotype(tRNAs t, aaRSs a,struct Common_Variables * common_variables)
  :trnas(t),aarss(a){
  kmax = common_variables->kmax;
  epsilon = common_variables->epsilon;
  matches.resize(trnas.iis.cols(),aarss.iis.cols());
  Get_Code();
}

void Genotype::Get_Code(){
  kd=Get_kd();
  Eigen::MatrixXd c(trnas.iis.cols(),aarss.iis.cols());
  double hx;
  for(int i=0;i<trnas.iis.cols();i++){
    hx = inverted_hsum(kd,i,aarss.iis.cols());
    for(int j=0;j<aarss.iis.cols();j++){
      c(i,j) = hx/kd(i,j);
    }
  }
  code = c;
}

Eigen::MatrixXd Genotype::Get_kd(){
  Eigen::MatrixXd krate(trnas.iis.cols(),aarss.iis.cols());
  krate.setZero();
  
  for(int i=0;i<trnas.iis.cols();i++){
    for(int j=0;j<aarss.iis.cols();j++){
      matches(i,j) = __builtin_popcount((trnas.iis(1,i)&aarss.iis(1,j))&(~(trnas.iis(0,i)^aarss.iis(0,j))));
      krate(i,j) = kmax*exp(-1*epsilon*(matches(i,j)));
    }
  }    

  return krate;
}

//Genotype.cpp is part of atinflate_3
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
