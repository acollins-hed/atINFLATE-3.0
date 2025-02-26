#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"
#include "tRNAs.hpp"
#include "aaRSs.hpp"
#include "fitness_functions.hpp"

struct Genotype
{
  tRNAs trnas;
  aaRSs aarss;
  Eigen::MatrixXd code;
  Eigen::MatrixXd kd;
  Eigen::MatrixXd matches;

  Genotype();
  
  Genotype(tRNAs t, aaRSs a,struct Common_Variables * common_variables);

  void Get_Code(struct Common_Variables * common_variables);
  
  Eigen::MatrixXd Get_kd(struct Common_Variables * common_variables);
};

Genotype::Genotype(){
  tRNAs t;
  aaRSs a;
  trnas=t;
  aarss=a;
  matches.resize(trnas.iis.cols(),aarss.iis.cols());
  struct Common_Variables ec;
  ec.kmax = 10000;
  ec.epsilon = 0.954178206405955;
  Get_Code(&ec);
}

Genotype::Genotype(tRNAs t, aaRSs a,struct Common_Variables * common_variables)
  :trnas(t),aarss(a){
  matches.resize(trnas.iis.cols(),aarss.iis.cols());
  Get_Code(common_variables);
}

void Genotype::Get_Code(struct Common_Variables * common_variables){
  kd=Get_kd(common_variables);
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

Eigen::MatrixXd Genotype::Get_kd(struct Common_Variables * common_variables){
  Eigen::MatrixXd krate(trnas.iis.cols(),aarss.iis.cols());
  krate.setZero();
  
  for(int i=0;i<trnas.iis.cols();i++){
    for(int j=0;j<aarss.iis.cols();j++){
      matches(i,j) = __builtin_popcount((trnas.iis(1,i)&aarss.iis(1,j))&(~(trnas.iis(0,i)^aarss.iis(0,j))));
      krate(i,j) = common_variables->kmax*exp(-1*common_variables->epsilon*(matches(i,j)));
    }
  }    

  return krate;
}
