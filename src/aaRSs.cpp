#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include<random>
#include<map>
#include "common_variables.hpp"

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
