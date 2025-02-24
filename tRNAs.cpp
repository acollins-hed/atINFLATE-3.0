#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "common_variables.hpp"

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
