#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Genotype.hpp"
#include "Mutant.hpp"

struct Population
{
  Genotype genotype;
  double fitness;
  Eigen::MatrixXd codon_frequency;
  Population(struct Common_Variables * common_variables);

  Population(struct Genotype g, struct Common_Variables * common_variables);

  void Get_Codon_Freq(struct Common_Variables * common_variables);

  void Sum_to_one(Eigen::MatrixXd * emat);

  Eigen::MatrixXd Diagonal(struct Common_Variables * common_variables);

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals);
};

Population::Population(struct Common_Variables * common_variables){
  Genotype genotype;
  Get_Codon_Freq(common_variables);
  //fitness = fitness_func(&codon_frequency, &genotype.code, &genotype.kd);
  fitness = 0;
}
Population::Population(struct Genotype g, struct Common_Variables * common_variables)
  :genotype(g){
  Get_Codon_Freq(common_variables);
  //fitness = fitness_func(&codon_frequency, &genotype.code, &genotype.kd);
  fitness = 0;
}
void Population::Get_Codon_Freq(struct Common_Variables * common_variables){
  Eigen::MatrixXd Q(genotype.trnas.iis.cols(),genotype.trnas.iis.cols());
  Eigen::MatrixXd wfit(genotype.trnas.iis.cols(),genotype.trnas.iis.cols());
  Eigen::MatrixXd diag(common_variables->N_site_type,genotype.trnas.iis.cols());
  diag = Diagonal(common_variables);
  codon_frequency.resize(common_variables->N_site_type,genotype.trnas.iis.cols());
  for(int i=0;i<common_variables->N_site_type;i++){
    wfit.setZero();
    wfit.diagonal() = diag.row(i).array();
    Q = common_variables->mutation_mat*wfit;
    Eigen::EigenSolver<Eigen::MatrixXd> ev(Q);
    //I have seen claims on the internet that Eigen automatically puts the
    //largest eigenvalue into cell 0 (I've not seen this on Eigen's own
    //website) and I've found this to be untrue, hence the following line.
    codon_frequency.row(i) = ev.eigenvectors().col(max_eigenvalue(&ev.eigenvalues())).real();
    //The following makes sure the components sum to one. The eigenvectors are
    //already normalized but of course that just makes their norm 1, not their
    //components sum to one.
    Sum_to_one(&codon_frequency);
  }
}
void Population::Sum_to_one(Eigen::MatrixXd * emat){
  double sum=0;
  for(int i=0;i<emat->rows();i++){
    for(int j=0;j<emat->cols();j++){
      sum = emat->coeff(i,j) + sum;
    }
    emat->row(i) *= (1/sum);
    sum = 0;
  }
}

Eigen::MatrixXd Population::Diagonal(struct Common_Variables * common_variables){
  Eigen::MatrixXd diag(common_variables->N_site_type,genotype.trnas.iis.cols());
  for(int stype = 0;stype<common_variables->N_site_type;stype++){
    for(int codon=0;codon<genotype.trnas.iis.cols();codon++){
      diag(stype,codon) = 0;
      for(int alpha=0;alpha<common_variables->N_aaRS;alpha++){
	diag(stype,codon) += common_variables->selection_mat(stype,common_variables->aa_to_st[common_variables->amino_acids(alpha)])*genotype.code(codon,alpha);
      }
    }
  }
  return diag;
}

int Population::max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals){
  int max = 0;
  for(int i=1;i<evals->size();i++){
    if(evals->coeff(max).real() < evals->coeff(i).real()){
      max = i;
    }
  }
  return max;
}
