#include<iostream>
#include<map>
#include<Eigen/Dense>
#include<Eigen/Core>
#include "Genotype.hpp"
#include "Mutant.hpp"

struct Population
{
  Genotype genotype;
  double fitness;
  double mu_per_codon;
  double phi;
  int N_site_type;
  std::map<double, int> aa_to_st;
  Eigen::VectorXd site_types;
  Eigen::MatrixXd mutation_mat;
  Eigen::MatrixXd selection_mat;
  Eigen::MatrixXd codon_frequency;
  Population();

  Population(struct Genotype g, struct Common_Variables * common_variables);

  void Get_Codon_Freq();

  void Sum_to_one(Eigen::MatrixXd * emat);

  Eigen::MatrixXd Diagonal();

  int max_eigenvalue(const Eigen::EigenSolver<Eigen::MatrixXcd>::EigenvalueType * evals);
};

Population::Population(){
  Genotype genotype;
  N_site_type = genotype.aarss.iis.cols();
  mu_per_codon = 0.0001;
  phi = 0.99;
  site_types.resize(N_site_type);
  mutation_mat.resize(genotype.trnas.iis.cols(),genotype.trnas.iis.cols());
  selection_mat.resize(N_site_type,N_site_type);
  for(int i = 0;i<genotype.trnas.iis.cols();i++){
    for(int j = 0; j<genotype.trnas.iis.cols();j++){
      if(i == j)
	mutation_mat(i,j) = 1 - 2*mu_per_codon;
      else{
	if(j == i+1 || j == i-1 || (j == 0 && i == genotype.trnas.iis.cols() - 1) || (j == genotype.trnas.iis.cols() - 1 && i == 0))
	  mutation_mat(i,j) = mu_per_codon;
	else
	  mutation_mat(i,j) = 0;
      }
    }
  }

  for(int i=0;i<N_site_type;i++)
    site_types(i) = ((double) i)/((double) N_site_type-1);
  
  for(int i=0;i<N_site_type;i++){
    for(int j=0;j<N_site_type;j++){
      selection_mat(i,j) = pow(phi,abs(site_types(i)-site_types(j)));
    }
  }
  
  Get_Codon_Freq();
  fitness = 0;
}

Population::Population(struct Genotype g, struct Common_Variables * common_variables)
  :genotype(g){
  N_site_type = common_variables->N_site_type;
  mu_per_codon = common_variables->mu_per_codon;
  phi = common_variables->phi;
  aa_to_st = common_variables->aa_to_st;
  site_types.resize(N_site_type);
  site_types = common_variables->site_types;
  selection_mat.resize(N_site_type,N_site_type);
  selection_mat = common_variables->selection_mat;
  mutation_mat.resize(genotype.trnas.iis.rows(),genotype.trnas.iis.cols());
  mutation_mat = common_variables->mutation_mat;
  Get_Codon_Freq();
  fitness = 0;
}

void Population::Get_Codon_Freq(){
  Eigen::MatrixXd Q(genotype.trnas.iis.cols(),genotype.trnas.iis.cols());
  Eigen::MatrixXd wfit(genotype.trnas.iis.cols(),genotype.trnas.iis.cols());
  Eigen::MatrixXd diag(N_site_type,genotype.trnas.iis.cols());
  diag = Diagonal();
  codon_frequency.resize(N_site_type,genotype.trnas.iis.cols());
  for(int i=0;i<N_site_type;i++){
    wfit.setZero();
    wfit.diagonal() = diag.row(i).array();
    Q = mutation_mat*wfit;
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

Eigen::MatrixXd Population::Diagonal(){
  Eigen::MatrixXd diag(N_site_type,genotype.trnas.iis.cols());
  for(int stype = 0;stype<N_site_type;stype++){
    for(int codon=0;codon<genotype.trnas.iis.cols();codon++){
      diag(stype,codon) = 0;
      for(int alpha=0;alpha<genotype.aarss.iis.cols();alpha++){
	diag(stype,codon) += selection_mat(stype,aa_to_st[genotype.aarss.aas(alpha)])*genotype.code(codon,alpha);
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
