#include<Eigen/Dense>
#include<Eigen/Core>

double inverted_hsum(Eigen::MatrixXd k, int i, int A){
  double sum=0;
  for(int j=0;j<A;j++){
	sum += 1/k.coeff(i,j);
  }
  return (double(1)/(sum));
}
