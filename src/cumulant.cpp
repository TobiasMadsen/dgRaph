#include <Rcpp.h>
using namespace Rcpp;

double c0(double t, double p, double s){
  return log(1-p+p*exp(s*t));
}

double c1(double t, double p, double s){
  double num = s*p*exp(s*t);
  double den = 1-p+p*exp(s*t);
  return num/den;
}

double c2(double t, double p, double s){
  double num = s*s*p*exp(s*t)*(1-p);
  double sqrtden = (1-p+p*exp(s*t));
  return num/sqrtden/sqrtden;
}

// [[Rcpp::export]]
double cumulant(double t, NumericVector p, NumericVector s) {
  /* Return cumulant generating function*/
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c0(t, p[i], s[i]);
  return ret;
}

// [[Rcpp::export]]
double cumulantD1(double t, NumericVector p, NumericVector s) {
  /* Return first derivative of cumulant generating function */
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c1(t, p[i], s[i]);
  return ret;
}

// [[Rcpp::export]]
double cumulantD2(double t, NumericVector p, NumericVector s) {
  /* Return second derivative of cumulant generating function */
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c2(t, p[i], s[i]);
  return ret;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# Check cumulant
cumulant(1., 0.4, 1.)
log( 0.6+0.4*exp(1))

# Check cumulantD1
cumulantD1(1., 0.4, 1.)
0.4*exp(1)/(0.4*exp(1)+0.6*1)

# Check first derivative
(cumulant(1.1, 0.3, 1.) - cumulant(1., 0.3, 1.))/0.1

cumulantD1(1.05, 0.3, 1.)

(cumulant(1.6, 0.3, 2.) - cumulant(1.5, 0.3, 2.))/0.1

cumulantD1(1.55, 0.3, 2.)

# Check second derivative
(cumulantD1(1.6, 0.3, 2.) - cumulantD1(1.5, 0.3, 2.))/0.1

cumulantD2(1.55, 0.3, 2.)
*/
