#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector EMbayes(NumericVector x, double alpha = 0, double beta = 0, double lambda = 0, double eps = 1e-8, bool verbose = true) {
  int n = x.size();
  double pi = 0.5;
  double s1 = sd(x);
  double s2 = s1;
  double mu1 = mean(x) - s1/4.0;
  double mu2 = mu1 + s2/2.0;
  double rho = 1;
  int k = 0;
  if(verbose) {
    Rcout << "Starting at:\n";
    Rcout << "pi = " << pi << ", mu1 = " << mu1 << ", mu2 = " << mu2 << ", s1² = " << s1*s1 << ", s2² = " << s2*s2 << "\n";
  }
  std::vector<double> tau(n);
  while(true) {
    k++;
    if(k % 100 == 0) checkUserInterrupt();
    double piOld = pi;
    for(int i = 0; i < n; i++)
      tau[i] = 1/(1 + (1-pi)/pi*s1/s2*exp(-0.5*( (x[i]-mu2)*(x[i]-mu2)/(s2*s2) - (x[i]-mu1)*(x[i]-mu1)/(s1*s1) )));
    
    double Stau = std::accumulate(tau.begin(), tau.end(), 0.0);
    pi = (Stau + alpha) / ((double) n + alpha + beta);
    
    mu1 = mu2 = 0;
    for(int i = 0; i < n; i++) {
      mu1 += tau[i]*x[i];
      mu2 += (1 - tau[i])*x[i];
    }
    mu1 /= Stau;
    mu2 /= ((double) n - Stau);
    
    double S1 = -Stau*mu1*mu1, S2 = -((double) n - Stau)*mu2*mu2;
    for(int i = 0; i < n; i++) {
      double xx = x[i]*x[i];
      S1 += tau[i]*xx;
      S2 += (1 - tau[i])*xx;
    }

    if((S1 <= 0) || (S2 <= 0)) break;
    
    double gamma = ((double) n)/(S1 + rho*S2);
    rho = ((double) n - Stau + lambda)/(gamma * S2 + lambda);
    s1 = 1/sqrt(gamma);
    s2 = s1/sqrt(rho);
    
    if( pi != pi || fabs(1 - pi/piOld) < eps )
      break;
  }
  if(verbose)
    Rcout << "Algorithm stopped after " << k << " iterations\n";
  
  NumericVector R = NumericVector::create(_["pi"] = pi, _["mu1"] = mu1, _["mu2"] = mu2, _["s1"] = s1*s1, _["s2"] = s2*s2);
  return R;
}
