// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;
using namespace R;


// utility functions;

//' xlogx
//'
//' utility function
//'
//' @param x a positive number or zero
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
double xlogx(double x){
  double res=0.0;
  if (x!=0.0){
    res = x*log(x);
  }
  return(res);
}

// [[Rcpp::export]]
double logsumexp(arma::vec logv_arma)
{
  //int n = logv.size();
  //if(n<=1)
  //	cout<<"Warning in logsumexp"<<endl;
  double max = logv_arma.max();
  double answer = 0.0;
  // log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
  //			= max(logf) + log(sum(exp(logf - max(logf)))
  answer = max + log(sum(exp(logv_arma-max)));
  return answer;
}

//' logexpit to avoid numerical underflow
//'
//' @param x a number
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
double logexpit_cpp(double x)
{
  arma::vec tmp(2);tmp.zeros();
  tmp(1) = -x;
  return(-logsumexp(tmp));
}


// [[Rcpp::export]]
double logsumexp_row(arma::rowvec logv_arma)
{
  //int n = logv.size();
  //if(n<=1)
  //	cout<<"Warning in logsumexp"<<endl;
  double max = logv_arma.max();
  double answer = 0.0;
  // log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
  //			= max(logf) + log(sum(exp(logf - max(logf)))
  answer = max + log(sum(exp(logv_arma-max)));
  return answer;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//' Calculate variational moments during the updates
//'
//' Get all moments that need update when iterating over a total of p internal and leaf nodes
//'
//' @param prob variational probabilities for \code{s_u}; length p
//' @param prob_gamma should be fixed: \code{c(1,rep(0,p-1))}
//' @param mu_gamma variational Gaussian means (for \code{s_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p) array; In R, we used a list of p (J,K) matrices
//' @param sigma_gamma variational Gaussian variances (for \code{s_u=1} component)
//' for J*K logit(class-specific response probabilities); (J,K,p) array
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_u=1} component) -
//' this is a p by K-1 matrix; in R, we used a list of p vectors (each of length K-1)
//' @param Sigma_alpha variational Gaussian variances (for \code{s_u=1} component)
//' - this is an array of dimension (K-1, K-1, p); in R, we used a list of p matrices,
//' each of dimension K-1 by K-1.
//' @param anc a list of pL vectors, each vector has the node ids of the ancestors;
//' lengths may differ. The ancestors include the node concerned.
//' @param cardanc a numeric vector of length pL; integers. The number
//' of ancestors for each leaf node
//'
//' @return a List
//'
//' \describe{
//'   return List::create(Named("E_beta")=E_beta,
//'    Named("E_beta_sq")=E_beta_sq,
//'    Named("E_eta")=E_eta,
//'    Named("E_eta_sq")=E_eta_sq);
//'}
//'
//' @example
//' inst/example/variance_ss.R
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp(arma::vec prob,
                     arma::vec prob_gamma,
                     arma::cube mu_gamma,//J by K by p
                     arma::cube sigma_gamma,//J by K by p
                     arma::mat mu_alpha, // p by K-1
                     arma::mat Sigma_alpha,//  p by K-1
                     List anc, // length pL; each element is a vector.
                     arma::vec cardanc // length pL; integers
){
  int pL = cardanc.size();
  int J = mu_gamma.n_rows;
  int K = mu_gamma.n_cols;
  int p = mu_gamma.n_slices;
  arma::mat  E_eta(pL,K-1);E_eta.zeros(); // leaf level
  arma::mat  Sigma_xi_u(p,K-1);Sigma_xi_u.zeros();
  arma::mat  E_eta_sq(pL,K-1);E_eta_sq.zeros();// leaf level
  arma::cube E_beta(J,K,pL);E_beta.zeros();// leaf level
  arma::cube sigma_zeta_u(J,K,p);sigma_zeta_u.zeros();
  arma::cube E_beta_sq(J,K,pL);E_beta_sq.zeros(); // leaf level

  int n_anc=0;
  int uu=0;
  for (int v=0;v<pL;v++){
    // E_beta.slice(v)    += mu_gamma.slice(0);
    // E_beta_sq.slice(v) += sigma_gamma.slice(0); //not yet.

    arma::vec curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int u=0;u<n_anc;u++){
      uu =  (int) curr_anc(u)-1;

      E_beta.slice(v)    += prob_gamma(uu)*mu_gamma.slice(uu);
      E_beta_sq.slice(v) += prob_gamma(uu)*(sigma_gamma.slice(uu)+(1.0-prob_gamma(uu))*pow(mu_gamma.slice(uu),2.0)); //not yet.

      E_eta.row(v)        += prob(uu)*mu_alpha.row(uu);
      E_eta_sq.row(v)     += prob(uu)*(Sigma_alpha.row(uu)+(1.0-prob(uu))*pow(mu_alpha.row(uu),2.0));
    }
    E_beta_sq.slice(v)     += pow(E_beta.slice(v),2.0);
    E_eta_sq.row(v)        += pow(E_eta.row(v),2.0);
  }

  // return results:
  return List::create(Named("E_beta")=E_beta,
                      Named("E_beta_sq")=E_beta_sq,
                      Named("E_eta")=E_eta,
                      Named("E_eta_sq")=E_eta_sq
  );
}

//' Calculate variational moments during the updates (only for node u)
//'
//' update only selected moments that need update when iterating over u (except for \code{rmat})
//'
//' (one-node version of \code{\link{get_moments_cpp}})
//' @param leaves_u the leaf descendant node ids for node u
//' @param E_beta,E_beta_sq,E_eta,E_eta_sq moment updates produced by \code{\link{get_moments_cpp}}
//' @param prob variational probabilities for \code{s_u}; length p
//' @param prob_gamma should be fixed: \code{c(1,rep(0,p-1))}
//' @inheritParams get_moments_cpp
//' @return a List
//'
//' \describe{
//'   return List::create(Named("E_beta")=E_beta,
//'    Named("E_beta_sq")=E_beta_sq,
//'    Named("E_eta")=E_eta,
//'    Named("E_eta_sq")=E_eta_sq);
//'}
//'
//' @example
//' inst/example/variance_ss.R
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List get_moments_cpp_eco(arma::vec leaves_u,
                         arma::cube E_beta,// leaf level
                         arma::cube E_beta_sq, // leaf level
                         arma::mat  E_eta, // leaf level
                         arma::mat  E_eta_sq,
                         arma::vec prob,
                         arma::vec prob_gamma,
                         arma::cube mu_gamma,//J by K by p
                         arma::cube sigma_gamma,//J by K by p
                         arma::mat mu_alpha, // p by K-1
                         arma::mat Sigma_alpha,//  p by K-1
                         List anc,
                         arma::vec cardanc
){
  int v = 0; // iterator over the leaves (actual ids in leaves_u, not enumeration)
  int n_anc=0;
  int nv = leaves_u.size();
  int u=0;
  for (int vv=0;vv<nv;vv++){
    v = (int) leaves_u(vv)-1;
    E_beta.slice(v) *=0.0;
    E_eta.row(v)    *= 0.0;
    E_beta_sq.slice(v) *=0.0;
    E_eta_sq.row(v)    *=0.0;

    // E_beta.slice(v)    += mu_gamma.slice(0);
    // E_beta_sq.slice(v) += sigma_gamma.slice(0); //not yet.

    arma::vec curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int uu=0;uu<n_anc;uu++){
      u = (int) curr_anc(uu)-1;
      E_beta.slice(v)    += 0.0+ prob_gamma(u)*mu_gamma.slice(u);
      E_beta_sq.slice(v) += 0.0+prob_gamma(u)*(sigma_gamma.slice(u)+(1.0-prob_gamma(u))*mu_gamma.slice(u)%mu_gamma.slice(u)); //not yet.
      E_eta.row(v)       += 0.0+ prob(u)*mu_alpha.row(u);
      E_eta_sq.row(v)     += 0.0+prob(u)*(Sigma_alpha.row(u)+(1.0-prob(u))*(mu_alpha.row(u)%mu_alpha.row(u)));
    }
    E_beta_sq.slice(v)  += 0.0+E_beta.slice(v)%E_beta.slice(v);
    E_eta_sq.row(v)     += 0.0+E_eta.row(v)%E_eta.row(v);
  }

  // return results:
  return List::create(Named("E_beta")=E_beta,
                      Named("E_beta_sq")=E_beta_sq,
                      Named("E_eta")=E_eta,
                      Named("E_eta_sq")=E_eta_sq
  );
}

//' Summarize the posterior mean, sd and confidence interval (grouped or individual leaf nodes)
//'
//' @param prob a vector of variational probability. Length = p.
//' At the extremes, it can also be a vector of zeros and ones, indicating which nodes
//' are selected based on variational probability (prob > 0.5).
//' @param mu_gamma variational Gaussian means (for \code{s_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p) array; In R, we used a list of p (J,K) matrices
//' @param sigma_gamma variational Gaussian variances (for \code{s_u=1} component)
//' for J*K logit(class-specific response probabilities); (J,K,p) array
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_u=1} component) -
//' this is a p by K-1 matrix; in R, we used a list of p vectors (each of length K-1)
//' @param Sigma_alpha variational Gaussian variances (for \code{s_u=1} component)
//' - this is an array of dimension (K-1, K-1, p); in R, we used a list of p matrices,
//' each of dimension K-1 by K-1.
//' @param anc a list of pL vectors, each vector has the node ids of the ancestors;
//' lengths may differ. The ancestors include the node concerned.
//' @param cardanc a numeric vector of length pL; integers. The number
//' of ancestors for each leaf node
//' @param z  = \code{ ci_level+(1-ci_level)/2}
//'
//' @return a list
//' \describe{
//'
//' Named("beta_est")=beta_est,
//' Named("beta_sd_est")=beta_sd_est,
//' Named("beta_cil")=beta_cil,
//' Named("beta_ciu")=beta_ciu,
//'
//'
//' Named("eta_est")=eta_est,
//' Named("eta_var_est")=eta_var_est,
//' Named("eta_cil")=eta_cil,
//' Named("eta_ciu")=eta_ciu
//'
//' }
//' @export
// [[Rcpp::export]]
List get_est_cpp(arma::vec prob,
                 arma::cube mu_gamma,//J by K by p
                 arma::cube sigma_gamma,//J by K by p
                 arma::mat  mu_alpha, // p by K-1
                 arma::mat Sigma_alpha,//  p by K-1
                 List anc, // length pL; each element is a vector.
                 arma::vec cardanc, // length pL; integers
                 double z
){
  int pL = anc.size();
  int J = mu_gamma.n_rows;
  int K = mu_gamma.n_cols;
  int p = mu_gamma.n_slices;

  arma::cube beta_est(J,K,pL);beta_est.zeros();
  arma::cube beta_sd_est(J,K,pL);beta_sd_est.zeros();
  arma::cube beta_cil(J,K,pL);beta_cil.zeros();
  arma::cube beta_ciu(J,K,pL);beta_ciu.zeros();
  arma::mat eta_est(pL,K-1);eta_est.zeros();
  arma::mat eta_sd_est(pL,K-1);eta_sd_est.zeros();
  arma::mat eta_cil(pL,K-1);eta_cil.zeros();
  arma::mat eta_ciu(pL,K-1);eta_ciu.zeros();

  int n_anc=0;
  int uu=0;
  arma::vec prob_gamma(p);prob_gamma.zeros();
  prob_gamma(0) = 1;
  NumericVector curr_anc(p);
  for (int v=0;v<pL;v++){
    curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int u=0;u<n_anc;u++){
      uu = (int) curr_anc(u)-1;
      beta_est.slice(v)     += 0.0+prob_gamma(uu)*mu_gamma.slice(uu);
      beta_sd_est.slice(v)  += 0.0+prob_gamma(uu)*sigma_gamma.slice(uu);// this is sd, not sq. Var_q[gamma_u |s_u = 1]
      eta_est.row(v)        += 0.0+prob(uu)*mu_alpha.row(uu);
      eta_sd_est.row(v)     += 0.0+prob(uu)*(Sigma_alpha.row(uu)+(1.0-prob(uu))*pow(mu_alpha.row(uu),2.0));
    }
    beta_sd_est.slice(v) = 0.0+pow(beta_sd_est.slice(v),0.5);
    eta_sd_est.row(v)   = 0.0+pow(eta_sd_est.row(v),0.5);
    beta_cil.slice(v) = 0.0+beta_est.slice(v) - z*beta_sd_est.slice(v);
    beta_ciu.slice(v) = 0.0+beta_est.slice(v) + z*beta_sd_est.slice(v);
    eta_cil.row(v) = 0.0+eta_est.row(v) - z*eta_sd_est.row(v);
    eta_ciu.row(v) = 0.0+eta_est.row(v) + z*eta_sd_est.row(v);
  }
  // return results:
  return List::create(Named("beta_est")=beta_est,
                      Named("beta_sd_est")=beta_sd_est,
                      Named("beta_cil")=beta_cil,
                      Named("beta_ciu")=beta_ciu,
                      Named("eta_est")=eta_est,
                      Named("eta_sd_est")=eta_sd_est,
                      Named("eta_cil")=eta_cil,
                      Named("eta_ciu")=eta_ciu
  );
}


//' Update the variational probabilities of each observation in one of K classes
//'
//' This function updates the N by K matrix \code{rmat} in the package
//'
//' @param psi,g_psi,phi,g_phi local variational parameters
//' @param X transformed data: 2Y-1
//' @param E_beta,E_eta,E_beta_sq,E_eta_sq moment updates produced by \code{\link{get_moments_cpp}}
//' @param v_lookup a vector of length equal to the total number of rows in \code{X};
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//' @return  N by K variational multinomial probabilities; row sums are 1s.
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat update_rmat(arma::cube psi, arma::cube g_psi,arma::mat phi, arma::mat g_phi,
                      arma::mat X,
                      arma::cube E_beta,arma::mat E_eta, arma::cube E_beta_sq,arma::mat E_eta_sq,
                      arma::vec v_lookup){
  int n = X.n_rows, J = psi.n_cols, K = psi.n_slices;
  int v = 0;
  X = X+0.0;
  arma::mat res(n,K);res.zeros();
  arma::vec tmp(n);tmp.zeros();
  for (int i=0;i<n;i++){
    v = v_lookup(i)-1;// get leaf id.
    for (int k=0;k<K;k++){
      for (int j=0;j<J;j++){
        res(i,k) += -log(1.0+exp(-psi(v,j,k)))+(1.0*X(i,j)*E_beta(j,k,v)-psi(v,j,k))*0.5-g_psi(v,j,k)*(E_beta_sq(j,k,v)-pow(psi(v,j,k),2.0));
      }
      if (k<1){//first segment
        res(i,k)   += -log(1.0+exp(-phi(v,k)))+(E_eta(v,k)-phi(v,k))*0.5-g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else if (k< K-1){ // not the first, not the last.
        for (int m=0;m<k-1;m++){
          res(i,k)  += -log(1.0+exp(-phi(v,m)))+(-E_eta(v,m)-phi(v,m))*0.5-g_phi(v,m)*(E_eta_sq(v,m)-pow(phi(v,m),2.0));
        }
        res(i,k)  += -log(1.0+exp(-phi(v,k)))+(E_eta(v,k)-phi(v,k))*0.5-g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else{// k==K-1; the last segment:
        for (int s=0;s<K-1;s++){
          res(i,k)  += -log(1.0+exp(-phi(v,s)))+(-E_eta(v,s)-phi(v,s))*0.5-g_phi(v,s)*(E_eta_sq(v,s)-pow(phi(v,s),2.0));
        }
      }
    }
    tmp(i) = logsumexp_row(res.row(i));
  }
  for (int i=0;i<n;i++){
    for (int k=0;k<K;k++){
      res(i,k) = exp(res(i,k)-tmp(i));
    }
  }
  return(res);
}


//' Update the variational probabilities of each observation in one of K classes
//'
//' This function updates the N by K matrix \code{rmat} in the package
//'
//' @param unknown_ids a vector of integers representing subject ids with unkonwn class memberships
//' @param psi,g_psi,phi,g_phi local variational parameters
//' @param X transformed data: 2Y-1
//' @param E_beta,E_eta,E_beta_sq,E_eta_sq moment updates produced by \code{\link{get_moments_cpp}}
//' @param v_lookup a vector of length equal to the total number of rows in \code{X};
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//' @return  N by K variational multinomial probabilities; row sums are 1s.
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat update_rmat_partial(arma::vec unknown_ids,arma::cube psi, arma::cube g_psi,arma::mat phi, arma::mat g_phi,
                              arma::mat X,
                              arma::cube E_beta,arma::mat E_eta, arma::cube E_beta_sq,arma::mat E_eta_sq,
                              arma::vec v_lookup){
  int n = X.n_rows, J = psi.n_cols, K = psi.n_slices;
  int v = 0;
  X = X+0.0;
  arma::mat res(n,K);res.zeros();
  arma::vec tmp(n);tmp.zeros();
  int n_unknown = unknown_ids.size();
  int i_unknown = 0;
  int i = 0;
  for (int i_unknown=0;i_unknown<n_unknown;i_unknown++){
    i = unknown_ids(i_unknown) -1;
    v = v_lookup(i)-1;// get leaf id.
    // bool is_known = std::find(known_ids.begin(), known_ids.end(),i+1)!=known_ids.end();
    //if (is_known){ continue;}
    // if (is_known){res(i,(int) known_Z(i_known,1)-1) = 1.0; ++i_known; continue;}
    for (int k=0;k<K;k++){
      for (int j=0;j<J;j++){
        res(i,k) += -log(1.0+exp(-psi(v,j,k)))+(1.0*X(i,j)*E_beta(j,k,v)-psi(v,j,k))*0.5-g_psi(v,j,k)*(E_beta_sq(j,k,v)-pow(psi(v,j,k),2.0));
      }
      if (k<1){//first segment
        res(i,k)   += -log(1.0+exp(-phi(v,k)))+(E_eta(v,k)-phi(v,k))*0.5-g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else if (k< K-1){ // not the first, not the last.
        for (int m=0;m<k-1;m++){
          res(i,k)  += -log(1.0+exp(-phi(v,m)))+(-E_eta(v,m)-phi(v,m))*0.5-g_phi(v,m)*(E_eta_sq(v,m)-pow(phi(v,m),2.0));
        }
        res(i,k)  += -log(1.0+exp(-phi(v,k)))+(E_eta(v,k)-phi(v,k))*0.5-g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else{// k==K-1; the last segment:
        for (int s=0;s<K-1;s++){
          res(i,k)  += -log(1.0+exp(-phi(v,s)))+(-E_eta(v,s)-phi(v,s))*0.5-g_phi(v,s)*(E_eta_sq(v,s)-pow(phi(v,s),2.0));
        }
      }
    }
    tmp(i) = logsumexp_row(res.row(i));
  }
  for (int i=0;i<n;i++){
    for (int k=0;k<K;k++){
      res(i,k) = exp(res(i,k)-tmp(i)); // all 1s for rows that correspond to known subjects.
    }
  }
  return(res);
}


//' Update gamma and alpha together. Update the variational mean and variance for logit of
//' class-specific response probabilities (for the \code{s_u=1} component)
//'
//' shared tau's
//'
//' @param u node id (internal or leaf node
//' @param g_psi,g_phi g of local variational parameters
//' @param tau_2_t_u,tau_1_t_u variational Gaussian variances for gamma and alpha
//' @param E_beta,E_zeta_u moment updates produced by \code{\link{get_moments_cpp}};
//' \code{E_zeta_u} is directly calculated: \code{prob[u]*sigma_gamma[u,,]}
//' @param X transformed data: 2Y-1
//' @param E_eta leaves' expected eta
//' @param E_xi_u node u's expected xi
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent
//' @param levels a vector of possibly repeating integers from 1 to Fg, or L,
//' @param subject_ids the ids of subjects in the leaf descendants of node u
//' @param v_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//' @return  a list
//' \describe{
//'   \item{resA}{actually 1/A in the paper, this is variance}
//'   \item{resB}{}
//'   \item{logresBsq_o_A}{}
//'   \item{resC}{actually 1/C in the paper, this is variance}
//'   \item{resD}{}
//'   \item{logresDsq_o_C}{}
//' }
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_gamma_alpha_subid(int u,
                              arma::cube g_psi,arma::mat g_phi,
                              double tau_2_t_u,double tau_1_t_u,
                              arma::cube E_beta,arma::mat E_zeta_u,arma::mat X,
                              arma::mat E_eta, arma::vec E_xi_u,
                              arma::mat rmat,arma::vec h_pau,arma::vec levels,
                              arma::vec subject_ids,
                              arma::vec v_lookup
){
  int n = subject_ids.size(), J = g_psi.n_cols, K = g_psi.n_slices;
  int p = h_pau.size(), pL = g_psi.n_rows;
  int ii = 0;
  int vv = 0;
  int uu = (int) u-1;
  X = X+0.0;
  arma::mat resA(J,K);resA.zeros(); // inv A, or variance
  arma::mat resB(J,K);resB.zeros();
  arma::mat logresBsq_o_A(J,K);logresBsq_o_A.zeros();
  arma::vec pre_resA(2);pre_resA.zeros();

  arma::vec resC(K-1);resC.zeros(); // inv C, or variance
  arma::vec resD(K-1);resD.zeros();
  arma::vec logresDsq_o_C(K-1);logresDsq_o_C.zeros();
  arma::vec pre_resC(2);pre_resC.zeros();

  for (int k=0;k<K;k++){
    if (k<(K-1)){ resC(k) = -log(tau_1_t_u)-log(h_pau(uu));}
    for (int j=0;j<J;j++){
      resA(j,k) = -log(tau_2_t_u)-log(h_pau(uu));
      for (int i=0;i<n;i++){
        ii = (int) subject_ids(i)-1;
        vv = (int) v_lookup(ii)-1;
        // update gamma:
        pre_resA(0) = resA(j,k);
        pre_resA(1) = log(2.0) + log(std::min(std::max(rmat(ii,k),1.0e-200),1.0-1e-200)) + log(g_psi(vv,j,k));
        resA(j,k)   = logsumexp(pre_resA);
        // resA(j,k) += 2.0*rmat(ii,k)*g_psi(vv,j,k);
        resB(j,k) += rmat(ii,k)*(X(ii,j)*0.5- 2.0*g_psi(vv,j,k)*(E_beta(j,k,vv)-E_zeta_u(j,k)));

        if (j<1 && k<(K-1)){// only do it once.
          // update alpha:
          for (int m=k;m<K;m++){ // currently does not deal with K=2 case. NB: need fixing.
            pre_resC(0) = resC(k);
            pre_resC(1) = log(2.0)+log(rmat(ii,m))+log(g_phi(vv,k));
            resC(k) =  logsumexp(pre_resC);
            if (m<k+1){
              resD(k) += rmat(ii,m)*0.5-2.0*rmat(ii,m)*g_phi(vv,k)*(E_eta(vv,k)-E_xi_u(k));
            }else{
              resD(k) += -rmat(ii,m)*0.5-2.0*rmat(ii,m)*g_phi(vv,k)*(E_eta(vv,k)-E_xi_u(k));
            }
          }
        }
      }
      logresBsq_o_A(j,k) = 2.0*log(abs(resB(j,k)))-resA(j,k);
      resA(j,k) = exp(-resA(j,k));
      // resA(j,k) += pow(tau_2_t_u*h_pau(uu),-1.0);
    }
    if ( k<(K-1)){// only do it once.
      logresDsq_o_C(k) = 2.0*log(abs(resD(k)))-resC(k);
      resC(k) = exp(-resC(k));
      // resC(k) += pow(tau_1_t_u*h_pau(uu),-1.0);
    }
  }
  return List::create(Named("resA")=resA,
                      Named("resB")=resB,
                      Named("logresBsq_o_A")=logresBsq_o_A,
                      Named("resC")=resC,
                      Named("resD")=resD,
                      Named("logresDsq_o_C")=logresDsq_o_C);

}


//' Update gamma and alpha together. Update the variational mean and variance for logit of
//' class-specific response probabilities (for the \code{s_u=1} component)
//'
//' separate tau's
//'
//' @param u node id (internal or leaf node
//' @param g_psi,g_phi g of local variational parameters
//' @param tau_2_t_u,tau_1_t_u variational Gaussian variances for gamma and alpha
//' @param E_beta,E_zeta_u moment updates produced by \code{\link{get_moments_cpp}};
//' \code{E_zeta_u} is directly calculated: \code{prob[u]*sigma_gamma[u,,]}
//' @param X transformed data: 2Y-1
//' @param E_eta leaves' expected eta
//' @param E_xi_u node u's expected xi
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent
//' @param levels a vector of possibly repeating integers from 1 to Fg, or L,
//' @param subject_ids the ids of subjects in the leaf descendants of node u
//' @param v_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//' @return  a list
//' \describe{
//'   \item{resA}{actually 1/A in the paper, this is variance}
//'   \item{resB}{}
//'   \item{logresBsq_o_A}{}
//'   \item{resC}{actually 1/C in the paper, this is variance}
//'   \item{resD}{}
//'   \item{logresDsq_o_C}{}
//' }
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_gamma_alpha_subid_separate_tau(int u,
                                           arma::cube g_psi,arma::mat g_phi,
                                           arma::mat tau_2_t_u,arma::vec tau_1_t_u,
                                           arma::cube E_beta,arma::mat E_zeta_u,arma::mat X,
                                           arma::mat E_eta, arma::vec E_xi_u,
                                           arma::mat rmat,arma::vec h_pau,arma::vec levels,
                                           arma::vec subject_ids,
                                           arma::vec v_lookup
){
  int n = subject_ids.size(), J = g_psi.n_cols, K = g_psi.n_slices;
  int p = h_pau.size(), pL = g_psi.n_rows;
  int ii = 0;
  int vv = 0;
  int uu = (int) u-1;
  X = X+0.0;
  arma::mat resA(J,K);resA.zeros(); // inv A, or variance
  arma::mat resB(J,K);resB.zeros();
  arma::mat logresBsq_o_A(J,K);logresBsq_o_A.zeros();
  arma::vec pre_resA(2);pre_resA.zeros();

  arma::vec resC(K-1);resC.zeros(); // inv C, or variance
  arma::vec resD(K-1);resD.zeros();
  arma::vec logresDsq_o_C(K-1);logresDsq_o_C.zeros();
  arma::vec pre_resC(2);pre_resC.zeros();

  for (int k=0;k<K;k++){
    if (k<(K-1)){ resC(k) = -log(tau_1_t_u(k))-log(h_pau(uu));}
    for (int j=0;j<J;j++){
      resA(j,k) = -log(tau_2_t_u(j,k))-log(h_pau(uu));
      for (int i=0;i<n;i++){
        ii = (int) subject_ids(i)-1;
        vv = (int) v_lookup(ii)-1;
        // update gamma:
        pre_resA(0) = resA(j,k);
        pre_resA(1) = log(2.0) + log(std::min(std::max(rmat(ii,k),1.0e-208),1.0-1e-208)) + log(g_psi(vv,j,k));
        resA(j,k)   = logsumexp(pre_resA);
        // resA(j,k) += 2.0*rmat(ii,k)*g_psi(vv,j,k);
        resB(j,k) += rmat(ii,k)*(X(ii,j)*0.5- 2.0*g_psi(vv,j,k)*(E_beta(j,k,vv)-E_zeta_u(j,k)));

        if (j<1 && k<(K-1)){// only do it once.
          // update alpha:
          for (int m=k;m<K;m++){ // currently does not deal with K=2 case. NB: need fixing.
            pre_resC(0) = resC(k);
            pre_resC(1) = log(2.0)+log(rmat(ii,m))+log(g_phi(vv,k));
            resC(k) =  logsumexp(pre_resC);
            if (m<k+1){
              resD(k) += rmat(ii,m)*0.5-2.0*rmat(ii,m)*g_phi(vv,k)*(E_eta(vv,k)-E_xi_u(k));
            }else{
              resD(k) += -rmat(ii,m)*0.5-2.0*rmat(ii,m)*g_phi(vv,k)*(E_eta(vv,k)-E_xi_u(k));
            }
          }
        }
      }
      // logresBsq_o_A(j,k) = 2.0*log(abs(resB(j,k)))-resA(j,k);
      // resA(j,k) = exp(-resA(j,k));
      // resA(j,k) += pow(tau_2_t_u*h_pau(uu),-1.0);
    }
    // if ( k<(K-1)){// only do it once.
    //   logresDsq_o_C(k) = 2.0*log(abs(resD(k)))-resC(k);
    //   resC(k) = exp(-resC(k));
    //   // resC(k) += pow(tau_1_t_u*h_pau(uu),-1.0);
    // }
  }

  logresBsq_o_A = 2.0*log(abs(resB))-resA;
  resA = exp(-resA);

  logresDsq_o_C = 2.0*log(abs(resD))-resC;
  resC = exp(-resC);
  return List::create(Named("resA")=resA,
                      Named("resB")=resB,
                      Named("logresBsq_o_A")=logresBsq_o_A,
                      Named("resC")=resC,
                      Named("resD")=resD,
                      Named("logresDsq_o_C")=logresDsq_o_C);

}

//' Initialize sigma alpha for shared tau
//'
//' @param u,g_phi,rmat,tau_1_t,h_pau,subject_ids,v_lookup NB: lazy now, see other functions.
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::rowvec getC(int u, arma::mat g_phi,
                  arma::mat rmat, arma::vec tau_1_t,//K-1
                  arma::vec h_pau,
                  arma::vec subject_ids,
                  arma::vec v_lookup
){
  int n = subject_ids.size(), K = rmat.n_cols;
  int p = tau_1_t.size(), pL = g_phi.n_rows;
  int ii = 0;
  int vv = 0;
  u=(int) u-1;
  arma::rowvec resC(K-1);resC.zeros();
  for (int k=0;k<K-1;k++){
    for (int i=0;i<n;i++){
      ii = (int) subject_ids(i)-1;
      vv = (int) v_lookup(ii)-1;
      for (int m=k;m<K;m++){ // currently does not deal with K=2 case. NB: need fixing.
        resC(k) += 2.0*rmat(ii,m)*g_phi(vv,k);
      }
    }
    resC(k) += pow(tau_1_t(u)*h_pau(u),-1.0);
  }
  return resC;
}

//' Initialize sigma alpha for distinct tau
//'
//' @param u,g_phi,rmat,tau_1_t,h_pau,subject_ids,v_lookup see \code{\link{getC}}
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::rowvec getC_separate_tau(int u, arma::mat g_phi,
                               arma::mat rmat, arma::mat tau_1_t,//p by K-1
                               arma::vec h_pau,
                               arma::vec subject_ids,
                               arma::vec v_lookup
){
  int n = subject_ids.size(), K = rmat.n_cols;
  int p = tau_1_t.n_rows, pL = g_phi.n_rows;
  int ii = 0;
  int vv = 0;
  u=u-1;
  arma::rowvec resC(K-1);resC.zeros();
  for (int k=0;k<K-1;k++){
    for (int i=0;i<n;i++){
      ii = subject_ids(i)-1;
      vv = v_lookup(ii)-1;
      for (int m=k;m<K;m++){ // currently does not deal with K=2 case. NB: need fixing.
        resC(k) += 2.0*rmat(ii,m)*g_phi(vv,k);
      }
    }
    resC(k) += pow(tau_1_t(u,k)*h_pau(u),-1.0);
  }
  return resC;
}

//' calculate line 1 and 2 and 13 of ELBO to assess convergence
//' and choose among converged estimates from many restarts
//'
//' @param psi,g_psi,phi,g_phi see \code{\link{update_hyperparams}}
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param E_beta,E_beta_sq,E_eta,E_eta_sq moments during
//' VI updates from \code{\link{get_moments_cpp}}
//' @param X transformed data: 2Y-1
//' @param v_lookup a vector of indicators; of size N, each indicating the leaf id (from 1 to pL)
//' for each sample.
//'
//' @return line 1, 2, 13 of the ELBO
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List get_line1_2_13_subid(arma::cube psi, arma::cube g_psi,
                          arma::mat phi, arma::mat g_phi,
                          arma::mat rmat,
                          arma::cube E_beta, arma::cube E_beta_sq,
                          arma::mat E_eta, arma::mat E_eta_sq,
                          arma::mat X,
                          arma::vec v_lookup){
  int n = X.n_rows, J = psi.n_cols, K = rmat.n_cols;
  int v = 0;
  double res1 = 0.0;
  double res2 = 0.0;
  double res3 = 0.0;
  X = X+0.0;
  for (int k=0;k<K;k++){
    for (int ii=0;ii<n;ii++){
      v = (int) v_lookup(ii)-1;
      for (int j=0;j<J;j++){
        res1 += -rmat(ii,k)*log(1.0+exp(-psi(v,j,k)))+rmat(ii,k)*(1.0*X(ii,j)*E_beta(j,k,v)-psi(v,j,k))*0.5-rmat(ii,k)*g_psi(v,j,k)*(E_beta_sq(j,k,v)-pow(psi(v,j,k),2.0));
      }
      if (k<1){//first segment
        res2   += -rmat(ii,k)*log(1.0+exp(-phi(v,k)))+rmat(ii,k)*(E_eta(v,k)-phi(v,k))*0.5-rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else if (k< K-1){ // not the first, not the last.
        for (int m=0;m<k-1;m++){
          res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,m)))+rmat(ii,k)*(-E_eta(v,m)-phi(v,m))*0.5-rmat(ii,k)*g_phi(v,m)*(E_eta_sq(v,m)-pow(phi(v,m),2.0));
        }
        res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,k)))+rmat(ii,k)*(E_eta(v,k)-phi(v,k))*0.5-rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
      } else{// k==K-1; the last segment:
        for (int s=0;s<K-1;s++){
          res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,s)))+rmat(ii,k)*(-E_eta(v,s)-phi(v,s))*0.5-rmat(ii,k)*g_phi(v,s)*(E_eta_sq(v,s)-pow(phi(v,s),2.0));
        }
      }
      res3 += - xlogx(rmat(ii,k)); // line 13
    }
  }
  return List::create(Named("res1")=res1,
                      Named("res2")=res2,
                      Named("res3")=res3);
}



