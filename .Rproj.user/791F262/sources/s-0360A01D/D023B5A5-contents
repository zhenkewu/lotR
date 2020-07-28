// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;
using namespace R;

//' Get all moments that need update when iterating over a total of p internal and leaf nodes
//'
//' This is slower than the skinny version (\code{\link{get_moments_cpp}})
//'
//' @param prob variational probabilities for \code{s_u}; length p
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
//' @param D_k_array a utility array of dimension (K-1,K,K); for each slice,
//' multiplying it with a vector of length K produces \code{(alpha_k-alpha_m), for m not k}
//' @param idnotk K-1 by K, utility matrix, column k indicates the indices that are not k
//' @param submat K-1 by K, utility matrix with a K-1-identity, right padded with
//' a column of zeros.
//'
//' @return a List
//'
//' \describe{
//'     Named("E_beta")=E_beta,
//'     Named("E_zeta")=E_zeta,
//'     Named("E_beta_sq")=E_beta_sq,
//'     Named("E_xi")=E_xi,
//'     Named("E_eta")=E_eta,
//'     Named("E_eta_sq")=E_eta_sq
//' }
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp(arma::vec prob,
                     arma::cube mu_gamma,//J by K by p
                     arma::cube sigma_gamma,//J by K by p
                     arma::mat mu_alpha, // p by K-1
                     arma::mat Sigma_alpha,//  p by K-1
                     List anc, // length pL; each element is a vector.
                     arma::vec cardanc // length pL; integers
){
  int pL = cardanc.size();
  int J = sigma_gamma.n_rows;
  int K = sigma_gamma.n_cols;
  int p = sigma_gamma.n_slices;
  arma::cube E_zeta(J,K,p);E_zeta.zeros();
  arma::cube sigma_zeta_u(J,K,p);sigma_zeta_u.zeros();
  arma::cube E_beta(J,K,pL);E_beta.zeros();// leaf level
  arma::cube E_beta_sq(J,K,pL);E_beta_sq.zeros(); // leaf level

  arma::mat  E_xi(p,K-1);E_xi.zeros();
  arma::mat  Sigma_xi_u(p,K-1);Sigma_xi_u.zeros();
  arma::mat  E_eta(pL,K-1);E_eta.zeros(); // leaf level
  arma::mat  E_eta_sq(pL,K-1);E_eta_sq.zeros();// leaf level

  int n_anc=0;
  int uu=0;
  for (int u=0;u<p;u++){
    // response probabilities:
    E_zeta.slice(u) =  0.0+prob(u)*mu_gamma.slice(u);
    sigma_zeta_u.slice(u) = 0.0+prob(u)*(sigma_gamma.slice(u)+(1.0-prob(u))*pow(mu_gamma.slice(u),2.0));

    // class probabilities:
    E_xi.row(u)           = 0.0+prob(u)*mu_alpha.row(u);
    Sigma_xi_u.row(u)     = 0.0+prob(u)*(Sigma_alpha.row(u)+(1.0-prob(u))*pow(mu_alpha.row(u),2.0));
  }
  for (int v=0;v<pL;v++){
    arma::vec curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int u=0;u<n_anc;u++){
      uu =  (int) curr_anc(u)-1;
      // Rcout << "ancectors: " << uu << std::endl;
      E_beta.slice(v)    += 0.0+E_zeta.slice(uu);
      E_beta_sq.slice(v) += 0.0+sigma_zeta_u.slice(uu); //not yet.

      E_eta.row(v)        += 0.0+E_xi.row(uu);
      E_eta_sq.row(v)     += 0.0+Sigma_xi_u.row(uu);
    }
    E_beta_sq.slice(v)     += 0.0+pow(E_beta.slice(v),2.0);
    E_eta_sq.row(v)        += 0.0+pow(E_eta.row(v),2.0);
  }

  // return results:
  return List::create(Named("E_beta")=E_beta,
                      Named("E_zeta")=E_zeta,
                      Named("E_beta_sq")=E_beta_sq,
                      Named("E_xi")=E_xi,
                      Named("E_eta")=E_eta,
                      Named("E_eta_sq")=E_eta_sq
  );
}


//' get moments that need update when iterating over u (except for \code{rmat})
//'
//'
//' (skinny version of \code{\link{get_moments_cpp}})
//'
//'
//' @param prob variational probabilities for \code{s_u}; length p
//' @param mu_gamma variational Gaussian means (for \code{s_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p) array; In R, we used a list of p (J,K) matrices
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_u=1} component) -
//' this is a p by K-1 matrix; in R, we used a list of p vectors (each of length K-1)
//' @param anc a list of pL vectors, each vector has the node ids of the ancestors;
//' lengths may differ. The ancestors include the node concerned.
//' @param cardanc a numeric vector of length pL; integers. The number
//' of ancestors for each leaf node
//'
//' @return a List
//'
//' \describe{
//' Named("E_beta")=E_beta,
//' Named("E_zeta")=E_zeta,
//' # Named("E_beta_sq")=E_beta_sq,
//' Named("E_xi")=E_xi,
//' Named("E_eta")=E_eta,
//' Named("E_xi_diff")=E_xi_diff,
//' Named("E_eta_diff")=E_eta_diff,
//' # Named("E_eta_diff_sq")=E_eta_diff_sq
//' }
//'
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp_skinny(arma::vec prob,
                            arma::cube mu_gamma,//J by K by p
                            arma::mat mu_alpha, // p by K-1
                            List anc, // length pL; each element is a vector.
                            arma::vec cardanc // length pL; integers
){
  int pL = anc.size();
  int J = mu_gamma.n_rows;
  int K = mu_gamma.n_cols;
  int p = mu_gamma.n_slices;
  arma::cube E_zeta(J,K,p);E_zeta.zeros();
  arma::cube E_beta(J,K,pL);E_beta.zeros();// leaf level

  arma::mat  E_xi(p,K-1);E_xi.zeros();
  arma::mat  E_eta(pL,K-1);E_eta.zeros(); // leaf level

  int n_anc=0;
  int uu=0;
  NumericVector curr_anc(p);
  for (int u=0;u<p;u++){
    E_zeta.slice(u) = 0.0+prob(u)*mu_gamma.slice(u);
    E_xi.row(u) = 0.0+prob(u)*mu_alpha.row(u);
  }
  for (int v=0;v<pL;v++){
    curr_anc = anc[v];
    n_anc = cardanc(v);
    for (int u=0;u<n_anc;u++){
      uu = (int) curr_anc(u)-1;
      E_beta.slice(v)     += 0.0+E_zeta.slice(uu);
      E_eta.row(v)        += 0.0+E_xi.row(uu);
    }
  }

  // return results:
  return List::create(Named("E_beta")=E_beta,
                      Named("E_zeta")=E_zeta,
                      Named("E_xi")=E_xi,
                      Named("E_eta")=E_eta);
}


//' update only selected moments that need update when iterating over u (except for \code{rmat})
//'
//'
//' (skinny version of \code{\link{get_moments_cpp}})
//'
//'
//' @param prob variational probabilities for \code{s_u}; length p
//' @param mu_gamma variational Gaussian means (for \code{s_u=1} component) for J*K
//' logit(class-specific response probabilities); (J,K,p) array; In R, we used a list of p (J,K) matrices
//' @param mu_alpha variational Gaussian mean vectors (for \code{s_u=1} component) -
//' this is a p by K-1 matrix; in R, we used a list of p vectors (each of length K-1)
//' @param anc a list of pL vectors, each vector has the node ids of the ancestors;
//' lengths may differ. The ancestors include the node concerned.
//' @param cardanc a numeric vector of length pL; integers. The number
//' of ancestors for each leaf node
//'
//' @return a List
//'
//' \describe{
//' Named("E_beta")=E_beta,
//' Named("E_zeta")=E_zeta,
//' # Named("E_beta_sq")=E_beta_sq,
//' Named("E_xi")=E_xi,
//' Named("E_eta")=E_eta,
//' Named("E_xi_diff")=E_xi_diff,
//' Named("E_eta_diff")=E_eta_diff,
//' # Named("E_eta_diff_sq")=E_eta_diff_sq
//' }
//'
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_moments_cpp_eco(arma::vec leaves_u,
                         arma::cube E_beta,// leaf level
                         arma::cube E_beta_sq, // leaf level
                         arma::mat  E_eta, // leaf level
                         arma::mat  E_eta_sq,
                         arma::vec prob,
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
    arma::vec curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int uu=0;uu<n_anc;uu++){
      u = (int) curr_anc(uu)-1;
      E_beta.slice(v)    += 0.0+ prob(u)*mu_gamma.slice(u);
      E_eta.row(v)       += 0.0+ prob(u)*mu_alpha.row(u);
      E_beta_sq.slice(v) += 0.0+prob(u)*(sigma_gamma.slice(u)+(1.0-prob(u))*mu_gamma.slice(u)%mu_gamma.slice(u)); //not yet.
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


//' Summarize the posterior mean, sd and confidence interval
//'
//' @param node_select a vector of zeros and ones, indicating which nodes
//' are selected based on variational probability (prob > 0.5). Length = p
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
//' @param z  = \code{ ci_level+(1-ci_level)/2}
//'
//' @return a list
//' \describe{
//'
//' Named("beta_est")=beta_est,
//' Named("beta_sd_est")=beta_sd_est,
//' Named("beta_cil")=beta_cil,
//' Named("beta_ciu")=beta_ciu,
//' Named("eta_est")=eta_est,
//' Named("eta_var_est")=eta_var_est,
//' Named("eta_cil")=eta_cil,
//' Named("eta_ciu")=eta_ciu
//'
//'
//' }

// [[Rcpp::export]]
List get_est_cpp(arma::vec node_select,
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
  NumericVector curr_anc(p);
  for (int v=0;v<pL;v++){
    curr_anc = anc[v];
    n_anc = (int) cardanc(v);
    for (int u=0;u<n_anc;u++){
      uu = (int) curr_anc(u)-1;
      beta_est.slice(v)     += 0.0+node_select(uu)*mu_gamma.slice(uu);
      beta_sd_est.slice(v)  += 0.0+node_select(uu)*sigma_gamma.slice(uu);
      eta_est.row(v)        += 0.0+node_select(uu)*mu_alpha.row(uu);
      eta_sd_est.row(v)     += 0.0+node_select(uu)*Sigma_alpha.row(uu);
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

//' xlogx
//'
//' utility function
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
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
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


//' Z rmat update
//'
//' @param psi, g_psi, phi, g_phi local variational parameters
//' @param X transformed data: 2Y-1
//' @param E_beta, E_eta_diff, E_beta_sq, E_eta_diff_sq
//' intermediate moment updates produced by \code{\link{get_moments_cpp}}
//' @param v_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to.
//' @param idnotk utility matrix; column k indicates the indices that are not k
//'
//' @return  n by K variational multinomial probabilities; row sums are 1s.
//'
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::mat update_rZ(arma::cube psi, arma::cube g_psi,
                    arma::mat phi, arma::mat g_phi,
                    arma::mat X, arma::cube E_beta,
                    arma::mat E_eta, arma::cube E_beta_sq,
                    arma::mat E_eta_sq,arma::vec v_lookup){
  int n = X.n_rows, J = psi.n_cols, K = psi.n_slices;
  int v = 0;
  X = X+0.0;
  arma::mat res(n,K);res.zeros();
  arma::vec tmp(n);tmp.zeros();
  //arma::vec incre(2);incre.zeros();
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
  // for (int i=0;i<n;i++){
  //   incre(0) = log(1e-10)+tmp(i); // log scale
  //   for (int k=0;k<K;k++){
  //   incre(1) = res(i,k); // log scale
  //     res(i,k) = logsumexp(incre);
  //   }
  // }
  for (int i=0;i<n;i++){
    //tmp(i) = logsumexp_row(res.row(i));
    for (int k=0;k<K;k++){
      res(i,k) = exp(res(i,k)-tmp(i));
    }
  }
  return(res);
}



//' Update the variational mean and variance for logit of
//' class-specific response probabilities (for the \code{s_u=1} component)
//'
//' @param u node id
//' @param psi a (pL,J,K) array of local variational parameters for approximating
//' the \code{expit(X^v_ij * beta^v_jk)}
//' @param g_psi g transformed psi
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param tau_2_t variational Gaussian variances for gamma; array, dimension:(J,K,p)
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent
//' @param levels a vector of possibly repeating integers from 1 to Fg, or L,
//' where L indicates the number of node subsets, each of which share the hyperparameters
//' \code{rho} - the prior probability of \code{s_u=1}, \code{tau_1[l,]} - K-1, \code{tau_2[l,,]}, J by K,
//' the prior variances for alpha and gamma.
//' @param E_beta, E_zeta the moments obtained from \code{\link{get_moments_cpp}};
//' E_beta are sums of ancestral elements of E_zeta
//' @param X transformed data: 2Y-1
//' @param subject_ids the ids of subjects in the leaf descendants of node u
//' @param v_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//'
//' @return  a list
//' \describe{
//' Named("resMu")=resMu, J by K
//'
//' Named("resSigma")=resSigma, J by K
//'
//' }
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_gamma_subid(int u,
                        arma::cube g_psi,
                        arma::mat rmat,
                        double tau_2_t_u,
                        arma::vec h_pau,
                        arma::vec levels,
                        arma::cube E_beta,
                        arma::mat E_zeta_u, // J by K
                        arma::mat X,
                        arma::vec subject_ids, //subjects that are in the one of leaf descendants of node u.
                        arma::vec v_lookup){
  int n = subject_ids.size(), J = g_psi.n_cols, K = g_psi.n_slices;
  int p = h_pau.size(), pL = g_psi.n_rows;
  int ii = 0;
  int vv = 0;
  u = (int) u-1;
  X = X+0.0;
  arma::mat resA(J,K);resA.zeros();
  arma::mat resB(J,K);resB.zeros();
  arma::mat resBsqA(J,K);resBsqA.zeros();
  arma::vec pre_resA(2);pre_resA.zeros();
  for (int j=0;j<J;j++){
    for (int k=0;k<K;k++){
      resA(j,k) = -log(tau_2_t_u)-log(h_pau(u));
      for (int i=0;i<n;i++){
        ii = (int) subject_ids(i)-1;
        vv = (int) v_lookup(ii)-1;
        pre_resA(0) = resA(j,k);
        pre_resA(1) = log(2.0) + log(rmat(ii,k)) + log(g_psi(vv,j,k));
        resA(j,k)  = logsumexp(pre_resA);
        //resA(j,k) += 2.0*rmat(ii,k)*g_psi(vv,j,k);

        resB(j,k) += rmat(ii,k)*(X(ii,j)*0.5- 2.0*g_psi(vv,j,k)*(E_beta(j,k,vv)-E_zeta_u(j,k)));
        //Rcout<<E_beta(j,k,vv)-E_zeta(j,k,u)<<std::endl;
      }
      resBsqA(j,k) = 2.0*log(abs(resB(j,k)))-resA(j,k);
      resA(j,k) = exp(-resA(j,k));
      //resA(j,k) += pow(tau_2_t_u(j,k)*h_pau(u),-1.0);
    }
  }
  return List::create(Named("resA")=resA,
                      Named("resB")=resB,
                      Named("resBsqA")=resBsqA);
}

//' Update the variational mean and variance for logit of
//' class-specific response probabilities (for the \code{s_u=1} component)
//'
//' @param u node id
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param tau_1_t variational Gaussian variances for alpha; matrix, dimension:(p by K-1)
//' @param h_pau a numeric vector of length p indicating the branch length
//' between a node and its parent
//' @param levels a vector of possibly repeating integers from 1 to Fg, or L,
//' @param E_eta, E_xi the moments obtained from \code{\link{get_moments_cpp}};
//' E_eta are sums of ancestral elements of E_xi
//' @param X transformed data: 2Y-1
//' @param D_k_array, GvK_array, I_tilde,submat,idnotk utility quantities
//' @param subject_ids the ids of subjects in the leaf descendants of node u
//' @param v_lookup a vector of length equal to the total number of rows in X;
//' each element is an integer, indicating which leaf does the observation belong to.
//'
//'
//' @return  a list
//' \describe{
//' Named("resMu")=resMu, K-1 column vector
//'
//' Named("resSigma")=resSigma, K-1 by K-1 matrix
//'
//' }
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List update_alpha_subid(int u, arma::mat g_phi,
                        arma::mat rmat, double tau_1_t_u,//K-1
                        arma::vec h_pau,
                        arma::vec levels,
                        arma::mat E_eta, arma::vec E_xi_u,
                        arma::vec subject_ids,
                        arma::vec v_lookup
){
  int n = subject_ids.size(), K = rmat.n_cols;
  int p = h_pau.size(), pL = g_phi.n_rows;
  arma::vec pre_resC(2);pre_resC.zeros();
  int ii = 0;
  int vv = 0;
  int uu= (int) u-1; // the index in Rcpp is one smaller than the actual index.
  arma::vec resC(K-1);resC.zeros();
  arma::vec resD(K-1);resD.zeros();
  arma::vec resDsqC(K-1);resDsqC.zeros();

  for (int k=0;k<K-1;k++){ // iterate over K-1 alpha's.
    resC(k) = -log(tau_1_t_u)-log(h_pau(uu));
    for (int i=0;i<n;i++){
      ii = (int) subject_ids(i)-1;
      vv = (int) v_lookup(ii)-1;
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
    resDsqC(k) = 2.0*log(abs(resD(k)))-resC(k);
    resC(k) = exp(-resC(k));
    //resC(k) += pow(tau_1_t_u(k)*h_pau(uu),-1.0);
  }
  return List::create(Named("resC")=resC,
                      Named("resDsqC")=resDsqC,
                      Named("resD")=resD);
}

//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
arma::rowvec getC(int u, arma::mat g_phi,
                  arma::mat rmat, arma::vec tau_1_t,//p by K-1
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



//' calculate line 1 and 2 and 13 of ELBO to assess convergence
//' and choose among converged estimates from many restarts
//'
//' @param psi, g_psi see \code{\link{update_gamma_subid}}
//' @param phi a (pL,K,K) array of local variational parameters for approximating
//' the \code{expit(eta_vk-eta_vm)}
//' @param g_phi g transformed phi
//' @param rmat a matrix of variational probabilities of all observations
//' belong to K classes; N by K; each row sums to 1
//' @param outcomes_units a list of length pL, each is a vector indicating the
//' observations in each leaf (aka outcome)
//' @param cardleaf a vector of legnth pL; counts the number of observations in each leaf
//' @param E_beta, E_beta_sq, E_beta_diff, E_beta_diff_sq moments during
//' VI updates from \code{\link{get_moments_cpp}}
//' @param X transformed data: 2Y-1
//' @param idnotk utility matrix
//'
//' @return ELBO value; negative
//'
//' @useDynLib lotR
//' @importFrom Rcpp sourceCpp
//' @export
// [[Rcpp::export]]
List get_line1_2_subid(arma::cube psi, arma::cube g_psi,
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
        // Rcout<<"term in res1 "<<rmat(ii,k)*g_psi(v,j,k)*(E_beta_sq(j,k,v)-pow(psi(v,j,k),2.0))<<std::endl;
        res1 += -rmat(ii,k)*log(1.0+exp(-psi(v,j,k)))+rmat(ii,k)*(1.0*X(ii,j)*E_beta(j,k,v)-psi(v,j,k))*0.5-rmat(ii,k)*g_psi(v,j,k)*(E_beta_sq(j,k,v)-pow(psi(v,j,k),2.0));
      }
      if (k<1){//first segment
        res2   += -rmat(ii,k)*log(1.0+exp(-phi(v,k)))+rmat(ii,k)*(E_eta(v,k)-phi(v,k))*0.5-rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
        // Rcout<<"term in res2 "<<rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0))<<std::endl;
      } else if (k< K-1){ // not the first, not the last.
        for (int m=0;m<k-1;m++){
          res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,m)))+rmat(ii,k)*(-E_eta(v,m)-phi(v,m))*0.5-rmat(ii,k)*g_phi(v,m)*(E_eta_sq(v,m)-pow(phi(v,m),2.0));
        // Rcout<<"term in res2 "<<rmat(ii,k)*g_phi(v,m)*(E_eta_sq(v,m)-pow(phi(v,m),2.0))<<std::endl;

        }
        res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,k)))+rmat(ii,k)*(E_eta(v,k)-phi(v,k))*0.5-rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0));
        // Rcout<<"term in res2 "<<rmat(ii,k)*g_phi(v,k)*(E_eta_sq(v,k)-pow(phi(v,k),2.0))<<std::endl;
      } else{// k==K-1; the last segment:
        for (int s=0;s<K-1;s++){
          res2  += -rmat(ii,k)*log(1.0+exp(-phi(v,s)))+rmat(ii,k)*(-E_eta(v,s)-phi(v,s))*0.5-rmat(ii,k)*g_phi(v,s)*(E_eta_sq(v,s)-pow(phi(v,s),2.0));
          // Rcout<<"term in res2 "<<rmat(ii,k)*g_phi(v,s)*(E_eta_sq(v,s)-pow(phi(v,s),2.0))<<std::endl;
        }
      }
      res3 += - xlogx(rmat(ii,k)); // line 13
    }
  }
  return List::create(Named("res1")=res1,
                      Named("res2")=res2,
                      Named("res3")=res3);
}

// //' calculate line 8
// //' @return ELBO value; negative
// //'
// //' @useDynLib lotR
// //' @importFrom Rcpp sourceCpp
// //' @export
// // [[Rcpp::export]]
// List get_line7_8(arma::vec prob, List tau_2_t, arma::vec h_pau,arma::cube sigma_gamma,int J, int K){
//   double res7 = (J+0.0)*(K+0.0)*sum(prob)*0.5;
//   double res8 = (J+0.0)*(K+0.0)*sum(1-prob)*0.5;
//   int p = prob.size();
//   static const double pi = 3.14159265;
//   for (int u=0; u<p; u++){
//     arma::mat curr_tau_2_t = tau_2_t[u];
//     for (int j=0; j<J; j++){
//       for (int k=0; k<K; k++){
//         res7 += prob(u)*log(2*pi*sigma_gamma(u,j,k))*0.5;
//         res8 += (1-prob(u))*log(2*pi*curr_tau_2_t(j,k))*0.5;
//       }
//     }
//   }
//   return List::create(Named("res7")=res7,
//                       Named("res8")=res8);
// }

//get_line7_8(J*K*sum(prob)*(1+log(2*pi))+sum(sapply(1:p,function(u) sum(prob[u]*log(sigma_gamma[u,,])))))/2

// #J*K*sum(1-prob)/2 + sum(mapply(FUN=function(pp,mat,hh){sum(pp*log(2*pi*mat*hh))},pp=1-prob,mat=tau_2_t,hh=h_pau))/2


//
// // [[Rcpp::export]]
// double gcpp(double eta){
//   double res;
//   if (abs(eta)<1e-6) {
//     res = 0.125;
//   } else{
//     res = (1/(2*eta))*(1/(1+exp(-eta))-0.5);
//   }
//   return res;
// }
//
