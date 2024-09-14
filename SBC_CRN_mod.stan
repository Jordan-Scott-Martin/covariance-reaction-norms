functions{

//functions are used fom prior work by
//Dan Schrage (https://gitlab.com/dschrage/rcovreg)
  
  real sum_square_x(matrix x, int i, int j) {
    int j_prime;
    real sum_x = 0;
    if(j==1) return(sum_x);

    j_prime = 1;
    while(j_prime < j) {
      sum_x += x[i,j_prime]^2;
      j_prime += 1;
    }
    return(sum_x);
  }
  
  matrix lkj_to_chol_corr(row_vector constrained_reals, int ntrait) {
    int z_counter;
    matrix[ntrait,ntrait] x;

    z_counter = 1;
    x[1,1] = 1;
    for(j in 2:ntrait) {
      x[1,j] = 0;
    }
    for(i in 2:ntrait) {
      for(j in 1:ntrait) {
        if(i==j) {
          x[i,j] = sqrt(1 - sum_square_x(x, i, j));
        } else if(i > j) {
          x[i,j] = constrained_reals[z_counter]*sqrt(1 - sum_square_x(x, i, j));
          z_counter += 1;
        } else { 
          x[i,j] = 0;
        }
      }
    }
    return(x);
  }
}

data{
  int<lower=1> N;
  int<lower=1> C;
  int<lower=1> I;
  int<lower=1> D;
  int<lower=0> P;
  
  array[N] int<lower=0> id;
  array[N] int<lower=0> c_id;
  array[N] int<lower=0> idc;
  
  matrix[C,P] X;
  matrix[I,I] A;
  
  int<lower=1> cm;
  array[C, cm] int cmat; //array of integers
  array[C] int<lower=0> cn;
  int<lower=1> cnt;
  
  array[N] vector[D] z;
}

transformed data{
  matrix[I, I] LA = cholesky_decompose(A);
  int ncor = (D*(D-1))/2;
}

parameters{
  //fixed effects
  matrix[P, D] B_m; //RN of means
  matrix[P, D] B_v; //RN of variances
  matrix[P, ncor] B_cpc; //RN of canonical partial corrs

  //random effects
  matrix[cnt, D] Z_G; //all context-specific additive genetic values
  cholesky_factor_corr[D] L_E; //Cholesky corr matrix for residuals
  vector<lower=0>[D] sd_E; //residual standard deviations
}

model{
  //predicted values from reaction norms
  //means
  matrix[C, D] mu = X * B_m;
                       
  //variances
  matrix[C, D] sd_G = sqrt(exp(X * B_v));
  
  //correlations (expressed as canonical partial correlations)
  matrix[C, ncor] cpc_G = tanh(X * B_cpc);

  //scale context-specific multivariate additive genetic effects
  matrix[cnt, D] mat_G;
  int pos = 1; //keep track of position 1:cnt
  for(c in 1:C){
      mat_G[pos:(pos+cn[c]-1)] = 
        LA[cmat[c,1:cn[c]],cmat[c,1:cn[c]]] * Z_G[pos:(pos+cn[c]-1)] * 
                diag_pre_multiply(sd_G[c],lkj_to_chol_corr(cpc_G[c], D))';
      pos = pos + cn[c];                             
  }
                  
  //likelihood
  for(n in 1:N){
  row_vector[3] lin_pred = mu[c_id[n]] + mat_G[idc[n]];
  z[n] ~ multi_normal_cholesky(lin_pred, 
                              //Chol factorized lower-tri residual cor matrix
                              diag_pre_multiply(sd_E, L_E)); 
  }
  
  //priors
  to_vector(B_m) ~ normal(0,1);
  to_vector(B_v) ~ normal(0,1);
  to_vector(B_cpc) ~ normal(0,1);
  to_vector(Z_G) ~ std_normal();
  
  sd_E ~ exponential(2);
  L_E ~ lkj_corr_cholesky(2);
}


