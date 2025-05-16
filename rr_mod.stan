data{
  int<lower=1> N; //total number of observations
  int<lower=1> I;
  array[N] vector[2] z;
  array[N] int id;
  vector[N] x;
  matrix[I,I] A; //relatedness matrix
}

transformed data{
  matrix[I, I] LA = cholesky_decompose(A);
}

parameters{
  //fixed effects
  vector[2] u; //intercepts
  vector[2] p; //slopes
  

  //random effects
  matrix[I, 4] Z_G; //all context-specific additive genetic values
  cholesky_factor_corr[4] L_G; //Cholesky corr matrix for g
  vector<lower=0>[4] sd_G; //g standard deviations
  cholesky_factor_corr[2] L_E; //Cholesky corr matrix for residuals
  vector<lower=0>[2] sd_E; //residual standard deviations
}

transformed parameters {
  matrix[I, 4] mat_G = LA * Z_G * diag_pre_multiply(sd_G,L_G)';
}


model{
  //scale context-specific multivariate additive genetic effects
  vector[I] u1i = col(mat_G,1);
  vector[I] p1i = col(mat_G,2);
  vector[I] u2i = col(mat_G,3);
  vector[I] p2i = col(mat_G,4);
  //
  
                  
  //likelihood
  for(n in 1:N){
  row_vector[2] lin_pred = [(u[1] + u1i[id[n]] + (p[1] + p1i[id[n]]) * x[n]),
                            (u[2] + u2i[id[n]] + (p[2] + p2i[id[n]]) * x[n])] ;
  z[n] ~ multi_normal_cholesky(lin_pred, diag_pre_multiply(sd_E, L_E)); 
  }
  
  //priors
  u ~ normal(0,1);
  p ~ normal(0,1);
  to_vector(Z_G) ~ std_normal();
  sd_G ~ exponential(2);
  L_G ~ lkj_corr_cholesky(2);
  sd_E ~ exponential(2);
  L_E ~ lkj_corr_cholesky(2);
}

generated quantities{
  matrix[4,4] Gr = L_G * L_G';
  matrix[4,4] G = diag_matrix(sd_G) * Gr * diag_matrix(sd_G);
}


