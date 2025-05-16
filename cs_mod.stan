data{
  int<lower=1> I;
  array[I] vector[4] zcs;
  matrix[I,I] A; //relatedness matrix
}

transformed data{
  matrix[I, I] LA = cholesky_decompose(A);
}

parameters{
  //fixed effects
  vector[4] B_m; //intercepts

  //random effects
  matrix[I, 4] Z_G; //all context-specific additive genetic values
  cholesky_factor_corr[4] L_G; //Cholesky corr matrix for g
  vector<lower=0>[4] sd_G; //g standard deviations
  cholesky_factor_corr[4] L_E; //Cholesky corr matrix for residuals
  vector<lower=0>[4] sd_E; //residual standard deviations
}

model{
  //scale context-specific multivariate additive genetic effects
  matrix[I, 4] mat_G = LA * Z_G * diag_pre_multiply(sd_G,L_G)';
                  
  //likelihood
  for(i in 1:I){
  row_vector[4] lin_pred = B_m' + mat_G[i];
  zcs[i] ~ multi_normal_cholesky(lin_pred, diag_pre_multiply(sd_E, L_E)); 
  }
  
  //priors
  B_m ~ normal(0,1);
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


