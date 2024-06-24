functions {

  //functions are used fom prior work by

  //Dan Schrage (https://gitlab.com/dschrage/rcovreg)

  

  real sum_square_x(matrix x, int i, int j) {

    int j_prime;

    real sum_x = 0;

    if (j == 1) {

      return sum_x;

    }

    

    j_prime = 1;

    while (j_prime < j) {

      sum_x += x[i, j_prime] ^ 2;

      j_prime += 1;

    }

    return sum_x;

  }

  

  matrix lkj_to_chol_corr(row_vector constrained_reals, int D) {

    int z_counter;

    matrix[D, D] x;

    z_counter = 1;

    x[1, 1] = 1;

    for (j in 2 : D) {

      x[1, j] = 0;

    }

    for (i in 2 : D) {

      for (j in 1 : D) {

        if (i == j) {

          x[i, j] = sqrt(1 - sum_square_x(x, i, j));

        } else if (i > j) {

          x[i, j] = constrained_reals[z_counter]

                    * sqrt(1 - sum_square_x(x, i, j));

          z_counter += 1;

        } else {

          x[i, j] = 0;

        }

      }

    }

    return x;

  }

}

data {

  int<lower=0> D;

  int<lower=0> P;

  int<lower=0> nref;

  int<lower=0> nfd;

  int<lower=0> nbs;

  int<lower=0> ngd;

  int<lower=0> nfocal;

  int<lower=0> C;

  int<lower=0> ngroup_fd;

  int<lower=0> ngroup_bs;

  int<lower=0> ngroup_gd;

  int<lower=0> nseason_fd;

  int<lower=0> nseason_bs;

  int<lower=0> nseason_gd;

  int<lower=0> nfxs_fd;

  int<lower=0> nfxs_bs;

  int<lower=0> nfxs_gd;

  

  int<lower=0> cm;

  int<lower=0> cnt;

  array[C] int<lower=0> cn;

  array[nfd] int<lower=0> cfid_fd;

  array[nbs] int<lower=0> cfid_bs;

  array[ngd] int<lower=0> cfid_gd;

  array[C, cm] int cmat;

  

  array[nfd] int<lower=0> focal_fd;

  array[nbs] int<lower=0> focal_bs;

  array[ngd] int<lower=0> focal_gd;

  array[nfd] int<lower=0> context_fd;

  array[nbs] int<lower=0> context_bs;

  array[ngd] int<lower=0> context_gd;

  array[nfd] int<lower=0> group_fd;

  array[nbs] int<lower=0> group_bs;

  array[ngd] int<lower=0> group_gd;

  array[nfd] int<lower=0> season_fd;

  array[nbs] int<lower=0> season_bs;

  array[ngd] int<lower=0> season_gd;

  array[nfd] int<lower=0> fxs_fd;

  array[nbs] int<lower=0> fxs_bs;

  array[ngd] int<lower=0> fxs_gd;

  

  matrix[nfocal, nfocal] A;

  matrix[C, P] X;

  vector[nfd] offset_fd;

  vector[ngd] offset_gd;

  

  array[nfd] int feed;

  array[nbs] int babysit;

  array[nbs] int bs_tot;

  array[ngd] int guard;

}

transformed data {

  matrix[nfocal, nfocal] LA = cholesky_decompose(A);

  int ncor = (D * (D - 1)) ;

  // Compute, thin, and then scale QR decomposition

  matrix[C, P] Q = qr_thin_Q(X) * sqrt(C - 1);

  matrix[P, P] R = qr_thin_R(X) / sqrt(C - 1);

  matrix[P, P] R_inv = inverse(R);

}

parameters {

  matrix[P, D] B_mq;

  matrix[P, D] B_vq;

  matrix[P, ncor] B_cpcq;

  

  real b_offset_fd;

  real b_offset_gd;

  

  matrix[cnt, D] Z_G;

  matrix[nfocal, D] Z_E;

  cholesky_factor_corr[D] L_E;

  vector<lower=0>[D] sd_E;

  

  vector<lower=0>[nref] sd_fd;

  vector<lower=0>[nref] sd_bs;

  vector<lower=0>[nref] sd_gd;

  vector[ngroup_fd] z_group_fd;

  vector[ngroup_bs] z_group_bs;

  vector[ngroup_gd] z_group_gd;

  vector[nseason_fd] z_season_fd;

  vector[nseason_bs] z_season_bs;

  vector[nseason_gd] z_season_gd;

  vector[nfxs_fd] z_fxs_fd;

  vector[nfxs_bs] z_fxs_bs;

  vector[nfxs_gd] z_fxs_gd;

}

model {

  //environmental reaction norms

  //means

  matrix[C, D] mu = Q * B_mq;

  

  //variances

  matrix[C, D] sd_G = sqrt(exp(Q * B_vq));

  

  //correlations (expressed as canonical partial correlations)

  //matrix form is more efficient to code due to size of ncor

  matrix[C, ncor] cpc_G = tanh(Q * B_cpcq);

  

  //scale non-focal univariate random effects

  vector[ngroup_fd] re_group_fd = z_group_fd * sd_fd[1];

  vector[ngroup_bs] re_group_bs = z_group_bs * sd_bs[1];

  vector[ngroup_gd] re_group_gd = z_group_gd * sd_gd[1];

  vector[nseason_fd] re_season_fd = z_season_fd * sd_fd[2];

  vector[nseason_bs] re_season_bs = z_season_bs * sd_bs[2];

  vector[nseason_gd] re_season_gd = z_season_gd * sd_gd[2];

  vector[nfxs_fd] re_fxs_fd = z_fxs_fd * sd_fd[3];

  vector[nfxs_bs] re_fxs_bs = z_fxs_bs * sd_bs[3];

  vector[nfxs_gd] re_fxs_gd = z_fxs_gd * sd_gd[3];

  

  //scale focal multivariate permanent environmental effects

  matrix[nfocal, D] mat_E = Z_E * diag_pre_multiply(sd_E, L_E)';

  

  //initialize mean linear predictors

  vector[nfd] mu_fd = mu[context_fd, 1] + offset_fd * b_offset_fd

                      + col(mat_E, 1)[focal_fd] + re_group_fd[group_fd]

                      + re_season_fd[season_fd] + re_fxs_fd[fxs_fd];

  

  vector[nbs] mu_bs = mu[context_bs, 2] + col(mat_E, 2)[focal_bs]

                      + re_group_bs[group_bs] + re_season_bs[season_bs]

                      + re_fxs_bs[fxs_bs];

  

  vector[nfd] mu_gd = mu[context_gd, 3] + offset_gd * b_offset_gd

                      + col(mat_E, 3)[focal_gd] + re_group_gd[group_gd]

                      + re_season_gd[season_gd] + re_fxs_gd[fxs_gd];

  

  //scale focal context-specific multivariate additive genetic effects

  matrix[cnt, D] mat_G;

  int pos = 1;

  for (c in 1 : C) {

    mat_G[pos : pos + cn[c] - 1] = LA[cmat[c, 1 : cn[c]], cmat[c, 1 : cn[c]]]

                                   * Z_G[pos : (pos + cn[c] - 1)]

                                   * diag_pre_multiply(sd_G[c],

                                                       lkj_to_chol_corr(

                                                       cpc_G[c], D))';

    pos = pos + cn[c];

  }

  

  //add context-specific genetic effects to linear predictors

  for (n in 1 : nfd) {

    mu_fd[n] += col(mat_G, 1)[cfid_fd[n]];

  }

  

  for (n in 1 : nbs) {

    mu_bs[n] += col(mat_G, 2)[cfid_bs[n]];

  }

  

  for (n in 1 : ngd) {

    mu_gd[n] += col(mat_G, 3)[cfid_gd[n]];

  }

  

  //likelihoods

  feed ~ poisson(exp(mu_fd));

  babysit ~ binomial(bs_tot, inv_logit(mu_bs));

  guard ~ poisson(exp(mu_gd));

  

  //priors

  to_vector(B_mq) ~ normal(0, 1);

  to_vector(B_vq) ~ normal(0, 1);

  to_vector(B_cpcq) ~ normal(0, 1);

  b_offset_fd ~ normal(0, 1);

  b_offset_gd ~ normal(0, 1);

  

  to_vector(Z_G) ~ std_normal();

  to_vector(Z_E) ~ std_normal();

  z_group_fd ~ std_normal();

  z_group_bs ~ std_normal();

  z_group_gd ~ std_normal();

  z_season_fd ~ std_normal();

  z_season_bs ~ std_normal();

  z_season_gd ~ std_normal();

  z_fxs_fd ~ std_normal();

  z_fxs_bs ~ std_normal();

  z_fxs_gd ~ std_normal();

  

  sd_E ~ exponential(2);

  sd_fd ~ exponential(2);

  sd_bs ~ exponential(2);

  sd_gd ~ exponential(2);

  L_E ~ lkj_corr_cholesky(2);

}

generated quantities {

  matrix[P, D] B_m;

  matrix[P, D] B_v;

  matrix[P, ncor] B_cpc;

  

  for (i in 1 : D) {

    B_m[ : , i] = R_inv * B_mq[ : , i];

    B_v[ : , i] = R_inv * B_vq[ : , i];

  }

  

  for (i in 1 : ncor) {

    B_cpc[ : , i] = R_inv * B_cpcq[ : , i];

  }

}




