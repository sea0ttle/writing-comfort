// ============================================================
// hb_maxdiff.stan
//
// Hierarchical Bayes model for MaxDiff / best-worst scaling.
// Each respondent j has individual-level utilities u[j, k]
// for K attributes, derived from population-level means mu
// and standard deviations sigma via a non-centered
// parameterization.
//
// For each task t, the respondent chose attribute best[t]
// as most important and worst[t] as least important, with
// mid[t] as the middle (unchosen) option. The likelihood
// factorizes as:
//   P(best | best, mid, worst) * P(worst | mid, worst)
// using conditional logit structure.
// ============================================================

data {
  int<lower=1> J;                          // number of respondents
  int<lower=1> K;                          // number of attributes
  int<lower=1> T;                          // number of tasks (all respondents combined)
  array[T] int<lower=1,upper=K> best;      // best-chosen attribute index
  array[T] int<lower=1,upper=K> worst;     // worst-chosen attribute index
  array[T] int<lower=1,upper=K> mid;       // middle (unchosen) attribute index
  array[T] int<lower=1,upper=J> id;        // respondent index for each task
}

parameters {
  matrix[J, K] z;          // non-centered individual deviations
  vector[K] mu;            // population-level mean utilities
  vector<lower=0>[K] sigma; // population-level standard deviations
}

transformed parameters {
  matrix[J, K] u_raw;
  matrix[J, K] u;
  row_vector[K] mu_r = to_row_vector(mu);
  row_vector[K] sd_r = to_row_vector(sigma);
  for (j in 1:J) {
    u_raw[j] = mu_r + z[j] .* sd_r;
    u[j] = u_raw[j] - rep_row_vector(mean(u_raw[j]), K);
  }
}

model {
  to_vector(z) ~ normal(0, 1);
  mu ~ normal(0, 3);
  sigma ~ normal(0, 1.5);
  for (t in 1:T) {
    int j = id[t];
    real denom2 = log_sum_exp(u[j, mid[t]], u[j, worst[t]]);
    real denom1 = log_sum_exp(u[j, best[t]], denom2);
    target += u[j, best[t]] - denom1
            + u[j, worst[t]] - denom2;
  }
}
