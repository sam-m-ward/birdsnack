data {
  int<lower=1> Nc;
}

parameters {
  real dummy;
}

model {
  dummy ~ std_normal();
}

generated quantities {
  cholesky_factor_corr[Nc] L_cint_eta;
  L_cint_eta     = lkj_corr_cholesky_rng(Nc, 1);
}
