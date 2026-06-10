rm(list = ls())

library(Rcpp)
library(mclust)
library(mvtnorm)
library(parallel)

Rcpp::sourceCpp("utils_.cpp")
source("utils_.R")

set.seed(123)

TT <- 1000
P  <- 10
K  <- 3
out_frac <- 0.05

rel_ <- list(
  c(1,2,7),
  c(3,4,7),
  c(5,6,7)
)

sim0 <- sim_data_stud_t(
  seed  = 123,
  TT    = TT,
  P     = P,
  Ktrue = K,
  mu    = 2,
  rho   = 0.2,
  nu    = 4,
  pers  = 0.95
)

sim <- simulate_sparse_hmm(
  Y         = sim0$SimData,
  rel_      = rel_,
  true_stat = sim0$mchain,
  perc_out  = out_frac,
  seed      = 123,
  Out_bound = 30
)

Y       <- sim$Y
truth   <- sim0$mchain
W_truth <- sim$W_truth
idx_out <- sim$out_indices

fit <- feat_weight_jump(
  Y       = Y,
  zeta    = 10,
  lambda  = 0.5,
  K       = K,
  n_init  = 5,
  n_iter  = 20,
  tol     = 1e-6,
  verbose = TRUE,
  truth   = truth,
  tukey   = TRUE,
  ncores  = NULL,
  mif=7
)

ARI=adjustedRandIndex(truth, fit$s)
cat("ARI:", round(ARI, 3), "\n")
cat("Final loss:          ", round(fit$loss, 3), "\n")

round(fit$W, 3)

x11()
par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

plot(truth, type = "s", lwd = 1,
     main = "True latent states, with outliers coded as 0",
     xlab = "time", ylab = "state")

plot(fit$s, type = "s", lwd = 1,
     main = "Estimated latent states",
     xlab = "time", ylab = "state")

plot(fit$loss_vec$iter, fit$loss_vec$loss, type = "b",
     main = "Loss path",
     xlab = "iteration", ylab = "loss")

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

image(t(W_truth[K:1, ]), axes = FALSE,
      main = "True relevant features",
      xlab = "state", ylab = "feature")
axis(1, at = seq(0, 1, length.out = K), labels = K:1)
axis(2, at = seq(0, 1, length.out = P), labels = 1:P)

image(t(fit$W[K:1, ]), axes = FALSE,
      main = "Estimated feature weights",
      xlab = "state", ylab = "feature")
axis(1, at = seq(0, 1, length.out = K), labels = K:1)
axis(2, at = seq(0, 1, length.out = P), labels = 1:P)

x11()
plot(truth,
     type = "s",
     lwd = 3,
     ylim = c(0.5, K + 0.5),
     xlab = "Time",
     ylab = "State",
     main = "True vs estimated state sequence")

lines(fit$s,
      type = "s",
      lwd = 2,
      lty = 2,
      col="red")

legend("topright",
       legend = c("True", "Estimated"),
       lwd = c(3,2),
       lty = c(1,2),
       bty = "n",
       col=c("black","red"))
