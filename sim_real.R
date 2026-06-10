rm(list = ls())

library(Rcpp)
library(parallel)
library(cluster)
library(quantmod)
library(xts)
library(zoo)

Rcpp::sourceCpp("utils_.cpp")
source("utils_.R")

set.seed(123)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

# DATA --------------------------------------------------------------------

assets <- c(
  SPY = "SPY",   # US equity
  EFA = "EFA",   # developed ex-US equity
  EEM = "EEM",   # emerging markets equity
  TLT = "TLT",   # long-term US Treasuries
  IEF = "IEF",   # intermediate US Treasuries
  LQD = "LQD",   # investment grade corporate bonds
  HYG = "HYG",   # high-yield bonds
  GLD = "GLD",   # gold
  DBC = "DBC",   # broad commodities
  UUP = "UUP"    # ETF Invesco DB US Dollar Index Bullish Fund, increases when the US dollar gets stronger
)

getSymbols(
  Symbols = assets,
  src = "yahoo",
  from = "2020-01-01",
  to = Sys.Date(),
  auto.assign = TRUE
)

prices <- do.call(
  merge,
  lapply(assets, function(sym) Ad(get(sym)))
)

colnames(prices) <- names(assets)

Y_xts <- diff(log(prices))
Y_xts <- na.omit(Y_xts)

Y <- as.data.frame(scale(coredata(Y_xts)))
colnames(Y) <- colnames(Y_xts)

dates <- index(Y_xts)

TT <- nrow(Y)
P  <- ncol(Y)

# GRID SELECTION ----------------------------------------------------------

zeta0s  <- c(0.1, 1, 5, 10, 50, 100, 200)
lambdas <- c(0, 0.25, 0.5, 0.75, 1)
Ks      <- 2:5

ex <- expand.grid(
  K      = Ks,
  zeta   = zeta0s,
  lambda = lambdas,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

di <- dist(Y)

parallel_try <- TRUE
n_cores <- 8

if (parallel_try) {
  
  sls <- mclapply(seq_len(nrow(ex)), function(i) {
    
    fit_i <- feat_weight_jump(
      Y       = Y,
      zeta    = ex$zeta[i],
      lambda  = ex$lambda[i],
      K       = ex$K[i],
      tol     = 1e-16,
      n_init  = 5,
      n_iter  = 10,
      verbose = FALSE,
      mif     = 1,
      tukey   = TRUE
    )
    
    median(silhouette(fit_i$s, di)[, 3])
    
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  
} else {
  
  sls <- lapply(seq_len(nrow(ex)), function(i) {
    
    fit_i <- feat_weight_jump(
      Y       = Y,
      zeta    = ex$zeta[i],
      lambda  = ex$lambda[i],
      K       = ex$K[i],
      tol     = 1e-16,
      n_init  = 5,
      n_iter  = 10,
      verbose = FALSE,
      mif     = 1,
      tukey   = TRUE
    )
    
    median(silhouette(fit_i$s, di)[, 3])
  })
}

ex$silhouette <- unlist(sls)

best_row <- ex[which.max(ex$silhouette), ]

cat("\nSelected hyperparameters\n")
print(best_row)

# FINAL FIT ---------------------------------------------------------------

fit_real <- feat_weight_jump(
  Y       = Y,
  zeta    = best_row$zeta,
  lambda  = best_row$lambda,
  K       = best_row$K,
  tol     = 1e-16,
  n_init  = 10,
  n_iter  = 50,
  verbose = TRUE,
  mif     = 1,
  tukey   = TRUE
)

cat("\nFinal loss:", round(fit_real$loss, 3), "\n")

# RESULTS -----------------------------------------------------------------

res_ <- data.frame(
  date    = dates,
  Y,
  cluster = fit_real$s
)

state_means <- aggregate(. ~ cluster, data = res_[, c(colnames(Y), "cluster")], mean)
state_sds   <- aggregate(. ~ cluster, data = res_[, c(colnames(Y), "cluster")], sd)

cat("\nState means\n")
print(round(state_means, 3))

cat("\nState standard deviations\n")
print(round(state_sds, 3))

est_feat <- fit_real$W
colnames(est_feat) <- colnames(Y)
rownames(est_feat) <- paste0("State_", seq_len(nrow(est_feat)))

cat("\nEstimated feature weights\n")
print(round(est_feat, 3))

# PLOTS -------------------------------------------------------------------

x11()
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

plot(dates,
     fit_real$s,
     type = "s",
     lwd = 2,
     xlab = "Time",
     ylab = "State",
     main = "Estimated state sequence")

plot(fit_real$loss_vec$iter,
     fit_real$loss_vec$loss,
     type = "b",
     xlab = "Iteration",
     ylab = "Loss",
     main = "Loss path")

library(fields)
cols <- colorRampPalette(c("white", "blue"))(100)

x11()
image.plot(
  z = fit_real$W[nrow(fit_real$W):1, ],
  col = cols,
  axes = FALSE,
  main = "Estimated feature weights",
  xlab = "State",
  ylab = "Asset"
)

axis(
  1,
  at = seq(0, 1, length.out = nrow(fit_real$W)),
  labels = nrow(fit_real$W):1
)

axis(
  2,
  at = seq(0, 1, length.out = P),
  labels = colnames(Y),
  las = 2
)