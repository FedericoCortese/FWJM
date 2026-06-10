# Feature-Weighted Jump Model (FWJM)

This repository contains an implementation of the **Feature-Weighted Jump Model (FWJM)** for time-dependent clustering with state-specific feature relevance.

The method jointly estimates:

* a latent state sequence,
* state-specific medoids,
* state-specific feature weights,

while encouraging temporal persistence through a transition penalty.

Robustness to outliers is achieved through the use of **Tukey's biweight loss**.

## Main Function

```r
fit <- feat_weight_jump(
  Y       = data,
  zeta    = 10,
  lambda  = 0.5,
  K       = 3,
  n_init  = 10,
  n_iter  = 50,
  tukey   = TRUE
)
```

The function returns:

```r
fit$W         # feature weights
fit$s         # estimated state sequence
fit$medoids   # state medoids
fit$loss      # final objective value
fit$loss_vec  # optimization path
```

---

## Toy Example

A complete toy example is provided using simulated data generated from a hidden Markov model with Student-(t) emissions.

The example demonstrates:

1. Simulation of latent states.
2. Construction of state-specific relevant features.
3. Injection of outliers.
4. FWJM estimation.
5. Recovery of state sequences.
6. Recovery of feature relevance patterns.

---

## Real Data Example

A real-data example is included using financial assets downloaded from Yahoo Finance.

The dataset contains several asset classes:

| Ticker | Asset Class          |
| ------ | -------------------- |
| SPY    | US Equity            |
| EFA    | International Equity |
| EEM    | Emerging Equity      |
| TLT    | Government Bonds     |
| LQD    | Corporate Bonds      |
| HYG    | High Yield Bonds     |
| VNQ    | Real Estate          |
| GLD    | Gold                 |
| DBC    | Commodities          |
| UUP    | US Dollar            |

Hyperparameters ((K,\lambda,\zeta)) are selected using a grid search maximizing the median silhouette coefficient.

---

## Hyperparameters

### Number of States

```r
K
```

Number of latent states.

### Temporal Penalty

```r
lambda
```

Controls persistence of the latent state sequence.

* Small values allow frequent regime changes.
* Large values encourage longer-lasting states.

### Feature-Weight Concentration

```r
zeta
```

Controls the distribution of feature weights.

* Small values yield sparse state-specific weights.
* Large values yield more uniform weights.

---

## Dependencies

Required R packages:

```r
Rcpp
parallel
cluster
mclust
mvtnorm
quantmod
xts
zoo
```

---

## References

TO DO
