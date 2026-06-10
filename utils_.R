feat_weight_jump <- function(Y,
                             zeta,
                             lambda,
                             K,
                             tol     = FALSE,
                             n_init  = 10,
                             n_iter  = 10,
                             verbose = FALSE,
                             mif     = NULL,
                             truth   = NULL,
                             ncores  = NULL,
                             tukey   = TRUE) {
  
  # Fits the Feature-Weighted Jump Model (FWJM) for time-dependent clustering.
  # The method jointly estimates a latent state sequence, state-specific medoids,
  # and state-specific feature weights under a temporal smoothness penalty.
  # Continuous variables can be handled robustly through Tukey's biweight
  # transformation.
  #
  # Arguments
  #   Y        : data frame or matrix of observations (T x P).
  #   zeta     : positive regularization parameter controlling the concentration
  #              of feature weights. Smaller values produce sparser weights.
  #   lambda   : non-negative temporal penalty controlling state persistence.
  #   K        : number of latent states.
  #   tol      : optional convergence tolerance on the objective function.
  #   n_init   : number of random initializations.
  #   n_iter   : maximum number of iterations per initialization.
  #   verbose  : logical; if TRUE, prints progress information.
  #   mif      : optional variable used to reorder state labels.
  #   truth    : optional true state sequence used for ARI computation.
  #   ncores   : number of cores for parallel execution (only available on Mac and Linux)
  #   tukey    : logical; if TRUE, uses Tukey's biweight transformation for
  #              continuous-feature dissimilarities.
  #
  # Value
  #   A list containing:
  #     W        : estimated state-specific feature weights (K x P).
  #     s        : estimated latent state sequence.
  #     medoids  : estimated state medoids.
  #     loss     : final objective value.
  #     loss_vec : objective function trajectory.
  #     elapsed  : computation time.

  P  <- ncol(Y)
  TT <- nrow(Y)
  
  start <- Sys.time()
  
  feat_type <- sapply(Y, class)
  cont_feat <- which(feat_type %in% c("numeric","integer"))
  P_cont    <- length(cont_feat)
  cat_ord_feat <- which(feat_type %in% c("factor","ordinal"))
  P_cat_ord <- length(cat_ord_feat)
  
  Gamma <- lambda * (1 - diag(K))
  
  run_one <- function(init_id) {
    
    W        <- matrix(1 / P, nrow = K, ncol = P)
    W_old    <- W
    loss_old <- Inf
    ARI      <- NA
    loss <- Inf
    eps_loss <- Inf
    medoids <- rep(NA_integer_, K)
    
    loss_vec <- data.frame(iter = 0, 
                           loss = Inf,
                           ARI = NA)
    
    Y2   <- data.frame(lapply(Y, function(x) if (is.character(x)) factor(x) else x))
    Ymat <- data.matrix(Y2)
    
    cls <- vapply(Y, function(x) class(x)[1], character(1))
    feat_type_int <- ifelse(cls %in% c("numeric","integer"), 0, 1)
    
    s <- initialize_states(Ymat, K, feat_type = feat_type_int,
                           reps = 50, scale = "m")
    
    # precompute dttp 
    dttp <- array(NA_real_, dim = c(TT, TT, P))
    
    if (P_cont > 0) {
      if(tukey){
        for (p_idx in seq_len(P_cont)) {
          p <- cont_feat[p_idx]
          x <- Y[, p]
          
          temp <- abs(outer(x, x, "-"))
          sc_mad <- mad(x, constant = 1, na.rm = TRUE)
          if (sc_mad <= 0) sc_mad <- 1
          temp <- temp / sc_mad
          
          #if (tukey) {
          temp_vec <- tukey_biw_vec_cpp(as.numeric(temp), 4.685)
          dim(temp_vec) <- dim(temp)
          temp <- temp_vec
          #}
          
          temp <- safe_scale_slice_median_cpp(temp, scale = "m")
          dttp[, , p] <- temp
        }
      }
      else{
        Y_cont_mat <- data.matrix(Y[, cont_feat, drop = FALSE])
        gows <- gower_dist_array(Y_cont_mat, Y_cont_mat, scale = "m")
        dttp[, , cont_feat] <- gows
      }
    }
    
    if (P_cat_ord > 0) {
      Y_cat_mat <- data.matrix(Y[, cat_ord_feat, drop = FALSE])
      gows <- gower_dist_array(Y_cat_mat, Y_cat_mat, scale = "m")
      dttp[, , cat_ord_feat] <- gows
    }
    
    # 
    # MAIN LOOP
    # 
    for (iter in seq_len(n_iter)) {
      
      #  D_max 
      D_left  <- matrix(0, TT, TT)
      D_right <- matrix(0, TT, TT)
      
      for (p in seq_len(P)) {
        Dp <- dttp[, , p]
        D_left  <- D_left  + Dp * W[s, p]
        D_right <- D_right + t(t(Dp) * W[s, p])
      }
      
      D_mean <- pmax(D_left, D_right)
      
      # Find medoids 
      medoids <- pam_fixed_cpp(D_mean, s, K)
      
      #  loss_by_state
      loss_by_state <- matrix(0, TT, K)
      for (k in seq_len(K)) {
        for (p in seq_len(P)) {
          loss_by_state[, k] <- loss_by_state[, k] +
            W[k, p] * dttp[, medoids[k], p]
        }
      }
      
      #  E-step 
      Estep <- E_step(loss_by_state, Gamma)
      V     <- Estep$V
      s     <- Estep$s
      
      if (length(unique(s)) == 1) {
        break
      }
      
      tab_s <- tabulate(s, nbins = K)
      empty <- which(tab_s == 0)
      
      # fix empty clusters by moving one point from the closest cluster (in terms of loss increase)
      if (length(empty) > 0) {
        for (k_empty in empty) {
          
          donor_clusters <- which(tab_s > 1)
          
          if (length(donor_clusters) == 0) {
            break
          }
          
          cand_idx <- which(s %in% donor_clusters)
          
          current_loss <- loss_by_state[cbind(cand_idx, s[cand_idx])]
          new_loss     <- loss_by_state[cand_idx, k_empty]
          
          cost_increase <- new_loss - current_loss
          
          candidate <- cand_idx[which.min(cost_increase)]
          
          old_k <- s[candidate]
          s[candidate] <- k_empty
          
          tab_s[old_k] <- tab_s[old_k] - 1
          tab_s[k_empty] <- tab_s[k_empty] + 1
        }
      }
      
      idx_by_k <- lapply(seq_len(K), function(k) which(s == k))
      nk <- tabulate(s, nbins = K)
      
      #  update W 
      Spk <- matrix(0, K, P)
      for (k in seq_len(K)) {
        idx <- idx_by_k[[k]]
        if (length(idx) > 0) {
          for (p in seq_len(P)) {
            Spk[k,p] <- sum(dttp[idx, medoids[k], p])
          }
        }
      }
      
      wcd <- exp(-Spk / zeta)
      wcd[!is.finite(wcd)] <- 0
      
      rs <- rowSums(wcd)
      bad_rs <- !is.finite(rs) | rs <= 0
      
      if (any(bad_rs)) {
        wcd[bad_rs, ] <- 1 / P
        rs[bad_rs] <- 1
      }
      
      W <- wcd / rs
      
      W_old <- W
      
      # loss 
      loss_by_state_final <- matrix(0, TT, K)
      
      for (k in seq_len(K)) {
        for (p in seq_len(P)) {
          loss_by_state_final[, k] <- loss_by_state_final[, k] +
            W[k, p] * dttp[, medoids[k], p]
        }
      }
      
      loss <- sum(loss_by_state_final[cbind(seq_len(TT), s)]) +
        lambda * sum(s[-1] != s[-TT]) +
        zeta * sum(W * log(pmax(W, 1e-16)))
      
      # check for non-finite loss
      if (!is.finite(loss)) {
        loss <- Inf
        break
      }
      
      eps_loss <- abs(loss - loss_old)
      if (!is.finite(eps_loss)) eps_loss <- Inf
      loss_old <- loss
      
      if (!is.null(truth)) {
        ARI <- mclust::adjustedRandIndex(truth, s)
      }
      
      loss_vec <- rbind(loss_vec,
                        c(iter, loss, ARI))
      
      #  PRINT
      if (verbose) {
        cat(sprintf(
          "init %2d | iter %3d/%3d | loss=%.4e | dloss=%.3e | ARI=%.4f\n",
          init_id, iter, n_iter, loss, eps_loss, ARI
        ))
      }
      
      if (!is.null(tol) && is.numeric(tol) && tol > 0 &&
          is.finite(eps_loss) && eps_loss < tol) {
        break
      }
      
    }
    
    list(W=W, s=s, medoids=medoids,
         loss=loss, loss_vec=loss_vec[-1,], ARI=ARI)
  }
  
  #  multi init 
  if (!is.null(ncores)) {
    res_list <- mclapply(seq_len(n_init), run_one, mc.cores = ncores)
  } else {
    res_list <- lapply(seq_len(n_init), run_one)
  }
  
  losses   <- vapply(res_list, `[[`, numeric(1), "loss")
  best_run <- res_list[[which.min(losses)]]
  best_s     <- best_run$s
  best_loss  <- best_run$loss
  loss_vec   <- best_run$loss_vec
  best_W     <- best_run$W
  best_medoids <- Y[best_run$medoids, ]
  #
  if (!is.null(mif)) {
    new_best_s <- order_states_condMed(Y[, mif], best_s, 
                                       decreasing = F)
    tab <- table(
      factor(best_s,    levels = 1:K),
      factor(new_best_s, levels = 1:K)
    )
    
    perm   <- apply(tab, 2, which.max)  # per ogni nuova etichetta, la vecchia
    best_W <- best_W[perm, , drop = FALSE]
    
    best_medoids <- best_medoids[perm, , drop = FALSE]
    
  } else {
    new_best_s <- best_s
  }
  
  end <- Sys.time()
  
  return(list(
    W       = best_W,
    s       = new_best_s,
    medoids = best_medoids,
    loss    = best_loss,
    loss_vec = best_run$loss_vec,
    elapsed = end - start
  ))
}

order_states_condMed <- function(y, s, decreasing = FALSE) {
  
  # Relabels latent states according to the ordering of their conditional
  # medians for a given variable.
  #
  # Arguments
  #   y          : numeric vector used for ordering states.
  #   s          : state sequence.
  #   decreasing : logical; if TRUE, orders states from largest to smallest
  #                conditional median.
  #
  # Value
  #   Integer vector containing the relabeled state sequence.
  
  condMed <- sort(tapply(y, s, median, na.rm = TRUE), decreasing = decreasing)
  
  states_temp <- match(s, names(condMed))
  
  return(states_temp)
}

sim_data_stud_t=function(seed=123,
                         TT,
                         P,
                         Pcat=NULL,
                         Ktrue=3,
                         mu=1.5,
                         rho=0,
                         nu=4,
                         phi=.8,
                         pers=.95){
  
  # Simulates a hidden Markov model with multivariate Student-t emissions.
  #
  # Arguments
  #   seed   : random seed.
  #   TT     : number of observations.
  #   P      : number of continuous variables.
  #   Pcat   : number of variables converted to categorical variables.
  #   Ktrue  : number of latent states.
  #   mu     : separation parameter controlling state means.
  #   rho    : pairwise correlation among variables.
  #   nu     : degrees of freedom of the Student-t distribution.
  #   phi    : parameter used when generating categorical variables.
  #   pers   : self-transition probability of the latent Markov chain.
  #
  # Value
  #   A list containing:
  #     SimData : simulated dataset.
  #     mchain  : true latent state sequence.
  #     TT      : number of observations.
  #     P       : number of variables.
  #     K       : number of states.
  #     Ktrue   : true number of states.
  #     pers    : persistence parameter.
  #     seed    : random seed used.
  
  MU=seq(-mu, mu, length.out=Ktrue)
  
  # Markov chain simulation
  x <- numeric(TT)
  Q <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(Q)=rep(pers,Ktrue)
  init <- rep(1/Ktrue,Ktrue)
  set.seed(seed)
  x[1] <- sample(1:Ktrue, 1, prob = init)
  for(i in 2:TT){
    x[i] <- sample(1:Ktrue, 1, prob = Q[x[i - 1], ])
  }
  
  # Continuous variables simulation
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, TT, P * Ktrue)
  SimData = matrix(0, TT, P)
  
  set.seed(seed)
  for(k in 1:Ktrue){
    u = mvtnorm::rmvt(TT, sigma = (nu-2)*Sigma/nu, df = nu, delta = rep(MU[k],P))
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:TT) {
    k = x[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
  }
  
  SimData=data.frame(SimData)
  
  if(!is.null(Pcat)){
    for (j in 1:Pcat) {
      SimData[, j] <- get_cat_t(SimData[, j], x, MU, phi=phi, df = nu)
      SimData[, j]=factor(SimData[, j],levels=1:Ktrue)
    }  
  }
  
  
  
  return(list(
    SimData=SimData,
    mchain=x,
    TT=TT,
    P=P,
    K=Ktrue,
    Ktrue=Ktrue,
    pers=pers, 
    seed=seed))
  
}

simulate_sparse_hmm <- function(Y,
                                rel_,
                                true_stat,
                                perc_out   = 0.02,
                                seed       = NULL,
                                infl_fact=1,
                                Out_bound   = 30) {
  
  
  # Transforms a simulated HMM dataset into a sparse state-dependent structure by
  # assigning subsets of informative variables to each latent state and
  # optionally introducing outliers.
  #
  # Arguments
  #   Y          : simulated dataset.
  #   rel_       : list specifying which variables are relevant for each state.
  #   true_stat  : true latent state sequence.
  #   perc_out   : proportion of observations replaced by outliers.
  #   seed       : random seed.
  #   infl_fact  : scaling factor for injected noise.
  #   Out_bound  : magnitude of outlier contamination.
  #
  # Value
  #   A list containing:
  #     Y           : contaminated dataset.
  #     truth       : true state sequence with outliers coded as 0.
  #     out_indices : indices of injected outliers.
  #     W_truth     : true binary state-feature relevance matrix.
  
  if(!is.null(seed)) set.seed(seed)
  Y <- as.matrix(Y)
  TT  <- nrow(Y)
  P  <- ncol(Y)
  K=length(rel_)
  # 1) invert the rel_ list: for each feature p, which states mention p?
  inv_rel <- invert_rel(rel_, P)
  
  # 2) Irrelevant features = those never mentioned in rel_
  irrelevant <- which(vapply(inv_rel, length, integer(1)) == 0)
  if(length(irrelevant) > 0) {
    # permute their rows globally
    # Y[, irrelevant] <- Y[sample(TT), irrelevant]
    
    # Uniform noise
    Y[, irrelevant] <- infl_fact*matrix(
      runif(length(irrelevant) * TT,
            min = min(Y),
            max = max(Y)),
      nrow = TT, ncol = length(irrelevant)
    )
    
    # all_obs=as.vector(Y)
    # Y[, irrelevant]=matrix(
    #   sample(all_obs,length(irrelevant) * TT,replace=T),
    #   nrow = TT, ncol = length(irrelevant)
    # )
    
  }
  
  # 3) Relevant features = those that appear in at least one state
  relevant <- which(vapply(inv_rel, length, integer(1)) > 0)
  for(p in relevant) {
    relevant_states <- inv_rel[[p]]
    # indices of rows belonging to any of those states
    idx_in_state <- which(true_stat %in% relevant_states)
    # rows not in those states:
    idx_out_state <- setdiff(seq_len(TT), idx_in_state)
    if(length(idx_out_state) > 1) {
      # permute only those rows of column p
      # Y[idx_out_state, p] <- sample(Y[idx_out_state, p])
      
      # Uniform noise
      Y[idx_out_state, p] =
        infl_fact*runif(length(idx_out_state),
                        min = min(Y),
                        max = max(Y))
      # Y[idx_out_state, p]=sample(all_obs,
      #                            length(idx_out_state),
      #                            replace=T)
      
      
    }
  }
  
  # 4) Introduce outliers
  if (perc_out > 0) {
    
    n_out <- ceiling(TT * perc_out)
    idx_out <- sample(seq_len(TT), n_out)
    
    mins <- apply(Y, 2, min, na.rm = TRUE)
    maxs <- apply(Y, 2, max, na.rm = TRUE)
    
    Y_out <- matrix(NA_real_, nrow = n_out, ncol = P)
    
    for (j in seq_len(P)) {
      side <- sample(c("low", "high"), n_out, replace = TRUE)
      
      Y_out[side == "low",  j] <- runif(sum(side == "low"),
                                        mins[j] - Out_bound,
                                        mins[j])
      
      Y_out[side == "high", j] <- runif(sum(side == "high"),
                                        maxs[j],
                                        maxs[j] + Out_bound)
    }
    
    Y[idx_out, ] <- Y_out
    # 5) update truth: set outlier rows to 0
    new_truth <- true_stat
    new_truth[idx_out] <- 0L
  }
  else{
    new_truth <- true_stat
    idx_out=NULL
  }
  W_truth <- matrix(FALSE, nrow = K, ncol = P)
  
  for (k in seq_len(K)) {
    W_truth[k, rel_[[k]]] <- TRUE
  }
  
  list(
    Y          = data.frame(Y),
    truth      = new_truth,
    out_indices = idx_out,
    W_truth = W_truth
  )
}

invert_rel <- function(rel_, P) {
  
  # Constructs the inverse representation of a state-feature relevance
  # specification.
  #
  # Arguments
  #   rel_ : list where rel_[[k]] contains the variables relevant to state k.
  #   P    : total number of variables.
  #
  # Value
  #   A list of length P where the p-th element contains all states for which
  #   variable p is relevant.

  inv <- vector("list", P)
  for(i in seq_len(P)) inv[[i]] <- integer(0)
  
  # for each group k, append k to every member i in rel_[[k]]
  for(k in seq_along(rel_)) {
    members <- rel_[[k]]
    for(i in members) {
      inv[[i]] <- c(inv[[i]], k)
    }
  }
  
  inv
}