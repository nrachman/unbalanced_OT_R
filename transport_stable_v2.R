#transport_stablev2 comes from https://github.com/broadinstitute/wot/blob/master/wot/ot/ot_model.py
# If I remember correctly, this is  algorithm 2 of Chizat 2016 - Scaling Algorithms for Unbalanced Transport Problems
# https://arxiv.org/abs/1607.05816

# parameters come from
# https://github.com/broadinstitute/wot/blob/master/wot/ot/ot_model.py
# self.ot_config = {'local_pca': 30, 'growth_iters': 1, 'epsilon': 0.05, 'lambda1': 1, 'lambda2': 50,
#   'epsilon0': 1, 'tau': 10000, 'scaling_iter': 3000, 'inner_iter_max': 50, 'tolerance': 1e-8,
#   'max_iter': 1e7, 'batch_size': 5, 'extra_iter': 1000}

# epsilon <- 0.05
# lambda1 <- 1
# lambda2 <- 50
# epsilon0 <- 1 
# tau <- 10000
# scaling_iter<- 3000
# inner_iter_max <- 50
# tolerance <- 1e-8
# max_iter <-  1e7
# batch_size <- 5
# extra_iter <- 1000

transport_stablev2 <- function(C, lambda1 = 1, lambda2 = 50, epsilon = .05, scaling_iter = 3000, G = NULL, tau = 10000, epsilon0 = 1,
                               extra_iter = 1000, inner_iter_max = , ...) {
  
  # """
  #   Compute the optimal transport with stabilized numerics.
  #   Args:
  # 
  #       C: cost matrix to transport cell i to cell j
  #       lambda1: regularization parameter for marginal constraint for p.
  #       lambda2: regularization parameter for marginal constraint for q.
  #       epsilon: entropy parameter
  #       scaling_iter: number of scaling iterations
  #       G: growth value for input cells
  #   """
  

  # This function computes an entropy-regularized optimal transport plan between source and target distributions
  # with stabilization for numerical robustness.
  
  if(is.null(G)){
    # If G (the source distribution) is not provided, initialize it as uniform over rows of C
    G <- rep(1, nrow(C))
  }
  
  # Helper function that computes a decreasing epsilon schedule
  get_reg <- function(n) {
    (epsilon0 - epsilon_final) * exp(-n) + epsilon_final
  }
  
  warm_start <- !is.null(tau)  # Enable warm start if tau is given (used for stabilization check)
  epsilon_final <- epsilon     # Final target epsilon
  
  # Initial entropy regularization value
  epsilon_i <- if (warm_start) epsilon0 else epsilon
  
  # Uniform weights for marginal distributions
  dx <- rep(1 / nrow(C), nrow(C))  # uniform mass for rows
  dy <- rep(1 / ncol(C), ncol(C))  # uniform mass for columns
  
  # Source distribution p and target distribution q
  p <- G
  q <- rep(mean(G), ncol(C))  # q is initialized to the mean of G repeated to match target size
  
  # Dual variables (log-scaling vectors)
  u <- rep(0, length(p))
  v <- rep(0, length(q))
  
  # Initialize scaling factor for target side
  b <- rep(1, length(q))
  
  # Kernel matrix K initialized with entropy-regularized cost
  K <- exp(-C / epsilon_i)
  
  # Compute regularization scaling factors
  alpha1 <- lambda1 / (lambda1 + epsilon_i)
  alpha2 <- lambda2 / (lambda2 + epsilon_i)
  
  epsilon_index <- 0  # counter for how many times epsilon has been updated
  iterations_since_epsilon_adjusted <- 0  # counter for stabilization
  
  # Main Sinkhorn scaling iterations
  for (i in 1:scaling_iter) {
    # Update scaling vector a (source side)
    a <- (p / (K %*% (b * dy)))^alpha1 * exp(-u / (lambda1 + epsilon_i))
    
    # Update scaling vector b (target side)
    b <- (q / (t(K) %*% (a * dx)))^alpha2 * exp(-v / (lambda2 + epsilon_i))
    
    # Count iteration since last stabilization
    iterations_since_epsilon_adjusted <- iterations_since_epsilon_adjusted + 1
    
    # If scaling factors become too large or small, stabilize
    if (max(abs(a)) > tau || max(abs(b)) > tau) {
      # Absorb scalings into dual potentials
      u <- u + epsilon_i * log(a)
      v <- v + epsilon_i * log(b)
      
      # Recompute kernel K with updated duals
      K <- exp((matrix(u, nrow = length(u), ncol = length(v)) - C +
                  matrix(v, nrow = length(u), ncol = length(v), byrow = TRUE)) / epsilon_i)
      
      # Reset scaling vectors
      a <- rep(1, length(p))
      b <- rep(1, length(q))
    }
    
    # If warm_start is enabled and enough iterations passed, reduce epsilon
    if (warm_start && iterations_since_epsilon_adjusted == inner_iter_max) {
      epsilon_index <- epsilon_index + 1
      iterations_since_epsilon_adjusted <- 0
      
      # Absorb scalings into dual potentials
      u <- u + epsilon_i * log(a)
      v <- v + epsilon_i * log(b)
      
      # Update epsilon based on schedule
      epsilon_i <- get_reg(epsilon_index)
      
      # Update alpha coefficients for new epsilon
      alpha1 <- lambda1 / (lambda1 + epsilon_i)
      alpha2 <- lambda2 / (lambda2 + epsilon_i)
      
      # Recompute kernel with new epsilon
      K <- exp((matrix(u, nrow = length(u), ncol = length(v)) - C +
                  matrix(v, nrow = length(u), ncol = length(v), byrow = TRUE)) / epsilon_i)
      
      # Reset scaling vectors
      a <- rep(1, length(p))
      b <- rep(1, length(q))
    }
  }
  
  # Extra refinement iterations after stabilization
  for (i in 1:extra_iter) {
    a <- (p / (K %*% (b * dy)))^alpha1 * exp(-u / (lambda1 + epsilon_i))
    b <- (q / (t(K) %*% (a * dx)))^alpha2 * exp(-v / (lambda2 + epsilon_i))
  }
  
  # Compute final transport plan R = diag(a) K diag(b)
  R <- sweep(K, 1, a, "*")
  R <- sweep(R, 2, b, "*")
  
  # Normalize result by number of columns (can be omitted depending on context)
  return(R / ncol(C))
}
