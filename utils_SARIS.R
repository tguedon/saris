rm(list=ls())

metropolis = function(alpha, sample_size, burn_in = 100, initial_state = NULL , sigma_prop = NULL, adaptative = TRUE, batch_size = 25){
  
  if (is.null(sigma_prop)){sigma_prop = 1} # Matrice de variances des distributions de proposition
  
  # Initialise la chaîne de Markov
  if (is.null(initial_state)) {
    current_state <- rnorm(1, 0, sigma_prop )
  } else {
    current_state <- initial_state
  }
  
  # Initialise le tableau d'échantillons et le compteur d'acceptation
  samples <- rep(0, sample_size)
  
  acceptance = rep(0, sample_size+1)
  
  sigma_liste = rep(0, sample_size+1)
  sigma_liste[1] = sigma_prop
  
  mean_liste = rep(0,sample_size-burn_in+1)
  
  # Itère sur le nombre souhaité d'échantillons
  for (i in 1:(sample_size )) {
    
    proposed_state = current_state
    proposed_state <- current_state + rnorm(1, 0, sigma_prop)
    
    # calcul de la probabilité d'acceptation
    
    acceptance_prob = min(1,
                          alpha(proposed_state) 
                          /  alpha(current_state)  )
    
    u <- runif(1)
    
    if (is.na(acceptance_prob)){acceptance_prob=0}
    
    if (u < acceptance_prob){
      current_state = proposed_state
      
      # Compte les acceptations
      acceptance[i+1] = acceptance[i] +1
    }
    else{
      
      acceptance[i+1] = acceptance[i] 
      
    }
    
    # Met à jour les variances de la loi de proposition
    if(i>burn_in){if (i%%batch_size==0){
      
      delta = exp(min(0.01, i **(-0.5)))
      
      if(adaptative){
        
        if ( acceptance[i+1]/i > 0.44  ){
          
          sigma_prop = sigma_prop * delta
          
        }
        
        else{
          
          sigma_prop = sigma_prop / delta
          
        }}
      
    }}
    
    
    
    if (i>burn_in){
      
      mean_liste[i-burn_in +1] = ((i-1)/i)*mean_liste[i-burn_in] + current_state / i
      
    }
    
    sigma_liste[i+1] = sigma_prop
    
    # On stock le nouvel état 
    samples[i ] = current_state
    
  }
  
  return(list("samples"=samples, "acceptance" = acceptance[2:length(acceptance)]/(1:sample_size), "sigma_prop" = sigma_liste, "mean_liste" = mean_liste ))
  
  
}

bissection <- function(f, a = 0.0001, b = 10, tol = 1e-6, max_iter = 10000) {
  if (f(a) * f(b) >= 0) {
    return(-10)  }
  
  iter <- 0
  while ((b - a) / 2 > tol && iter < max_iter) {
    midpoint <- (a + b) / 2
    
    
    if (f(midpoint) == 0) {
      return(midpoint) # Solution exacte trouvée
    }
    else if (f(a) * f(midpoint) < 0) {
      b <- midpoint
    } else {
      a <- midpoint
    }
    iter <- iter + 1
  }
  
  return((a + b) / 2) # Approximation de la solution
}

BridgeSARIS = function(sampletot, f0,f1,r0, stepsize,epsilon = 0){
  
  # sampletot simulé selon alpha
  
  K = length(sampletot)
  iteration = K
  k = 0
  
  
  grad = 0
  r = r0
  r_liste = rep(0,K)
  grad_liste = rep(0,K)
  rAV_liste = rep(0,K)
  
  
  r_liste[1] = r0
  grad_liste[1] = grad
  rAV_liste[1] = r0
  
  
  for (k in 1:(K-1)){
    
    f0_x = f0(sampletot[k])
    f1_x = f1(sampletot[k])
    
    den = f0_x+r*f1_x
    if (den == 0){grad=0}else{grad = 1-2*f0_x/den}
    
    r = max(r - stepsize(k) * grad,0)
    
    
    
    r_liste[k+1] = r
    grad_liste[k+1] = grad
    rAV_liste[k+1] = rAV_liste[k] + (1/(k+1))*(r-rAV_liste[k])
    
  } 
  
  return(list("r" = r_liste, "grad" = grad_liste, "rAV" = rAV_liste))
}

BridgeSARISlog = function(sampletot, f0,f1,r0, stepsize,epsilon = 0){
  
  # sampletot simulé selon alpha
  
  K = length(sampletot)
  iteration = K
  k = 0
  
  
  grad = 0
  r = r0
  r_liste = rep(0,K)
  grad_liste = rep(0,K)
  rAV_liste = rep(0,K)
  
  
  r_liste[1] = r0
  grad_liste[1] = grad
  rAV_liste[1] = r0
  
  
  for (k in 1:(K-1)){
    
    f0_x = f0(sampletot[k])
    f1_x = f1(sampletot[k])
    
    if(f0_x>f1_x){
      ratio = f1_x/f0_x
      grad = 1-2/(1+ratio*exp(r))
    }else if(f0_x<f1_x){
      ratio = f0_x/f1_x
      grad = 1-2*ratio/(ratio+exp(r))
    }else{grad=0} 
    
    r = r-stepsize(k)*grad
    
    
    
    r_liste[k+1] = r
    grad_liste[k+1] = grad
    rAV_liste[k+1] = rAV_liste[k] + (1/(k+1))*(r-rAV_liste[k])
    
  } 
  
  return(list("r" = r_liste, "grad" = grad_liste, "rAV" = rAV_liste))
}

BridgeSARISlogRM = function(K,sample0,sample1, f0,f1,r0, stepsize,epsilon = 0){
  
  # sampletot simulé selon alpha
  
  
  
  grad = 0
  r = r0
  r_liste = rep(0,K)
  grad_liste = rep(0,K)
  rAV_liste = rep(0,K)
  
  sample_liste = rep(0,K)
  r_liste[1] = r0
  grad_liste[1] = grad
  rAV_liste[1] = r0
  
  
  for (k in 1:K){
    
    proba = 1/(1+exp(r))
    
    if (runif(1)<proba){x=sample(sample1,1)}
    else{x=sample(sample0,1)}
    
    f0_x = f0(x)
    f1_x = f1(x)
    
    grad = (f0_x-f1_x*exp(r))/(f0_x+f1_x*exp(r))
    
    r = r+stepsize(k)*grad
    
    
    r_liste[k+1] = r
    grad_liste[k+1] = grad
    rAV_liste[k+1] = rAV_liste[k] + (1/(k+1))*(r-rAV_liste[k])
    sample_liste[k] = x
    
  } 
  
  return(list("r" = r_liste, "grad" = grad_liste, "rAV" = rAV_liste,sample_liste))
}

SARIS_RM = function(K, r0,f0,f1,stepsize, alpha,acc = 0.44){
  
  # on resoud E_alpha(xi,r)[h(xi,r)] = 0
  
  iteration = K
  k = 0
  
  grad = 0
  r = r0
  r_liste = rep(0,K)
  grad_liste = rep(0,K)
  rAV_liste = rep(0,K)
  
  
  #r_liste[1] = r0
  #grad_liste[1] = grad
  #rAV_liste[1] = r0
  sigma_prop=1
  
  sample_liste <- rep(0, K)
  
  acceptance = rep(0, K)
  
  sigma_liste = rep(0, K)
  switch_liste = rep(0,K)
  current_state = 0
  
  for (k in 1:K){
    
    
    
    
    # simulate : 
    
    
    
    proposed_state <- current_state + rnorm(1, 0, sigma_prop)
    
    # calcul de la probabilité d'acceptation
    
    acceptance_prob = min(1,
                          alpha(proposed_state, r) 
                          /  alpha(current_state,r)  )
    
    u <- runif(1)
    
    if (u < acceptance_prob){
      current_state = proposed_state
      
      # Compte les acceptations
      acceptance[k] = acceptance[max(k-1,1)] +1
      switch_liste[k] = 1
    }
    else{
      
      acceptance[k] = acceptance[max(k-1,1)] 
      
    }
    
    # Met à jour les variances de la loi de proposition
    
    
    # Met à jour les variances de la loi de proposition
    if (k%%25==0){
      
      delta = exp(min(0.01, k **(-0.5)))
      
      
      
      if ( acceptance[k]/k > acc  ){
        
        sigma_prop = sigma_prop * delta
        
      }
      
      else{
        
        sigma_prop = sigma_prop / delta
        
      }
      
    }
    
    
    sigma_liste[k] = sigma_prop
    sample_liste[k] = current_state
    # nouveau grad 
    grad = (f1(current_state)*r-f0(current_state))/alpha(current_state,r)
    
    
    
    
    # grad and update 
    
    
    r = r - stepsize(k) * grad
    
    
    
    r_liste[k] = r
    grad_liste[k] = grad
    
    
    rAV_liste[k] = rAV_liste[max(k-1,1)] + (1/(k))*(r-rAV_liste[max(k-1,1)])}
  
  
  
  return(list("r" = r_liste, "grad" = grad_liste, "rAV" = rAV_liste, "sigma_prop" = sigma_liste, "sample" = sample_liste, "acceptance" = acceptance/(1:length(acceptance)), "switch_liste"=switch_liste))
  
}

SARIS_RMlog = function(K, r0,f0,f1,stepsize, alpha,acc = 0.44){
  
  # on resoud E_alpha(xi,r)[h(xi,r)] = 0
  
  iteration = K
  k = 0
  
  grad = 0
  r = r0
  r_liste = rep(0,K)
  grad_liste = rep(0,K)
  rAV_liste = rep(0,K)
  
  
  #r_liste[1] = r0
  #grad_liste[1] = grad
  #rAV_liste[1] = r0
  sigma_prop=1
  
  sample_liste <- rep(0, K)
  
  acceptance = rep(0, K)
  
  sigma_liste = rep(0, K)
  switch_liste = rep(0,K)
  current_state = 0
  
  for (k in 1:K){
    
    
    
    
    # simulate : 
    
    
    
    proposed_state <- current_state + rnorm(1, 0, sigma_prop)
    
    # calcul de la probabilité d'acceptation
    
    acceptance_prob = min(1,
                          alpha(proposed_state, exp(r)) 
                          /  alpha(current_state,exp(r))  )
    
    u <- runif(1)
    
    if (u < acceptance_prob){
      current_state = proposed_state
      
      # Compte les acceptations
      acceptance[k] = acceptance[max(k-1,1)] +1
      switch_liste[k] = 1
    }
    else{
      
      acceptance[k] = acceptance[max(k-1,1)] 
      
    }
    
    # Met à jour les variances de la loi de proposition
    
    
    # Met à jour les variances de la loi de proposition
    if (k%%25==0){
      
      delta = exp(min(0.01, k **(-0.5)))
      
      
      
      if ( acceptance[k]/k > acc  ){
        
        sigma_prop = sigma_prop * delta
        
      }
      
      else{
        
        sigma_prop = sigma_prop / delta
        
      }
      
    }
    
    
    sigma_liste[k] = sigma_prop
    sample_liste[k] = current_state
    # nouveau grad 
    grad = (f1(current_state)*exp(r)-f0(current_state))/alpha(current_state,exp(r))
    
    
    
    
    # grad and update 
    
    
    r = r - stepsize(k) * grad
    
    
    
    r_liste[k] = r
    grad_liste[k] = grad
    
    
    rAV_liste[k] = rAV_liste[max(k-1,1)] + (1/(k))*(r-rAV_liste[max(k-1,1)])}
  
  
  
  return(list("r" = r_liste, "grad" = grad_liste, "rAV" = rAV_liste, "sigma_prop" = sigma_liste, "sample" = sample_liste, "acceptance" = acceptance/(1:length(acceptance)), "switch_liste"=switch_liste))
  
}

### joint proc 

# fonction qui calculent les apramètres de la loi à postérieori des X2

var_post = function(beta,gamma2,sigma2){return(gamma2*sigma2/((beta**2)*gamma2+sigma2))}

moy_post = function(y,x1,mu,beta0,beta1,beta2,gamma2,sigma2){
  
  deltay = y-beta0-beta1*x1
  v = var_post(beta=beta2, gamma2 = gamma2, sigma2=sigma2)
  
  return((v * (deltay*beta2/sigma2 + mu/gamma2)))
}

# fonctions qui calculent les gradients des différents paramètres

grad_beta0 = function(yi, beta0, beta, xi, sigma2){return((yi-beta0-t(xi)%*%beta)/sigma2)}
grad_beta = function(yi, beta0, beta, xi, sigma2){return(((yi-beta0-t(xi)%*%beta)/sigma2)*xi)}

grad_mu = function(xi,mu,gamma_vec){return((xi-mu)/gamma_vec)}
grad_gamma = function(xi,mu,gamma_vec){
  
  return( (xi-mu)**2/(2*gamma_vec**2) - 1/(2*gamma_vec)  )
  
}

# fonctions qui calclent la vraisemblance complète

complete_llik  = function(yi, beta0, beta, xi, sigma2, mu ,gamma_vec){return(log(dnorm(yi, mean = beta0+t(xi)%*%beta,sd = sqrt(sigma2))) +log( dnorm(xi[1], mean = mu[1], sd = sqrt(gamma_vec[1]))) +log( dnorm(xi[2], mean = mu[2], sd = sqrt(gamma_vec[2])) ))}

marg_lik = function(y, beta0, beta, x, sigma2, mu ,gamma_vec,n_mis){
  #x avec les na
  nindiv = length(y)
  res = 0
  nindiv = length(y)
  for (i in 1:n_mis){
    res = res + log( dnorm(y[i], mean = beta0+beta[1]*x[i,1] +mu[2]*beta[2], sd = sqrt(beta[2]**2*gamma_vec[2]+sigma2))  ) + log( dnorm(x[i,1], mean = mu[1], sd = sqrt(gamma_vec[1])))
    
  }
  
  for(i in (n_mis+1):nindiv){    
    res = res + log(dnorm(y[i], mean = beta0+t(x[i,])%*%beta,sd = sqrt(sigma2))) +log( dnorm(x[i,1], mean = mu[1], sd = sqrt(gamma_vec[1]))) +log( dnorm(x[i,2], mean = mu[2], sd = sqrt(gamma_vec[2])) )
  }
  
  res
  
}

