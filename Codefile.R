
#-------------------------------------------------------------------------------
library(readr)
library(tidyverse)
#----Perch approach-------------------------------------------------------------
#line 65 - 132 are not relevant when fake data is implemented
#-----Chunk 1 - generating fake data--------------------------------------------
#generate fake data
new_validation_example = function(filename,
                                  timestamp_offset,
                                  score, 
                                  is_pos, 
                                  bin, 
                                  bin_weight) {
  list(
    filename = filename,
    timestamp_offset = timestamp_offset,
    score = score,
    is_pos = is_pos,
    bin = bin,
    bin_weight = bin_weight
  )
}

#Generating the fake data
generate_fake_validation_data = function(n = 100, bins = 5) {
  scores = rexp(n,rate =6)
  scores = scores/max(scores) # to make sure it <= 1
  quantile_bounds = seq(0,
                        1,
                        length.out = bins + 1)
  bin_indices = findInterval(scores,
                             quantile_bounds,
                             rightmost.closed = TRUE)           #returns the bin number for each score
  bin_weights = table(bin_indices) / n                          #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)]) #gives the bin proportion associated with each score
  
  examples = vector("list",
                    n)
  
  for (i in seq_len(n)) {
    ex = new_validation_example(
      filename = paste0("file_", sample(1:10, 1), ".wav"),
      timestamp_offset = runif(1, 0, 60),
      score = scores[i],
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.20, 0.7, 0.1)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] = ex
  }
  return(examples)
}

generate_fake_validation_data2 = function(n = 100, bins = 5) {
  scores = rexp(n,rate = 6)
  scores = scores/max(scores) # to makes ure it is <= 1
  quantile_bounds = seq(0,
                        1,
                        length.out = bins + 1)
  bin_indices = findInterval(scores,
                             quantile_bounds,
                             rightmost.closed = TRUE)           #returns the bin number for each score
  bin_weights = table(bin_indices) / n                          #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)]) #gives the bin proportion associated with each score
  
  examples = vector("list",
                    n)
  
  for (i in seq_len(n)) {
    ex <- new_validation_example(
      filename = paste0("file_", sample(1:10, 1), ".wav"),
      timestamp_offset = runif(1, 0, 60),
      score = scores[i],
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.20, 0.35, 0.45)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] <- ex
  }
  return(examples)
}

#debug(generate_fake_validation_data)
fdata = generate_fake_validation_data(n = 1000,
                                      bins= 5)
fdata2 = generate_fake_validation_data2(n = 1000,
                                       bins= 5)

study_level_data = rbind(fdata,fdata2)


# -------Chunk 2 - Call density estimation-------------------------------------
#binning and weighting are calculalted within the function
Study_level_call_density = function(examples,num_beta_samples=1000,beta_prior=0.1){
  
  #loading in of examples data and converting it to usable format
  examples = unlist(examples)
  df = as.data.frame(do.call(rbind, study_level_data), stringsAsFactors = FALSE)  #This will be the one that i check my
  logit.scores = sapply(unlist(df$score), function(b) log(b)-log(1-b))            #generating logit scores
  df$score = logit.scores                                                         #turning scores to logit scores
  examples = as.data.frame(do.call(rbind, study_level_data), 
                           stringsAsFactors = FALSE)                              #Conversion to data-frame
  examples = rename(examples,binNumber=bin, outcome = is_pos)                     #renaming of columns

 
  #Binning and bin weighting 
  #Positive/negative bin counts
  bin_number  = as.matrix(distinct(examples,binNumber)) 
  bin_number = as.matrix(sort(as.numeric(as.vector(bin_number[,1]))))
  bin_weights = matrix(NA, nrow = nrow(bin_number),ncol=1)
  bin_positive = matrix(NA,nrow = nrow(bin_number), ncol = 2)
  bin_negative = matrix(NA,nrow = nrow(bin_number), ncol = 2)
  
  #Bin weights
  for(i in 1:nrow(bin_number)){
    bin_weights[i,1]=nrow(filter(examples,binNumber==i))/nrow(examples)
  }
  
  #Putting in the counts of each bin (positive/negative)
  for(i in 1:nrow(bin_number)){
    examples2 = filter(examples,binNumber==i)
    bin_positive[i,2] = as.numeric(count(filter(examples2,outcome==1)))
    bin_negative[i,2] = as.numeric(count(filter(examples2,outcome==-1)))
  }
  bin_positive[,1]=bin_number ; bin_negative[,1] = bin_number
  

#Density  
  # Create beta distributions for each bin.
  shape.parameters = matrix(NA,ncol=2,nrow=nrow(bin_number))
  
  for (b in bin_number){
    shape.parameters[b,] = c(
      bin_positive[b,2] + beta_prior, bin_negative[b,2] + beta_prior
    )
  }
  
  
  #density estimation
  shape.parameters
  density_matrix = matrix(NA,nrow = nrow(bin_number),ncol = 1)
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = rbeta(1,shape.parameters[i,1]+ beta_prior,
                             shape.parameters[i,2]+  beta_prior)*bin_weights[i,1]
  }
  
  density = sum(density_matrix)
    
#Bootstrapping: To get distribution of density estimate
  NoBetaSamples = 10000
  q_betas = numeric(NoBetaSamples)
  
  for (i in seq_len(NoBetaSamples)) {
    # draw a bin-level theta_b from the same Beta as the estimator
    betas = rbeta(
      n = nrow(shape.parameters),
      shape1 = shape.parameters[,1] + beta_prior,  # match point estimator
      shape2 = shape.parameters[,2] + beta_prior
    )
    # weighted sum across bins (MATCH the estimator)
    q_betas[i] <- sum(betas * bin_weights[,1])
  }
  
  
  
  return(list(Density_Estimate = density,
              Samples = q_betas, 
              beta_parameters = shape.parameters,
              Study_level_weights = bin_weights,
              densities = density_matrix,
              negative_counts = bin_negative,
              positive_counts = bin_number))
}

#Implementation of call density  
out = Study_level_call_density(study_level_data,1000,0.1)  
out$Density_Estimate
out$Samples # still wrapping my heat around what the samples are supposed to be.
hist(out$Samples)

#95% confidence interval from bootstrap distribution, Variance and standard deviation 
confidence_interval = quantile(out$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out$Samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)


cat("95% Confidence interval: [",round(confidence_interval,4),"]",
    "\nvariance: ",round(Variance_bootstrap_samples,4),
    "\nstandard deviation: ",round(Sd_bootstrap_samples,4))


#---Distribution shifts---------------------------------------------------------
#Strategy 1
#every thing is the same except that now in this case, we are looking to substitute into formulas the matter of having
#sum(P(+|b)*P_s(b))
Strat1_call_density <- function(examples, num_beta_samples=1000, beta_prior=0.1){

  #loading in of examples data and converting it to usable format
  examples = unlist(examples)
  df = as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE)  #This will be the one that i check my
  logit.scores = sapply(unlist(df$score), function(b) log(b)-log(1-b))            #generating logit scores
  df$score = logit.scores                                                         #turning scores to logit scores
  examples = as.data.frame(do.call(rbind, fdata), 
                           stringsAsFactors = FALSE)                              #Conversion to data-frame
  examples = rename(examples,binNumber=bin, outcome = is_pos)                     #renaming of columns
  
  
  #Positive/negative bin counts
  bin_number  = as.matrix(distinct(examples,binNumber)) 
  bin_number = as.matrix(sort(as.numeric(as.vector(bin_number[,1]))))
  bin_weights = matrix(NA, nrow = nrow(bin_number),ncol=1)
  
  
  #Bin weights
  for(i in 1:nrow(bin_number)){
    bin_weights[i,1]=nrow(filter(examples,binNumber==i))/nrow(examples)
  }


#Density  
  #makes use of study level data
  #shape.parameters
  density_matrix = matrix(NA,nrow=nrow(bin_number), ncol =1)
 
  
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix[i,1] =  rbeta(1,out$beta_parameters[i,1]+beta_prior,
                              out$beta_parameters[i,2]+beta_prior)*bin_weights[i,1]
  }
  density = sum(density_matrix)

  
  #Bootstrapping: To get distribution of density estimate
  NoBetaSamples = 10000
  q_betas = numeric(NoBetaSamples)
  
  for (i in seq_len(NoBetaSamples)) {
    # draw a bin-level theta_b from the same Beta as the estimator
    betas = rbeta(
      n = nrow(out$beta_parameters),
      shape1 = out$beta_parameters[,1] + beta_prior,  # match point estimator
      shape2 = out$beta_parameters[,2] + beta_prior
    )
    # weighted sum across bins (MATCH the estimator)
    q_betas[i] <- sum(betas * bin_weights[,1])
  }
  
  
  return(list(Density_Estimate = density,
              Samples = q_betas,
              bin_densities = density_matrix))
}

#Implementation of strategy 1 call density  

out1 = Strat1_call_density(fdata,1000,0.1)
out1$Density_Estimate
hist(out1$Samples)

#95% confidence interval from bootstrap distribution, Variance and standard deviation 
confidence_interval = quantile(out1$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out1$Samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)


cat("95% Confidence interval: [",round(confidence_interval,4),"]",
    "\nvariance: ",round(Variance_bootstrap_samples,4),
    "\nstandard deviation: ",round(Sd_bootstrap_samples,4))

#---------Strategy 2------------------------------------------------------------
Strat2_call_density = function(site_examples, study_fit = out, eps = 1e-12,beta_prior= 0.1) {
  #Site level:  empirical measurements
  #site_examples = fdata
  site_df = as.data.frame(do.call(rbind, site_examples),
                          stringsAsFactors = FALSE) |> 
    dplyr::rename(binNumber = bin)
  logit.scores = sapply(unlist(site_df$score), function(b) log(b)-log(1-b))            #generating logit scores
  site_df$score = logit.scores 

  B = nrow(study_fit$beta_parameters)
  site_bin_weight = sapply(1:B, function(b) mean(as.numeric(site_df$binNumber) == b)) #Site bin weights
  site_bin_weight[site_bin_weight < eps] = eps # to avoid 0s when calculating things like log
  site_bin_weight = site_bin_weight / sum(site_bin_weight)#proportioning to it to make sure it sums to 1 after some may be set to 0
  
  site_bin_positive = matrix(NA,ncol = 2, nrow = B)
  site_bin_negative = matrix(NA,ncol = 2, nrow = B)
  
  for(i in 1:B){
    dataframe = filter(site_df,binNumber==i)
    site_bin_positive[i,2] = as.numeric(count(filter(dataframe, dataframe$is_pos==1)))
    site_bin_negative[i,2] = as.numeric(count(filter(dataframe, dataframe$is_pos==-1)))
    
  }
  site_bin_negative[,1] = c(seq(1,B))
  site_bin_positive[,1] = c(seq(1,B))
  
    #site level density
  site_level_densities = matrix(NA,ncol=1,nrow = B)
  for(i in 1:B){
    site_level_densities[i,1] = rbeta(1,
                                      site_bin_positive[i,2]+beta_prior,
                                      site_bin_negative[i,2]+beta_prior)*site_bin_weight[i]
  }
  site_level_density = sum(site_level_densities)
  
  #Study level: expected P(b|+) and P(b|-)
  alpha = study_fit$beta_parameters[,1]             # = k_pos + c
  beta = study_fit$beta_parameters[,2]              # = k_neg + c
  E_pos_given_b = alpha / (alpha + beta)           # E[P(+|b)]
  P_b = as.vector(study_fit$Study_level_weights)    # study P(b)
  P_b = pmax(P_b, eps)
  P_b = P_b / sum(P_b)  #Also to make sure in case < eps came true
  
  P_pos = sum(E_pos_given_b * P_b)
  P_b_pos = (E_pos_given_b * P_b) / P_pos              # P(+|b)P(b)/P(+)
  P_b_neg = ((1 - E_pos_given_b) * P_b) / (1 - P_pos)  # P(-|b)P(b)/(1-P(+))
  
  P_b_pos = pmax(P_b_pos, eps)
  P_b_pos = P_b_pos / sum(P_b_pos) #in case < eps was true
  P_b_neg = pmax(P_b_neg, eps)
  P_b_neg = P_b_neg / sum(P_b_neg) #in case < eps was true
  
  
  #Grid search q in (0,1) to minimize KL( P_s(b) || q P(b|+) + (1-q) P(b|-) )
  Ps_b = P_b_pos*site_level_density + P_b_neg*(1-site_level_density)
  
  q_seq = seq(0.001, 0.999, length.out = 2001)   # avoid exact 0/1
  KL_vals = sapply(q_seq, function(q) {
    Q = q * P_b_pos + (1 - q) * P_b_neg
    Q = pmax(Q, eps)
    Q = Q / sum(Q) # in case < eps was true
    sum(Ps_b * (log(Ps_b) - log(Q)))
  })
  return(Stat2_density_estimate = q_seq[which.min(KL_vals)])
}

bootstrap_strat2_site = function(site_examples, study_fit,
                                  Bdraws = 2000, 
                                  conf = 0.95, 
                                  eps = 1e-12) {
  
  #site_examples = fdata
  #study_fit = out
  #eps = 1e-12
  
  n = length(site_examples)
  
  point = Strat2_call_density(site_examples, 
                              study_fit = study_fit,
                              eps = eps)
  
 
  
  boots = Bdraws
  for (b in seq_len(Bdraws)) {
    idx = sample.int(n, n, replace = TRUE)
    boots[b] = Strat2_call_density(site_examples[idx],
                                   study_fit = study_fit,
                                   eps = eps)
  }
  
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(boots, probs = c(alpha, 1 - alpha)))
  return(list(point = point, se = sd(boots),
              conf_int = ci, 
              samples = boots, 
              conf = conf))
}


Strat2_call_density(fdata,study_fit = out,eps = 1e-12)
out2 = bootstrap_strat2_site(fdata,out,Bdraws = 2000, conf = 0.95, eps = 1e-12)

hist(out2$samples)
out2$conf_int
out2$point
out$Density_Estimate



#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(density1, density2){
  density_estimate = (density1 * density2)^(1/2)
  return (Density_estimate = density_estimate)
}

#it works need to figure out why and how
Strat3_call_density(Strat1_call_density(fdata,
                                        10000,
                                        0.1)$Density_Estimate,
                    Strat2_call_density(fdata,study_fit = out,eps = 1e-12))



#Bootstrapping and confidence intervals
bootstrap_strat3_site = function(site_examples,
                                 study_fit,
                                 Bdraws = 2000,
                                 conf = 0.95,
                                 eps = 1e-12) {
  
  # --- deterministic S1 (no randomness) ---
  strat1_point_internal = function(site_examples, study_fit, eps) {
    
    site_df = as.data.frame(do.call(rbind, site_examples), stringsAsFactors = FALSE)
    names(site_df)[names(site_df) == "bin"] = "binNumber"
    site_df$binNumber = as.numeric(site_df$binNumber)
    
    B = nrow(study_fit$beta_parameters)
    Ps_b = sapply(1:B, function(b) mean(site_df$binNumber == b))
    Ps_b = pmax(Ps_b, eps)
    Ps_b = Ps_b / sum(Ps_b)
    
    alpha = study_fit$beta_parameters[,1]
    beta  = study_fit$beta_parameters[,2]
    E_pos_given_b = alpha / (alpha + beta)
    
    sum(E_pos_given_b * Ps_b)
  }
  
  # --- Strategy 3 point estimate ---
  s1 = strat1_point_internal(site_examples, study_fit, eps)
  s2 = Strat2_call_density(site_examples, study_fit = study_fit, eps = eps)
  point = sqrt(s1 * s2)
  
  # --- Bootstrap samples ---
  n = length(site_examples)
  boots = numeric(Bdraws)
  for (b in seq_len(Bdraws)) {
    idx = sample.int(n, n, replace = TRUE)
    s1_b = strat1_point_internal(site_examples[idx],
                                 study_fit,
                                 eps)
    s2_b = Strat2_call_density(site_examples[idx],
                               study_fit = study_fit,
                               eps = eps)
    boots[b] = sqrt(s1_b * s2_b)
  }
  
  # --- Percentile CI ---
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(boots, probs = c(alpha, 1 - alpha)))
  
  list(point = point, samples = boots, se = sd(boots), conf_int = ci, conf = conf)
}

out3 = bootstrap_strat3_site(fdata,out,Bdraws = 2000, conf = 0.95, eps = 1e-12)
hist(out3$samples)
out3$point











