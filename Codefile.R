
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
generate_fake_validation_data <- function(n = 100, bins = 5) {

  scores = runif(n, 
                 min = 0,
                 max = 1)
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
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.4, 0.5, 0.1)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] <- ex
  }
  return(examples)
}

generate_fake_validation_data2 <- function(n = 100, bins = 5) {
  
  scores = runif(n, 
                 min = 0,
                 max = 1)
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
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.55, 0.35, 0.1)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] <- ex
  }
  return(examples)
}

#debug(generate_fake_validation_data)
fdata = generate_fake_validation_data(n = 100,
                                      bins= 5)
fdata2 = generate_fake_validation_data2(n = 100,
                                       bins= 5)

study_level_data = rbind(fdata,fdata2)


# -------Chunk 2 - Call density estimation-------------------------------------
#the algorithm will work as long it is supplied with data that is not binned and weighted but is validated binning and weighting happen within the function.
Study_level_call_density = function(examples,num_beta_samples=1000,beta_prior=0.1){
  #examples = Data2
  examples = unlist(examples)
  df <- as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE) #This will be the one that i check my bin and bin_weight generation  
  examples = as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE)
  examples = rename(examples,binNumber=bin, outcome = is_pos)

  #infinix operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b 
  
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
    density_matrix[i,1] = rbeta(1,shape.parameters[i,1]+beta_prior,
                             shape.parameters[i,2]+beta_prior)*bin_weights[i,1]
  }
  
  density = sum(density_matrix)
    
#Bootstrapping: To get distribution of density estimate
  q_betas=c()
  NoBetaSamples = 10000
  for(i in 1:NoBetaSamples){
    for(j in 1:nrow(bin_number)){
      rvalue = 0
      rvalue = rvalue + rbeta(1,shape.parameters[j,1],
                              shape.parameters[j,2])
    }
    q_betas=c(q_betas,rvalue)
  }
  return(list(Density_Estimate = density,
              Samples = q_betas, 
              beta_parameters = shape.parameters,
              Study_level_weights = bin_weights,
              densities = density_matrix))
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
Strat1_call_density = function(examples, num_beta_samples=1000, beta_prior=0.1){
  
  examples = unlist(examples)
  df <- as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE) #This will be the one that i check my bin and bin_weight generation  
  examples = as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE)
  examples = rename(examples,binNumber=bin, outcome = is_pos)

  #infinix operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b 
  
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
  density = 0 
  
  for(i in 1:nrow(out$beta_parameters)){
    density = density + rbeta(1,out$beta_parameters[i,1]+beta_prior,
                              out$beta_parameters[i,2]+beta_prior)*bin_weights[i,1]
  }

    
#Bootstrapping: To get distribution of density estimate
  q_betas=c()
  NoBetaSamples = 10000
  for(i in 1:NoBetaSamples){
    for(j in 1:nrow(bin_number)){
      rvalue = 0
      rvalue = rvalue + rbeta(1,out$beta_parameters[j,1],
                              out$beta_parameters[j,2])
    }
    q_betas=c(q_betas,rvalue)
  }
  return(list(Density_Estimate = density,
              Samples = q_betas))
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
#P(b|+)
Strat2_call_density = function(examples, num_beta_samples = 1000,beta_prior = 0.1){
  examples = unlist(examples)
  df <- as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE) #This will be the one that i check my bin and bin_weight generation  
  examples = as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE)
  examples = rename(examples,binNumber=bin, outcome = is_pos)
  
  #infinix operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b 
  
  #Positive/negative bin counts
  bin_number  = as.matrix(distinct(examples,binNumber)) 
  bin_number = as.matrix(sort(as.numeric(as.vector(bin_number[,1]))))
  bin_weights = matrix(NA, nrow = nrow(bin_number),ncol=1)
  
  
  #Bin weights
  for(i in 1:nrow(bin_number)){
    bin_weights[i,1]=nrow(filter(examples,binNumber==i))/nrow(examples)
  }
  
  #Density  calculations
  #Calculation of P(+|b)P(b)/P(+)
  b_given_pos = out$densities/out$Density_Estimate
  
  #Calculation of P(-|b)P(b)/P(-)
  neg_densities = matrix(NA,nrow = nrow(bin_number),
                         ncol =1)
  
  for( i in 1:nrow(bin_number)){
    neg_densities[i,1] = rbeta(1,out$beta_parameters[i,2] + beta_prior,
                              out$beta_parameters[i,1] + beta_prior)
  }
  
  b_given_neg = neg_densities / sum(neg_densities)
  
  
  #Density  
  # Create beta distributions for each bin.
  shape.parameters_site = matrix(NA,ncol=2,nrow=nrow(bin_number))
  
  for (b in bin_number){
    shape.parameters_site[b,] = c(
      bin_positive[b,2] + beta_prior, bin_negative[b,2] + beta_prior
    )
  }
  
  
  # site density estimation P_s(+)
  shape.parameters
  density_matrix = matrix(NA,nrow = nrow(bin_number),ncol = 1)
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = rbeta(1,shape.parameters[i,1]+beta_prior,
                                shape.parameters[i,2]+beta_prior)*bin_weights[i,1]
  }
  
 # P_s(+)
  prob_b_s = b_given_pos*density_matrix + b_given_neg * (1- density_matrix)
  
  
  #Q_q(b)
  q_seq = seq(0.0001,0.9999,length.out = 100)
  
  Q_q_b = b_given_pos*q_seq[1] + b_given_neg * (1-q_seq[1])
  KL_div.old = sum(prob_b_s*(log(prob_b_s)-log(Q_q_b)))
  Site_density = NA
  
  
  for (i in 2:length(q_seq)){
    Q_q_b = b_given_pos*q_seq[i] + b_given_neg * (1-q_seq[i])
    KL_div = sum(prob_b_s*(log(prob_b_s)-log(Q_q_b)))
    if (KL_div <=KL_div.old){
      KL_div.old = KL_div
      Site_density = q_seq[i]
    }
  }
  
  
  return(Density_Estimate = Site_density)
}

Strat2_call_density(fdata,1000,0.1)

#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(density1, density2){
  density_estimate = (density1 * density2)^(1/2)
  return (Density_estimate = density_estimate)
}

Strat3_call_density(Strat1_call_density(fdata,
                                        10000,
                                        0.1)$Density_Estimate,
                    Strat2_call_density(fdata,
                                        10000,
                                        0.1))




















