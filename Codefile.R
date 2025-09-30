
library(readr)
library(tidyverse)
library(ggplot2)

################################################################################
#                                 The datasets                                 #
################################################################################

#-----Chunk 1 - generating fake data--------------------------------------------
#generate fake data 1
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

generate_fake_validation_data = function(n = 100, bins = 5) {
  scores = rexp(n,rate =6)
  scores = scores/max(scores)                                   # to make sure it <= 1
  quantile_bounds = seq(0,1,length.out = bins + 1)
  
  bin_indices = findInterval(scores,
                             quantile_bounds,
                             rightmost.closed = TRUE)           #returns the bin number for each score
  
  bin_weights = table(bin_indices) / n                          #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)]) #gives the bin proportion associated with each score
  
  examples = vector("list",n)
  
  for (i in seq_len(n)) {
    ex = new_validation_example(
      filename = paste0("file_", sample(1:10, 1), ".wav"),
      timestamp_offset = runif(1, 0, 60),
      score = scores[i],
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.10, 0.8, 0.1)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] = ex
  }
  return(examples)
}

generate_fake_validation_data2 = function(n = 100, bins = 5) {
  scores = rexp(n,rate = 6)
  scores = scores/max(scores)                                           # to makes sure it is <= 1
  quantile_bounds = seq(0,1,length.out = bins + 1)
  
  bin_indices = findInterval(scores,
                             quantile_bounds,
                             rightmost.closed = TRUE)                   #returns the bin number for each score
  
  bin_weights = table(bin_indices) / n                                  #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)])  #gives the bin proportion associated with each score
  
  examples = vector("list",n)
  
  for (i in seq_len(n)) {
    ex <- new_validation_example(
      filename = paste0("file_", sample(1:10, 1), ".wav"),
      timestamp_offset = runif(1, 0, 60),
      score = scores[i],
      is_pos = sample(c(1, -1, 0), 1, prob = c(0.10, 0.80, 0.10)),
      bin = bin_indices[i],
      bin_weight = bin_weights_vec[i]
    )
    examples[[i]] <- ex
  }
  return(examples)
}

#Fake datasets
fdata = generate_fake_validation_data(n = 1000,bins= 5)
fdata2 = generate_fake_validation_data2(n = 1000,bins= 5)
study_level_data = rbind(fdata,fdata2)

#------------------------------ real data set ----------------------------------
validation_data = read_csv("Abyssinian Nightjar validation results.csv") 

validation_data1 = read_csv("African Dusky Flycatcher validation results Grid 2.csv")



validation_data$outcome = sample(c(1,-1,0),nrow(validation_data), replace = T)
validation_data = rename(validation_data,score = confidence)


#--------visualization of the data set -------------------------------------
Confidence.histogram = ggplot(validation_data, aes(x = score))+
  geom_histogram(binwidth = 0.060, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of Confidence Scores",
    x = "Confidence",
    y = "Frequency"
    ) +
  theme_minimal()

Logit.Confidence.histogram = ggplot(validation_data, aes(x = logit_conf))+
  geom_histogram(binwidth = 0.60, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of logit Confidence Scores",
    x = "logit Confidence",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("_output/confidence_histogram.jpeg",
       plot = Confidence.histogram, 
       width = 6, 
       height = 4, 
       dpi = 300)

ggsave("_output/Logit Confidence histogram.jpeg", 
       plot = Logit.Confidence.histogram, 
       width = 6, 
       height = 4, 
       dpi = 300)



################################################################################
#                                 Quantiling the dataset                       #
################################################################################
#---function to return quantiles------------------------------------------------
Quantile_function = function(examples){
  quantiles = c(0.5,0.75,0.875)
  binning.values = c()
  
  #examples = validation_data
  logits = sort(examples$logit_conf)
  for (i in 1:length(quantiles)){
    binning.values = c(binning.values,quantile(logits,quantiles[i]))
  }
  return(Quantiles = binning.values)
}
Bin.bounds = Quantile_function(validation_data)

#bins visualized on real dataset
binned.Logit.Confidence.histogram = ggplot(validation_data, aes(x = logit_conf))+
  geom_histogram(binwidth = 0.60, fill = "maroon", colour = "white")+
  geom_vline(xintercept = c(Bin.bounds), colour = "black", linewidth = 0.9)+ 
  annotate("text",
           x = c(Bin.bounds), 
           y = 50,             # adjust this based on your y-axis scale
           label = round(c(Bin.bounds),3),
           angle = 90, vjust = -0.5, hjust = +4.9) +
  labs(
    title = "Binned Histogram of logit Confidence Scores",
    x = "logit confidence",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("Binned logit confidence_histogram.jpeg", plot = binned.Logit.Confidence.histogram, width = 6, height = 4, dpi = 300)




################################################################################
#                        Study level density estimation                        #
################################################################################

Study_level_call_density = function(examples,
                                    num_beta_samples=1000,
                                    small_constant=0.1,
                                    bins = 5){
  
  #loading in of examples data and converting it to usable format
#-------------------binning and weighting---------------------------------------
  examples = validation_data
  quantile_bounds = Bin.bounds
  
  
  nrow(filter(examples,examples$outcome==0))
  
  
  
  bin_indices = findInterval(examples$logit_conf ,
                             quantile_bounds,
                             rightmost.closed = TRUE)+1                   #returns the bin number for each score
  
  bin_weights = table(bin_indices) / nrow(validation_data)   #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)])  #gives the bin proportion associated with each score
  examples = cbind(examples, bin_indices,bin_weights_vec)
  
  examples = rename(examples,binNumber=bin_indices)     

 
#--------------------- Binning and bin weighting -------------------------------
  
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
  

#--------------------------Density  estimation----------------------------------
  # Create beta distributions for each bin.
  shape.parameters = matrix(NA,ncol=2,nrow=nrow(bin_number))
  
  for (b in bin_number){
    shape.parameters[b,] = c(bin_positive[b,2] + small_constant,
                             bin_negative[b,2] + small_constant)
  }
  
  #estimating the density
  density_matrix = matrix(NA,nrow = nrow(bin_number),ncol = 1)
  
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = ((shape.parameters[i,1]+ small_constant)/
      ((shape.parameters[i,1]+ small_constant)+(shape.parameters[i,2] + small_constant)))*bin_weights[i,1]
  }
  
  density = sum(density_matrix)
    
  
#-------------Bootstrapping: distribution of density estimate-------------------
  
  NoBetaSamples = 10000
  q_betas = numeric(NoBetaSamples)
  
  for (i in seq_len(NoBetaSamples)) {
    betas = rbeta(n = nrow(shape.parameters),
      shape.parameters[,1] + small_constant, 
      shape.parameters[,2] + small_constant)

    q_betas[i] = sum(betas * bin_weights[,1])
  }
  
#--------------------------function output -------------------------------------
  
  return(list(Density_Estimate = density,
              Samples = q_betas, 
              beta_parameters = shape.parameters,
              Study_level_weights = bin_weights,
              densities = density_matrix,
              negative_counts = bin_negative,
              positive_counts = bin_number))
}

########################## function implementation #############################

out = Study_level_call_density(study_level_data,1000,0.1,5)  
out$Density_Estimate
out$Samples
#--------distribution of bootstrap samples and density estimates----------------

histogram.bootstrap.estimates = ggplot(data.frame(out$Samples), aes(x = out.Samples))+
  geom_histogram(binwidth = 0.008, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out$Density_Estimate, colour = "black", linewidth = 0.9)+
  annotate("text",
           x = out$Density_Estimate, 
           y = 50,             # adjust this based on your y-axis scale
           label = paste0("Density Estimate = ", round(out$Density_Estimate, 3)),
           angle = 90, vjust = -0.5, hjust = 0) +
  labs(
    title = "Histogram of density bootstrap estimates",
    x = "Density estimates",
    y = "Frequency") +
  theme_minimal()

ggsave("Historgram of bootstrap estimates study.jpeg",
       plot = histogram.bootstrap.estimates, 
       width = 6, 
       height = 4,
       dpi = 300)


#--------------------95% bootstrap CI, std and variance  -----------------------

confidence_interval = quantile(out$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out$Samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)

cat("95% Confidence interval: [",round(confidence_interval,4),"]",
    "\nvariance: ",round(Variance_bootstrap_samples,4),
    "\nstandard deviation: ",round(Sd_bootstrap_samples,4))




################################################################################
#                Strategy 1 strata level density estimation                    #
################################################################################
# ---------------------Distribution shifts Strategy 1---------------------------


#sum(P(+|b)*P_s(b))
Strat1_call_density = function(examples,
                                num_beta_samples=1000,
                                small_constant=0.1){

#-----------Converting data set to a usable format------------------------------
  quantile_bounds = Bin.bounds
  
  
  bin_indices = findInterval(examples$logit_conf ,
                             quantile_bounds,
                             rightmost.closed = TRUE)+1                   #returns the bin number for each score
  
  bin_weights = table(bin_indices) / nrow(validation_data)   #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)])  #gives the bin proportion associated with each score
  examples = cbind(examples, bin_indices,bin_weights_vec)
  examples = rename(examples,binNumber=bin_indices)     
  
  
  #Positive/negative bin counts
  bin_number  = as.matrix(distinct(examples,binNumber)) 
  bin_number = as.matrix(sort(as.numeric(as.vector(bin_number[,1]))))
  bin_weights = matrix(NA, nrow = nrow(bin_number),ncol=1)
  
  #Bin weights
  for(i in 1:nrow(bin_number)){
    bin_weights[i,1]=nrow(filter(examples,binNumber==i))/nrow(examples)
  }

#----------------strata level density estimation--------------------------------
  density_matrix = matrix(NA,nrow=nrow(bin_number), ncol =1)
 
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
    ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*bin_weights[i,1]
  }
  
  density = sum(density_matrix)

  
#-----------------------------Bootstrapping distribution -----------------------
  NoBetaSamples = 10000
  q_betas = numeric(NoBetaSamples)
  
  for (i in seq_len(NoBetaSamples)) {
    # draw a bin-level theta_b from the same Beta as the estimator
    betas = rbeta(
      n = nrow(out$beta_parameters),
      shape1 = out$beta_parameters[,1] + small_constant,  
      shape2 = out$beta_parameters[,2] + small_constant
    )
    
    q_betas[i] = sum(betas * bin_weights[,1])
  }
  
#----------------------------------function output------------------------------
  return(list(Density_Estimate = density,
              Samples = q_betas,
              bin_densities = density_matrix))
}


########################## function implementation #############################
#Implementation of strategy 1 call density  

out1 = Strat1_call_density(validation_data,1000,0.1)
out1$Density_Estimate

#---------------Implementation of study_level call density----------------------  

histogram.bootstrap.estimates.c1= ggplot(data.frame(out1$Samples), aes(x = out1.Samples))+
  geom_histogram(binwidth = 0.008, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out1$Density_Estimate, colour = "black", linewidth = 0.9)+
  labs(
    title = "Histogram of density bootstrap estimates cov 1",
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("Historgram of bootstrap estimates cov 1.jpeg", 
       plot = histogram.bootstrap.estimates.c1,
       width = 6,
       height = 4,
       dpi = 300)


#--------------------95% bootstrap CI, std and variance  -----------------------
confidence_interval = quantile(out1$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out1$Samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)

cat("95% Confidence interval: [",round(confidence_interval,4),"]",
    "\nvariance: ",round(Variance_bootstrap_samples,4),
    "\nstandard deviation: ",round(Sd_bootstrap_samples,4))


################################################################################
#                Strategy 2 strata level density estimation                    #
################################################################################
#---------Strategy 2------------------------------------------------------------
Strat2_call_density = function(site_examples, 
                               study_fit = out,
                               eps = 1e-12,small_constant= 0.1) {
  

#---------------------------Binning data and format-----------------------------

  quantile_bounds = Bin.bounds
  
  bin_indices = findInterval(site_examples$logit_conf ,
                             quantile_bounds,
                             rightmost.closed = TRUE)+1                         #returns the bin number for each score
  
  bin_weights = table(bin_indices) / nrow(validation_data)                      #gives the proportions in the different bins
  bin_weights_vec = as.numeric(bin_weights[as.character(bin_indices)])          #gives the bin proportion associated with each score
  site_examples = cbind(site_examples, bin_indices,bin_weights_vec)
  site_examples = rename(site_examples,binNumber=bin_indices)     
  
  
  B = nrow(study_fit$beta_parameters)
  site_bin_weight = sapply(1:B, function(b) mean(as.numeric(site_examples$binNumber) == b)) #Site bin weights
  site_bin_weight[site_bin_weight < eps] = eps                                        # to avoid 0s when calculating things like log
  site_bin_weight = site_bin_weight / sum(site_bin_weight)                            #proportioning to it to make sure it sums to 1 after some may be set to 0
  
  site_bin_positive = matrix(NA,ncol = 2, nrow = B)
  site_bin_negative = matrix(NA,ncol = 2, nrow = B)
  
  for(i in 1:B){
    dataframe = filter(site_examples,binNumber==i)
    site_bin_positive[i,2] = as.numeric(count(filter(dataframe, dataframe$outcome==1)))
    site_bin_negative[i,2] = as.numeric(count(filter(dataframe, dataframe$outcome==-1)))
  }
  
  site_bin_negative[,1] = c(seq(1,B))
  site_bin_positive[,1] = c(seq(1,B))
  
#--------------------------site density estimation -----------------------------
  
  site_level_densities = matrix(NA,ncol=1,nrow = B)
  
  for(i in 1:B){
    site_level_densities[i,1] = ((site_bin_positive[i,2]+small_constant)/
    ((site_bin_positive[i,2]+small_constant)+(site_bin_negative[i,2]+small_constant)))*site_bin_weight[i]
  }
  
  site_level_density = sum(site_level_densities)
  
  #Study level: expected P(b|+) and P(b|-)
  alpha = study_fit$beta_parameters[,1]             # = k_pos + c
  beta = study_fit$beta_parameters[,2]              # = k_neg + c
  E_pos_given_b = alpha / (alpha + beta)            # E[P(+|b)]
  P_b = as.vector(study_fit$Study_level_weights)    # study P(b)
  P_b = pmax(P_b, eps)
  P_b = P_b / sum(P_b)                              # Also to make sure in case < eps came true
  
  P_pos = sum(E_pos_given_b * P_b)
  P_b_pos = (E_pos_given_b * P_b) / P_pos              # P(+|b)P(b)/P(+)
  P_b_neg = ((1 - E_pos_given_b) * P_b) / (1 - P_pos)  # P(-|b)P(b)/(1-P(+))
  
  P_b_pos = pmax(P_b_pos, eps)
  P_b_pos = P_b_pos / sum(P_b_pos)                     # in case < eps was true
  P_b_neg = pmax(P_b_neg, eps)
  P_b_neg = P_b_neg / sum(P_b_neg) 

  
#--------------------------gird iteration/search--------------------------------
  
  Ps_b = P_b_pos*site_level_density + P_b_neg*(1-site_level_density)
    
  q_seq = seq(0.001, 0.999, length.out = 2001)                # avoid exact 0/1
  KL_vals = sapply(q_seq, function(q) {
    Q = q * P_b_pos + (1 - q) * P_b_neg
    Q = pmax(Q, eps)
    Q = Q / sum(Q)                                            # in case < eps was true
    sum(Ps_b * (log(Ps_b) - log(Q)))
  })
  return(list(Stat2_density_estimate = q_seq[which.min(KL_vals)],
         KL_values = KL_vals,
         qs = q_seq))
}

#-----------------------------Bootstrapping ------------------------------------

bootstrap_strat2_site = function(site_examples, study_fit ,Bdraws = 2000,conf = 0.95,eps = 1e-12) {
  
  #site_examples=validation_data
  n = nrow(site_examples)
  point = Strat2_call_density(site_examples, study_fit = study_fit, eps = eps)$Stat2_density_estimate
  boots = Bdraws
  
  for (b in seq_len(Bdraws)) {
    idx = sample.int(n, n, replace = TRUE)
    boots[b] = Strat2_call_density(site_examples[idx,], study_fit = study_fit, eps = eps)$Stat2_density_estimate
  }
  
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(unlist(boots), probs = c(alpha, 1 - alpha)))
  
  
  return(list(Density_estimate = point, se = sd(boots),
              conf_int = ci, 
              Samples = boots, 
              conf = conf))
}


########################## function implementation #############################

out2.density = Strat2_call_density(validation_data,study_fit = out,eps = 1e-12)
out2 = bootstrap_strat2_site(validation_data,out,Bdraws = 2000, conf = 0.95, eps = 1e-12)

#---------implementation of strategy 1 density estimation (site)----------------

histogram.bootstrap.estimates.s2= ggplot(data.frame(out2$Samples), aes(x = out2.Samples))+
  geom_histogram(binwidth = 0.009, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out2.density$Stat2_density_estimate, colour = "black", linewidth = 0.9)+
  labs(
    title = "Histogram of density bootstrap estimates s2",
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("Historgram of bootstrap estimates s2.jpeg", 
       plot = histogram.bootstrap.estimates.s2, 
       width = 6,
       height = 4,
       dpi = 300)

#-------------------------------kl divergence plot------------------------------

kl.data = cbind(out2.density$qs,out2.density$KL_values)
kl.plot = ggplot(data = data.frame(kl.data), aes(x = X1,y = X2))+
  geom_line()+
  geom_vline(xintercept = out2.density$Stat2_density_estimate, colour = "red")+
  labs(title = "Plot of KL divergence",
  x = "q values",
  y = "KL divergence") +
  theme_minimal()

ggsave( "KL divergence.jpeg", plot = kl.plot, width = 6, height = 4, dpi = 300)


################################################################################
#                Strategy 3 strategy level density estimation                  #
################################################################################
#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(density1, density2){
  density_estimate = (density1 * density2)^(1/2)
  
  return (Density_estimate = density_estimate)
}



#-------------------Bootstrapping and confidence intervals----------------------
bootstrap_strat3_site = function(site_examples,study_fit,Bdraws = 2000,conf = 0.95,eps = 1e-12) {
  
# -------------------------------- Bootstrap samples ---------------------------
  #site_examples = validation_data
  study_fit = out
  
  n = nrow(site_examples)
  boots = numeric(Bdraws)
  for (b in seq_len(Bdraws)) {
    s1 = Strat1_call_density(site_examples,10000,0.1)$Density_Estimate
    s2 = Strat2_call_density(site_examples,study_fit = out,eps = 1e-12, small_constant = 0.1)$Stat2_density_estimate
    boots[b] = sqrt(s1 * s2)
  }
  
  
#------------------------95% Confidence interval--------------------------------
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(boots, probs = c(alpha, 1 - alpha)))
 
  
#------------------------output of the function---------------------------------    
  return(list(samples = boots, se = sd(boots), conf_int = ci, conf = conf))
}

########################## function implementation #############################

out3.1 = Strat3_call_density(Strat1_call_density(validation_data,num_beta_samples=1000,
                                                 small_constant=0.1)$Density_Estimate,
                    Strat2_call_density(validation_data,study_fit = out,
                                        eps = 1e-12, small_constant = 0.1)$Stat2_density_estimate)

out3 = bootstrap_strat3_site(validation_data,out,Bdraws = 200,conf = 0.95,eps = 1e-12)

histogram.bootstrap.estimates.s3= ggplot(data.frame(out3$samples), aes(x = out3.samples))+
  geom_histogram(binwidth = 0.005, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out3.1, colour = "black",linewidth = 0.9)+
  labs(
    title = "Histogram of density bootstrap estimates s3",
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("Historgram of bootstrap estimates s3.jpeg",
       plot = histogram.bootstrap.estimates.s3, 
       width = 6,
       height = 4,
       dpi = 300)







