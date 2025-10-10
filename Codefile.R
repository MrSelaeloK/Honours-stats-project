
library(readr)
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(grid)

################################################################################
#                                 The datasets                                 #
################################################################################

#-----------------------------The datasets loading--------------------------------------
# === SETTINGS ===
folder_path <- "32KHz_raw_data"   # 🔹 Change this
output_file <- "combinedhz32Khz.csv"

# === STEP 1: Get all .txt files ===
files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)

# === STEP 2: Read each file safely ===
read_file_safe <- function(f) {
  message("Reading file: ", f)
  df <- read_delim(f, delim = "\t", col_types = cols(.default = "c"))
  df$source_file <- basename(f)
  return(df)
}

data_all <- map(files, read_file_safe)

# === STEP 3: Combine all files ===
combined <- bind_rows(data_all)

#converting the columns to the right data
combined$Confidence = as.double(combined$Confidence)
write_csv(combined, output_file)


#------------------------------ real data set ----------------------------------
strata1 = read_csv("combinedhz48Khz.csv") 
strata2 = read_csv("combinedhz32Khz.csv")

strata1 = rename(strata1,score = Confidence)
strata2 = rename(strata2,score = Confidence)
strata1$SamplingRate = 48
strata2$SamplingRate = 32

study_level_data = rbind(strata1,strata2)
study_level_data$logit_conf = log(study_level_data$score) 


#######################  choosing species  #####################################
species = function(species_name){
  study_level_data = filter(study_level_data,`Common Name` == species_name)
  strata1 = filter(strata1, `Common Name` == species_name)
  strata2 = filter(strata2, `Common Name` == species_name)
  
  return(list(study_level_data = study_level_data,
              strata1 = strata1,
              strata2 = strata2))
}

species_choice = species("African Pipit")


#--------visualization of the data set -------------------------------------
Confidence.histogram = ggplot(species_choice$study_level_data, aes(x = score))+
  geom_histogram(binwidth = 0.060, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of all Confidence Scores - African Pipit",
    x = "Confidence",
    y = "Frequency"
    ) +
  theme_minimal()

Logit.Confidence.histogram = ggplot(species_choice$study_level_data, aes(x = logit_conf))+
  geom_histogram(binwidth = 0.190, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of all logit Confidence Scores- African Pipit",
    x = "logit Confidence",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("_output/confidence_histogram_African Pipit.jpeg",
       plot = Confidence.histogram, 
       width = 6, 
       height = 4, 
       dpi = 300)

ggsave("_output/Logit Confidence histogram_African Pipit.jpeg", 
       plot = Logit.Confidence.histogram, 
       width = 6, 
       height = 4, 
       dpi = 300)



################################################################################
#                                 Quantiling the dataset  by logits            #
################################################################################
quantiles.3 = c(0.60,0.90)

#quantiles.4 = c(0.3,0.25,0.125,0.125)
#quantiles.5 = c(0.3,0.25,0.20,0.15,0.10)
#quantiles.6 = c(0.3, 0.275, 0.175, 0.15, 7.5, 2.5)


Quantile_function = function(examples,quantiles){
  binning.values = c()                                  #initial empty binning dataset
  
  examples = species_choice$study_level_data
  logits = sort(examples$logit_conf)
  for (i in 1:length(quantiles)){
    binning.values = c(binning.values,quantile(logits,quantiles[i]))
  }
  return(Quantiles = binning.values)
}
Bin.bounds = Quantile_function(species_choice$study_level_data,quantiles.3)

#bins visualized on real dataset
binned.Logit.Confidence.histogram = ggplot(species_choice$study_level_data,aes(x = logit_conf))+
  geom_histogram(binwidth = 0.10, fill = "maroon", colour = "white")+
  geom_vline(xintercept = c(Bin.bounds), colour = "black", linewidth = 0.9)+ 
  annotate("text",
           x = c(Bin.bounds), 
           y = 50,             # adjust this based on your y-axis scale
           label = round(c(Bin.bounds),3),
           angle = 90, vjust = -0.5, hjust = +4.9) +
  labs(
    title = "Binned Histogram of logit Confidence Scores_African Pipit",
    x = "logit confidence",
    y = "Frequency"
  ) +
  theme_minimal()

ggsave("_output/Binned logit confidence_histogram_African Pipit.jpeg",
       plot = binned.Logit.Confidence.histogram,
       width = 6, 
       height = 4,
       dpi = 300)



################################################################################
#                            Validation dataset                                #
################################################################################
#-------------------------validation of the dataset-----------------------------
x <- sample(c(-1,0, 1), 
            size = nrow(species_choice$study_level_data),
            replace = TRUE, prob = c(0.1, 0.3, 0.6))

species_choice$study_level_data$outcome = x

#-------------------------Getting the samples per bin---------------------------
#Size = 25
validation.data = function(Bin.bounds,Size, dataset){

bounds = c(min(species_choice$study_level_data$logit_conf),Bin.bounds, max(species_choice$study_level_data$logit_conf))
#Size = 25

Validated_data = data.frame()

for (i in 2:length(bounds)) {
  lower = bounds[i - 1]
  upper = bounds[i]
  
  # Proper two-sided condition: use "&"
  #The dataset parameter comes from the parameters and is used to pick study level or strata level 
  subset_rows = subset(dataset, logit_conf > lower & logit_conf <= upper)
  if(nrow(subset_rows)<Size){
    sampled = subset_rows[sample(nrow(subset_rows), Size, replace = T), ]
  }else{
    sampled = subset_rows[sample(nrow(subset_rows), Size), ]
  }
  Validated_data = rbind(Validated_data,sampled)
}
return(Validated_data)
}
validated.data = validation.data(Bin.bounds,Size = 25, species_choice$study_level_data)


################################################################################
#                        Study level density estimation                        #
################################################################################

Study_level_call_density = function(examples,
                                    num_beta_samples = 1000,
                                    small_constant = 0.1,
                                    quantile_bounds){
  
  #loading in of examples data and converting it to usable format
#-------------------binning and weighting---------------------------------------
  #examples = validated.data
  quantile_bounds = Bin.bounds   # e.g. c("60%" = -0.93496400, "90%" = -0.09486028)
  #bounds come from the unvalidated data
  bounds = as.numeric(c(min(species_choice$study_level_data$logit_conf),
             Bin.bounds,
             max(species_choice$study_level_data$logit_conf)))
  
  #counts in each bin
  pos.neg.neu = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu) <- paste0("Bin", 1:(length(bounds) - 1))
  pos.neg.neu
  
  #allocating to pos.neg.neu
  for(i in 2:length(bounds)){
    subset_rows = subset(examples,
                          logit_conf >= bounds[i - 1] &
                            logit_conf <= bounds[i])
    
      pos.neg.neu["Positive",i-1] = sum(subset_rows$outcome==1)
      pos.neg.neu["Negative",i-1] = sum(subset_rows$outcome==-1)
      pos.neg.neu["Unknown",i-1] = sum(subset_rows$outcome==0)
  }
  
  #bin weights
  bin.weights = matrix(NA, nrow = ncol(pos.neg.neu), ncol = 1)
  rownames(bin.weights) = paste0("Bin", 1:(length(bounds) - 1))
  colnames(bin.weights) = "weights"
  bin.weights
  
  proportions.pos.neg.neu = rowSums(pos.neg.neu)/sum(rowSums(pos.neg.neu))
  matrix.proportions.row = pos.neg.neu/rowSums(pos.neg.neu)
  
  for(i in 1:ncol(pos.neg.neu)){
    bin.weights[i,] = matrix.proportions.row[1,i]*proportions.pos.neg.neu[1] +
      matrix.proportions.row[2,i]*proportions.pos.neg.neu[2]
  }


#--------------------------Density  estimation----------------------------------
  
  # Create beta distributions for each bin.
  shape.parameters = matrix(NA,ncol=2,nrow=ncol(pos.neg.neu))
  
  for (i in 1:ncol(pos.neg.neu)){
    shape.parameters[i,] = c(pos.neg.neu[1,i] + small_constant,
                             pos.neg.neu[2,i] + small_constant)
  }
  
  #estimating the density
  density_matrix = matrix(NA,nrow = nrow(bin.weights),ncol = 1)
  
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = ((shape.parameters[i,1]+ small_constant)/
      ((shape.parameters[i,1]+ small_constant)+(shape.parameters[i,2] + small_constant)))*bin.weights[i,1]
  }
  
  density = sum(density_matrix)
#--------------------------function output -------------------------------------
  
  return(list(Density_Estimate = density,
              beta_parameters = shape.parameters,
              Study_level_weights = bin.weights,
              bin_densities = density_matrix,
              negative_counts = pos.neg.neu))
}

########################## Ground truth at the study ###########################

sum(validated.data$outcome ==1)/nrow(validated.data)
sum(validated.data$outcome ==-1)/nrow(validated.data)
sum(validated.data$outcome ==0)/nrow(validated.data)

########################## function implementation #############################

out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)  
out$Density_Estimate

########################## Study level bootstrap estimates #####################

study.level.bootstraps = function(number.samples){

  samples = vector(length = number.samples)
  for(i in 1:number.samples){
  validated.data = validation.data(Bin.bounds, Size = 25)
  out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)  
  samples[i] = out$Density_Estimate
  }
  
  return(samples)
}

study.boots = study.level.bootstraps(100000)



#--------distribution of bootstrap samples and density estimates----------------

histogram.bootstrap.estimates = ggplot(data.frame(study.boots), aes(x = study.boots))+
  geom_histogram(binwidth = 0.025, fill = "maroon", colour = "white")+
  geom_vline(xintercept = mean(study.boots,na.rm = T), colour = "black", linewidth = 0.9)+
  annotate("text",
           x = out$Density_Estimate, 
           y = 50,             # adjust this based on your y-axis scale
           label = paste0("Density Estimate = ", round(out$Density_Estimate, 3)),
           angle = 90, vjust = 5, hjust = 0) +
  labs(
    title = "Histogram of density bootstrap estimates_African Pipit",
    x = "Density estimates",
    y = "Frequency") +
  theme_minimal()

ggsave("_output/Historgram of bootstrap estimates study_African Pipit.jpeg",
       plot = histogram.bootstrap.estimates, 
       width = 6, 
       height = 4,
       dpi = 300)



################################################################################
#                Strategy 1 strata level density estimation                    #
################################################################################
# ---------------------Distribution shifts Strategy 1---------------------------


#sum(P(+|b)*P_s(b))
Strat1_call_density = function(examples,
                                small_constant=0.1){

#-----------Converting data set to a usable format------------------------------
  #validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))
  #examples = validated.data
  quantile_bounds = Bin.bounds   # e.g. c("60%" = -0.93496400, "90%" = -0.09486028)
  #bounds come from the unvalidated data
  
  bounds = as.numeric(c(min(species_choice$study_level_data$logit_conf),
                        Bin.bounds,
                        max(species_choice$study_level_data$logit_conf)))
  
  #counts in each bin
  pos.neg.neu = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu) <- paste0("Bin", 1:(length(bounds) - 1))
  pos.neg.neu
  
  #allocating to pos.neg.neu
  for(i in 2:length(bounds)){
    subset_rows = subset(examples,
                         logit_conf >= bounds[i - 1] &
                           logit_conf <= bounds[i])
    
    pos.neg.neu["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu["Unknown",i-1] = sum(subset_rows$outcome==0)
  }
  
  #bin weights
  bin.weights = matrix(NA, nrow = ncol(pos.neg.neu), ncol = 1)
  rownames(bin.weights) = paste0("Bin", 1:(length(bounds) - 1))
  colnames(bin.weights) = "weights"
  bin.weights
   
  
  proportions.pos.neg.neu = (rowSums(pos.neg.neu)/sum(rowSums(pos.neg.neu)))
  matrix.proportions.row = pos.neg.neu/(rowSums(pos.neg.neu)+1e-12)
  
  for(i in 1:ncol(pos.neg.neu)){
    bin.weights[i,] = matrix.proportions.row[1,i]*proportions.pos.neg.neu[1] +
      matrix.proportions.row[2,i]*proportions.pos.neg.neu[2]
  }
  

#----------------strata level density estimation--------------------------------
  density_matrix = matrix(NA,nrow=nrow(bin.weights), ncol =1)
 
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
    ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*bin.weights[i,1]
  }
  
  density = sum(density_matrix)
  
#----------------------------------function output------------------------------
  return(list(Density_Estimate = density,
              bin_densities = density_matrix))
}

############################ Bootstrap data ####################################
strategy1.bootstraps = function(numberSamples){
  #numberSamples = 10
  bootstrap.estimates.32 = vector(length = numberSamples)
  bootstrap.estimates.48 = vector(length = numberSamples)
  
  #includes a function call for the sampling rate
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==32))
    bootstrap.estimates.32[i] = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==32),0.1)$Density_Estimate
  }
  
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))
    bootstrap.estimates.48[i] = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==48),0.1)$Density_Estimate
  }
  
  return(list(Khz32.estimates = as.numeric(bootstrap.estimates.32),
              Khz48.estimates = as.numeric(bootstrap.estimates.48)))
}


########################## function implementation #############################
#Implementation of strategy 1 call density  

validated.data = validation.data(Bin.bounds,Size = 25, species_choice$study_level_data)

out1 = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==32),0.1)
out1$Density_Estimate

out1.1 = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==48),0.1)
out1.1$Density_Estimate

strategy.bootstraps = strategy1.bootstraps(2000)
hist(as.numeric(strategy.bootstraps$Khz32.estimates))
hist(as.numeric(strategy.bootstraps$Khz48.estimates))

#---------------Implementation of study_level call density----------------------  

histogram.bootstrap.estimates.strata1= ggplot(data.frame(strategy.bootstraps$Khz32.estimates), aes(x = strategy.bootstraps$Khz32.estimates ))+
  geom_histogram(binwidth = 0.008, fill = "maroon", colour = "white")+
  geom_vline(xintercept = mean(strategy.bootstraps$Khz32.estimates), colour = "black", linewidth = 0.9)+
  labs(
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()


histogram.bootstrap.estimates.strata2= ggplot(data.frame(strategy.bootstraps$Khz48.estimates), aes(x = strategy.bootstraps$Khz48.estimates))+
  geom_histogram(binwidth = 0.018, fill = "maroon", colour = "white")+
  geom_vline(xintercept = mean(strategy.bootstraps$Khz48.estimates), colour = "black", linewidth = 0.9)+
  labs(
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()


 histogram.bootstrap.estimates.stratas = grid.arrange(histogram.bootstrap.estimates.strata1,
             histogram.bootstrap.estimates.strata2,
             top = textGrob("Strategy 1 bootstrap estimates_African Pipit", gp = gpar(fontsize = 16, fontface = "bold")),
             ncol = 2)

ggsave("_output/histogram bootstrap estimates stratas strategy1_African Pipit.jpeg", 
       plot =  histogram.bootstrap.estimates.stratas,
       width = 6,
       height = 4,
       dpi = 300)



###################### Ground truths at the strata level #######################
ground.truth.data.32 = filter(species_choice$study_level_data,species_choice$study_level_data$SamplingRate==32)
ground.truth.data.48 = filter(species_choice$study_level_data,species_choice$study_level_data$SamplingRate==48)

sum(ground.truth.data.32$outcome==1)/nrow(ground.truth.data.32)
sum(ground.truth.data.48$outcome==1)/nrow(ground.truth.data.48)



################################################################################
#                Strategy 2 strata level density estimation                    #
################################################################################
#---------Strategy 2------------------------------------------------------------
Strat2_call_density = function(site_examples, 
                               study_fit = out,
                               eps = 1e-12,
                               small_constant= 0.1) {
  

#---------------------------Binning data and format-----------------------------

  validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))
  quantile_bounds = Bin.bounds
  
  bin_indices = findInterval(site_examples$logit_conf ,
                             quantile_bounds,
                             rightmost.closed = TRUE)+1                         #returns the bin number for each score
  
  bin_weights = table(bin_indices) / nrow(site_examples)                      #gives the proportions in the different bins
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

bootstrap_strat2_site = function(site_examples, study_fit ,
                                 Bdraws = 2000,
                                 conf = 0.95,
                                 eps = 1e-12) {
  
  #site_examples=validation_data
  n = nrow(site_examples)
  point = Strat2_call_density(site_examples,
                              study_fit = study_fit,
                              eps = eps)$Stat2_density_estimate
  boots = Bdraws
  
  for (b in seq_len(Bdraws)) {
    idx = sample.int(n, n, replace = TRUE)
    boots[b] = Strat2_call_density(site_examples[idx,], study_fit = study_fit, eps = eps)$Stat2_density_estimate
  }
  
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(unlist(boots), probs = c(alpha, 1 - alpha)))
  
  
  return(list(Density_estimate = point, 
              se = sd(boots),
              conf_int = ci, 
              Samples = boots, 
              conf = conf))
}


########################## function implementation #############################

out2.density = Strat2_call_density(strata1,
                                   study_fit = out,
                                   eps = 1e-12)
out2.2.density = Strat2_call_density(strata2,
                                   study_fit = out,
                                   eps = 1e-12)

out2 = bootstrap_strat2_site(strata1,out,Bdraws = 2000,
                             conf = 0.95,
                             eps = 1e-12)
out2.2 = bootstrap_strat2_site(strata2,out,Bdraws = 2000,
                             conf = 0.95,
                             eps = 1e-12)

#---------implementation of strategy 1 density estimation (site)----------------

histogram.bootstrap.estimates.s2= ggplot(data.frame(out2$Samples), aes(x = out2.Samples))+
  geom_histogram(binwidth = 0.0014, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out2.density$Stat2_density_estimate, colour = "black", linewidth = 0.9)+
  labs(
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()

histogram.bootstrap.estimates.s2.2= ggplot(data.frame(out2.2$Samples), aes(x = out2.2.Samples))+
  geom_histogram(binwidth = 0.0014, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out2.2.density$Stat2_density_estimate, colour = "black", linewidth = 0.9)+
  labs(
    x = "Density estimates",
    y = "Frequency"
  ) +
  theme_minimal()


histogram.bootstrap.estimates.strategy2 = grid.arrange(histogram.bootstrap.estimates.s2,
                                                     histogram.bootstrap.estimates.s2.2,
                                                     top = textGrob("Strategy 2 bootstrap Distribution_African Pipit",gp = gpar(fontsize = 16, fontface = "bold")),
                                                     ncol = 2)


ggsave("_output/Historgram of bootstrap estimates strategy2_African Pipit.jpeg", 
       plot = histogram.bootstrap.estimates.strategy2, 
       width = 6,
       height = 4,
       dpi = 300)



#--------------------95% bootstrap CI, std and variance  -----------------------
confidence_interval = quantile(out2$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out2$Samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)


results = data.frame(
  Statistic = c("95% Confidence Interval (lower)",
                "95% Confidence Interval (upper)",
                "Variance",
                "Standard Deviation",
                "Density Estimate"),
  Value = c(round(confidence_interval[1], 4),
            round(confidence_interval[2], 4),
            round(Variance_bootstrap_samples, 4),
            round(Sd_bootstrap_samples, 4),
            round(out2.density$Stat2_density_estimate,4))
)


confidence_interval2.2 = quantile(out2.2$Samples, probs=c(0.025,0.975))
Variance_bootstrap_samples2.2 = var(out2.2$Samples)
Sd_bootstrap_samples2.2 = sqrt(Variance_bootstrap_samples2.2)


results2.2 = data.frame(
  Statistic = c("95% Confidence Interval (lower)",
                "95% Confidence Interval (upper)",
                "Variance",
                "Standard Deviation",
                "Density Estimate"),
  Value = c(round(confidence_interval1.1[1], 4),
            round(confidence_interval1.1[2], 4),
            round(Variance_bootstrap_samples1.1, 4),
            round(Sd_bootstrap_samples1.1, 4),
            round(out2.2.density$Stat2_density_estimate,4))
)


side_by_side <- data.frame(
  Statistic = results$Statistic,
  Stratum_1 = results$Value,
  Stratum_2 = results2.2$Value
)

side_by_side


# View as table
latex.table.strat1 = kable(side_by_side, format = "latex")
writeLines(latex.table.strat1, "_output/bootstrap_summary_strat2_African Pipit.tex")




#-------------------------------kl divergence plot------------------------------

kl.data.strata1 = cbind(out2.density$qs,out2.density$KL_values)
kl.plot.strata1 = ggplot(data = data.frame(kl.data.strata1), aes(x = X1,y = X2))+
  geom_line()+
  geom_vline(xintercept = out2.density$Stat2_density_estimate, colour = "red")+
  labs(title = "Plot of KL divergence strata 1_African Pipit",
  x = "q values",
  y = "KL divergence") +
  theme_minimal()

kl.data.strata2 = cbind(out2.2.density$qs,out2.2.density$KL_values)
kl.plot.strata2 = ggplot(data = data.frame(kl.data.strata2), aes(x = X1,y = X2))+
  geom_line()+
  geom_vline(xintercept = out2.2.density$Stat2_density_estimate, colour = "red")+
  labs(title = "Plot of KL divergence Strata 2_African Pipit",
       x = "q values",
       y = "KL divergence") +
  theme_minimal()


KL.strategy2 = grid.arrange(kl.plot.strata1,
                            kl.plot.strata2,
                            ncol = 2)

ggsave( "_output/KL divergences_African Pipit.jpeg", 
        plot = KL.strategy2, 
        width = 6,
        height = 4, 
        dpi = 300)


################################################################################
#                Strategy 3 strategy level density estimation                  #
################################################################################
#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(density1, density2){
  density_estimate = (density1 * density2)^(1/2)
  
  return (Density_estimate = density_estimate)
}



#-------------------Bootstrapping and confidence intervals----------------------
bootstrap_strat3_site = function(strategy1.samples,
                                 strategy2.samples,
                                 conf = 0.95) {
  
# -------------------------------- Bootstrap samples ---------------------------
  boots = (strategy1.samples*strategy2.samples)^(1/2)
  
#------------------------95% Confidence interval--------------------------------
  alpha = (1 - conf) / 2
  ci = as.numeric(quantile(boots, probs = c(alpha, 1 - alpha)))
 
  
#------------------------output of the function---------------------------------    
  return(list(samples = boots,
              se = sd(boots),
              conf_int = ci, 
              conf = conf))
}

########################## function implementation #############################

out3.strata1 = Strat3_call_density(Strat1_call_density(strata1,num_beta_samples=1000,
                                                       small_constant=0.1)$Density_Estimate,
                                   Strat2_call_density(strata1,study_fit = out,
                                                       eps = 1e-12, small_constant = 0.1)$Stat2_density_estimate)


out3.strata2 = Strat3_call_density(Strat1_call_density(strata2,num_beta_samples=1000,
                                                 small_constant=0.1)$Density_Estimate,
                                   Strat2_call_density(strata2,study_fit = out,
                                                       eps = 1e-12, small_constant = 0.1)$Stat2_density_estimate)


out3.B.strata1 = bootstrap_strat3_site(out1$Samples,out2$Samples,conf = 0.95)
out3.B.strata2 = bootstrap_strat3_site(out1.1$Samples,out2.2$Samples,conf = 0.95)


histogram.bootstrap.strategy3= ggplot(data.frame(out3.B.strata1$samples), aes(x = out3.B.strata1$samples))+
  geom_histogram(binwidth = 0.005, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out3.strata1, colour = "black",linewidth = 0.9)+
  labs(x = "Density estimates",
       y = "Frequency") +
  theme_minimal()

histogram.bootstrap.strategy3.1= ggplot(data.frame(out3.B.strata2$samples), aes(x = out3.B.strata2$samples))+
  geom_histogram(binwidth = 0.005, fill = "maroon", colour = "white")+
  geom_vline(xintercept = out3.strata2, colour = "black",linewidth = 0.9)+
  labs(x = "Density estimates",
       y = "Frequency") +
  theme_minimal()


Strategy3.distribution = grid.arrange(histogram.bootstrap.strategy3,
                            histogram.bootstrap.strategy3.1,
                            top = textGrob("Strategy 3 bootstrap distribution_African Pipit",
                                           gp = gpar(fontsize = 16, fontface = "bold")),
                            ncol = 2)

ggsave("_output/Historgram of bootstrap distribution strategy 3_African Pipit.jpeg",
       plot = Strategy3.distribution, 
       width = 6,
       height = 4,
       dpi = 300)
#------------------------------95% Confidence ----------------------------------

confidence_interval = quantile(out3.B.strata1$samples, probs=c(0.025,0.975))
Variance_bootstrap_samples = var(out3.B.strata1$samples)
Sd_bootstrap_samples = sqrt(Variance_bootstrap_samples)


results = data.frame(
  Statistic = c("95% Confidence Interval (lower)",
                "95% Confidence Interval (upper)",
                "Variance",
                "Standard Deviation",
                "Density Estimate"),
  Value = c(round(confidence_interval[1], 4),
            round(confidence_interval[2], 4),
            round(Variance_bootstrap_samples, 4),
            round(Sd_bootstrap_samples, 4),
            round(out3.strata1,4))
)


confidence_interval3.3 = quantile(out3.B.strata2$samples, probs=c(0.025,0.975))
Variance_bootstrap_samples3.3 = var(out3.B.strata2$samples)
Sd_bootstrap_samples3.3 = sqrt(Variance_bootstrap_samples3.3)


results3.3 = data.frame(
  Statistic = c("95% Confidence Interval (lower)",
                "95% Confidence Interval (upper)",
                "Variance",
                "Standard Deviation",
                "Density Estimate"),
  Value = c(round(confidence_interval1.1[1], 4),
            round(confidence_interval1.1[2], 4),
            round(Variance_bootstrap_samples1.1, 4),
            round(Sd_bootstrap_samples1.1, 4),
            round(out3.strata2,4))
)


side_by_side <- data.frame(
  Statistic = results$Statistic,
  Stratum_1 = results$Value,
  Stratum_2 = results3.3$Value
)

side_by_side


# View as table
latex.table.strat1 = kable(side_by_side, format = "latex")
writeLines(latex.table.strat1, "_output/bootstrap_summary_strat3_African Pipit.tex")


################################################################################
#                                Bin number play                               #
################################################################################
