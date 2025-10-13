
library(readr)
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(grid)

################################################################################
#                                 The datasets                                 #
################################################################################

#-----------------------------The datasets loading------------------------------
# ---- SETTINGS ----
#folder_path = "32KHz_raw_data"   # 🔹 Change this
#output_file = "combinedhz32Khz.csv"

# -------- STEP 1: Get all .txt files -------
#files = list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)

# --- STEP 2: Read each file safely ---
#read_file_safe <- function(f) {
#  message("Reading file: ", f)
#  df = read_delim(f, delim = "\t", col_types = cols(.default = "c"))
#  df$source_file = basename(f)
#  return(df)
#}

#data_all = map(files, read_file_safe)

# --- STEP 3: Combine all files ----
#combined = bind_rows(data_all)

#converting the columns to the right data
#combined$Confidence = as.double(combined$Confidence)
#write_csv(combined, output_file)


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

#species = as.vector(distinct(study_level_data,study_level_data$`Common Name`))
species_choice = species("African Gray Flycatcher")


#--------------- Some information on the dataset--------------------------------
summary(species_choice$study_level_data$score)
nrow(species_choice$study_level_data)
nrow(filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))

#--------visualization of the data set -------------------------------------
Confidence.histogram = ggplot(species_choice$study_level_data, aes(x = score))+
  geom_histogram(binwidth = 0.065, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of all Confidence Scores - Red-and-yellow Barbet",
    x = "Confidence",
    y = "Frequency"
  ) +
  theme_minimal()


ggsave("_output/confidence_histogram_Red-and-yellow Barbet.jpeg",
       plot = Confidence.histogram, 
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

################################################################################
#                            Validation dataset                                #
################################################################################
#-------------------------Getting the samples froms bins------------------------

validation.data = function(Bin.bounds,Size, dataset){

  #validation of the dataset
  x = sample(c(-1,0, 1), size = nrow(dataset),replace = TRUE, prob = c(0.1, 0.3, 0.6))
  dataset$outcome = x
  
  bounds = c(min(dataset$logit_conf),Bin.bounds, max(dataset$logit_conf))

  Validated_data = data.frame()

  for (i in 2:length(bounds)) {
    lower = bounds[i - 1]
    upper = bounds[i]
    
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

#validated.data = validation.data(Bin.bounds,Size = 25,species_choice$study_level_data)


################################################################################
#                        Study level density estimation                        #
################################################################################

Study_level_call_density = function(examples,
                                    num_beta_samples = 1000,
                                    small_constant = 0.1,
                                    quantile_bounds){
  
#-------------------binning and weighting---------------------------------------
  #examples = validated.data
  quantile_bounds = Bin.bounds   # as logit values
  #bounds come from the unvalidated data
  bounds = as.numeric(c(min(species_choice$study_level_data$logit_conf),
             Bin.bounds,
             max(species_choice$study_level_data$logit_conf)))
  
  #counts in each bin
  pos.neg.neu = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu) = paste0("Bin", 1:(length(bounds) - 1))
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

########################## Study level bootstrap estimates #####################

study.level.bootstraps = function(number.samples, out, validated.data){

  samples = vector(length = number.samples)
  for(i in 1:number.samples){
  validated.data = validation.data(Bin.bounds, Size = 25,species_choice$study_level_data)
  out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)  
  samples[i] = out$Density_Estimate
  }
  
  return(samples)
}


################################################################################
#                Strategy 1 strata level density estimation                    #
################################################################################

#sum(P(+|b)*P_s(b))
Strat1_call_density = function(examples,small_constant=0.1,out){

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
  
  bootstrap.estimates.32 = vector(length = numberSamples)
  bootstrap.estimates.48 = vector(length = numberSamples)
  
  #includes a function call for the sampling rate
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==32))
    bootstrap.estimates.32[i] = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==32),0.1,
                                                    Study_level_call_density(validated.data,1000,0.1,Bin.bounds))$Density_Estimate
  }
  
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))
    bootstrap.estimates.48[i] = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==48),0.1,
                                                    Study_level_call_density(validated.data,1000,0.1,Bin.bounds))$Density_Estimate
  }
  
  return(list(Khz32.estimates = as.numeric(bootstrap.estimates.32),
              Khz48.estimates = as.numeric(bootstrap.estimates.48)))
}



################################################################################
#                Strategy 2 strata level density estimation                    #
################################################################################
Strat2_call_density = function(examples, 
                               out,
                               eps = 1e-12,
                               small_constant= 0.1) {
  
  
########################### Binning data and format ############################

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
  
  
######################### site density estimation ##############################
  
  out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)
  density_matrix = matrix(NA,nrow=nrow(bin.weights), ncol =1)
  
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
                              ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*bin.weights[i,1]
  }
  
  site_level_density = sum(density_matrix)
  
  
  #Study level: expected P(b|+) and P(b|-)
  alpha = out$beta_parameters[,1]             # = k_pos + c
  beta = out$beta_parameters[,2]              # = k_neg + c
  E_pos_given_b = alpha / (alpha + beta)            # E[P(+|b)]
  P_b = as.vector(out$Study_level_weights)    # study P(b)
  P_b = pmax(P_b, eps)
  P_b = P_b / sum(P_b)                              # Also to make sure in case < eps came true
  
  P_pos = sum(E_pos_given_b * P_b)
  P_b_pos = (E_pos_given_b * P_b) / P_pos              # P(+|b)P(b)/P(+)
  P_b_neg = ((1 - E_pos_given_b) * P_b) / (1 - P_pos)  # P(-|b)P(b)/(1-P(+))
  
  P_b_pos = pmax(P_b_pos, eps)
  P_b_pos = P_b_pos / sum(P_b_pos)                     # in case < eps was true
  P_b_neg = pmax(P_b_neg, eps)
  P_b_neg = P_b_neg / sum(P_b_neg) 

  
######################### grid iteration/search ################################
  
  Ps_b = P_b_pos*site_level_density + P_b_neg*(1-site_level_density)
    
  q_seq = seq(0.001, 0.999, length.out = 2001)                # avoid exact 0/1
  KL_vals = sapply(q_seq, function(q) {
    Q = q * P_b_pos + (1 - q) * P_b_neg
    Q = pmax(Q, eps)
    Q = Q / sum(Q)                                            # in case < eps was true
    sum(Ps_b * (log(Ps_b) - log(Q)))
  })
  
  q_seq[which.min(KL_vals)]
  
  return(list(Stat2_density_estimate = q_seq[which.min(KL_vals)],
         KL_values = KL_vals,
         qs = q_seq))
}

############################ Bootstrapping #####################################
 
strategy2.bootstraps = function(numberSamples){
  bootstrap.estimates.32 = vector(length = numberSamples)
  bootstrap.estimates.48 = vector(length = numberSamples)
  
  #includes a function call for the sampling rate
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==32))
    bootstrap.estimates.32[i] = Strat2_call_density(filter(validated.data,validated.data$SamplingRate==32),
                                                    Study_level_call_density(validated.data,1000,0.1,Bin.bounds),
                                                    eps = 1e-12)$Stat2_density_estimate
  }
  
  for(i in 1:numberSamples){
    validated.data = validation.data(Bin.bounds,Size = 25,filter(species_choice$study_level_data, species_choice$study_level_data$SamplingRate==48))
    out.2 = Strat2_call_density(filter(validated.data,validated.data$SamplingRate==48),
                                Study_level_call_density(validated.data,1000,0.1,Bin.bounds),
                                eps = 1e-12)
    bootstrap.estimates.48[i] = out.2$Stat2_density_estimate
  }
  
  return(list(Khz32.estimates = as.numeric(bootstrap.estimates.32),
              Khz48.estimates = as.numeric(bootstrap.estimates.48)))
  
}

########################## function implementation #############################
strategy2.estimates = function(){
  list.of.species = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                      "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  estimates.groundtruths.2  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.2) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.2) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  for(i in 1:nrow(estimates.groundtruths.2)){
    species_choice = species(list.of.species[i])
    Bin.bounds = Quantile_function(species_choice$study_level_data,quantiles.3) #number of bins
    
    #x = sample(c(-1,0, 1),size = nrow(species_choice$study_level_data),replace = TRUE, prob = c(0.1, 0.3, 0.6)) #setting of validation
    #species_choice$study_level_data$outcome = x
    
    #validated.data = validation.data(Bin.bounds,Size = 25, species_choice$study_level_data) #validation of the bins
    #out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)  
    simulated.estimates = strategy2.bootstraps(numberSamples = 200)
    estimates.groundtruths.2[i,1] = mean(simulated.estimates$Khz32.estimates) 
    estimates.groundtruths.2[i,3] = mean(simulated.estimates$Khz48.estimates) 
    
    #ground truth
    Samp.32 = nrow(filter((filter(validated.data,validated.data$SamplingRate==32)),outcome == 1))
    Samp.48 = nrow(filter((filter(validated.data,validated.data$SamplingRate==48)),outcome == 1))
    
    ground.truth.32 = Samp.32/nrow((filter(validated.data,validated.data$SamplingRate==32)))
    ground.truth.48 = Samp.48/nrow((filter(validated.data,validated.data$SamplingRate==48)))
    estimates.groundtruths.2[i,2] = ground.truth.32
    estimates.groundtruths.2[i,4] = ground.truth.48
  }
  
  return(estimates.groundtruths.2)
}
strategy2.estimates()









density.estimates.strategy2 =  strategy2.bootstraps(numberSamples = 2000)
hist(density.estimates.strategy2$Khz32.estimates)
hist(density.estimates.strategy2$Khz48.estimates)




########## implementation of strategy 2 density estimation (site) ##############

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





################################################################################
#                Strategy 3 strategy level density estimation                  #
################################################################################
#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(density1, density2){
  strat1.32.density = mean(density.estimates.strategy2$Khz32.estimates)
  strat2.32.density = mean(strategy.bootstraps$Khz32.estimates)
  
  strat1.48.density = mean(density.estimates.strategy2$Khz48.estimates)
  strat2.48.density = mean(strategy.bootstraps$Khz32.estimates)
  
  density_estimate.32 = (strat1.32.density * strat2.32.density)^(1/2)
  density_estimate.48 = (strat1.48.density * strat2.48.density)^(1/2)
  
  return (list(density.32 = density_estimate.32,
               density.32 = density_estimate.48))
}



################################################################################
#                Estimates and tables                                          #
################################################################################
Estimates = function(){
  list.of.species = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                      "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  estimates.groundtruths  = matrix(NA, nrow = 2, ncol = 5)
  rownames(estimates.groundtruths) = c("Estimate","Ground_Truth")
  colnames(estimates.groundtruths) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                       "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  estimates.groundtruths.1  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.1) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.1) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  estimates.groundtruths.2  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.2) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.2) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-and-yellow Barbet")
  
  
  
  for( i in 1:ncol(estimates.groundtruths)){
    species_choice = species(list.of.species[i])
    Bin.bounds = Quantile_function(species_choice$study_level_data,quantiles.3) #number of bins
    validated.data = validation.data(Bin.bounds,Size = 25, species_choice$study_level_data) #validation of the bins
    
    #Study_level--------------------------
    out = Study_level_call_density(validated.data,1000,0.1,Bin.bounds)  
    repeat.estimates =  study.level.bootstraps(20,out,validated.data)
    estimates.groundtruths[1,i]= mean(repeat.estimates)
    
    #Strategy 1---------------------------
    out1 = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==32),0.1,
                               Study_level_call_density(validated.data,1000,0.1,Bin.bounds))
    repeat.estimates1.0 = strategy1.bootstraps(200)
    estimates.groundtruths.1[i,1]= mean(repeat.estimates1.0$Khz32.estimates)
    out1 = Strat1_call_density(filter(validated.data,validated.data$SamplingRate==48),0.1,
                               Study_level_call_density(validated.data,1000,0.1,Bin.bounds))
    repeat.estimates1.1 = strategy1.bootstraps(200)
    estimates.groundtruths.1[i,3]= mean(repeat.estimates1.1$Khz48.estimates)
    
    #ground truth
    Samp.32 = nrow(filter((filter(validated.data,validated.data$SamplingRate==32)),outcome == 1))
    Samp.48 = nrow(filter((filter(validated.data,validated.data$SamplingRate==48)),outcome == 1))
    
    ground.truth.32 = Samp.32/nrow((filter(validated.data,validated.data$SamplingRate==32)))
    ground.truth.48 = Samp.48/nrow((filter(validated.data,validated.data$SamplingRate==48)))
    estimates.groundtruths.1[i,2] = ground.truth.32
    estimates.groundtruths.1[i,4] = ground.truth.48
    
    #Strategy 2---------------------------
    simulated.estimates = strategy2.bootstraps(numberSamples = 200)
    estimates.groundtruths.2[i,1] = mean(simulated.estimates$Khz32.estimates) 
    estimates.groundtruths.2[i,3] = mean(simulated.estimates$Khz48.estimates) 
    
    #ground truth
    Samp.32 = nrow(filter((filter(validated.data,validated.data$SamplingRate==32)),outcome == 1))
    Samp.48 = nrow(filter((filter(validated.data,validated.data$SamplingRate==48)),outcome == 1))
    
    ground.truth.32 = Samp.32/nrow((filter(validated.data,validated.data$SamplingRate==32)))
    ground.truth.48 = Samp.48/nrow((filter(validated.data,validated.data$SamplingRate==48)))
    estimates.groundtruths.2[i,2] = ground.truth.32
    estimates.groundtruths.2[i,4] = ground.truth.48
    
    
  }
  #ground truths
  estimates.groundtruths[2,] = rep(0.6, ncol(estimates.groundtruths))
  
  #View as table
  estimates.groundtruths.table = kable(estimates.groundtruths)
  estimates.groundtruths.table.1 = kable(estimates.groundtruths.1)
  estimates.groundtruths.table.2 = kable(estimates.groundtruths.2)
  
  return(list(estimates.groundtruths.table,
              estimates.groundtruths.table.1,
              estimates.groundtruths.table.2))
}

Study.level.ground = kable(Estimates(), format = "latex")
writeLines(Study.level.ground, "_output/study.level.3bins.25.samples.tex")


debug(Estimates)
Estimates()

################################################################################
#                                Bin number play                               #
################################################################################
