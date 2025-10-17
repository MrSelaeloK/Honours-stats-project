
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


 
####################### Validated input data ###################################
strata1.v = read_csv("48KHz_validation.csv") 
strata2.v = read_csv("32KHz_validation.csv")

strata1.v = rename(strata1.v,`Common Name` = commonName)
strata2.v = rename(strata2.v,`Common Name` = commonName)

strata1.v$SamplingRate = 48
strata2.v$SamplingRate = 32


#######################  choosing species  #####################################
species = function(species_name){
  strata1 = filter(strata1, `Common Name` == species_name)
  strata2 = filter(strata2, `Common Name` == species_name)
  #---validated-------------------------------------------
  strata1.v = filter(strata1.v, `Common Name` == species_name)
  strata2.v = filter(strata2.v, `Common Name` == species_name)
  

 # x = sample(c(-1,0, 1), size = nrow(strata1),replace = TRUE, prob = c(0.1, 0.3, 0.6))
 # strata1$outcome = x
  
  
  
  study_level_data = rbind(strata1,strata2)
  study_level_data.v = rbind(strata1.v,strata2.v)
  study_level_data$logit_conf = log(study_level_data$score)
  strata1$logit_conf = log(strata1$score)
  strata2$logit_conf = log(strata2$score)
  
  
  
  return(list(study_level_data = study_level_data,
              strata1 = strata1,
              strata2 = strata2,
              study_level_data.v = study_level_data.v,
              strata1.v = strata1.v,
              strata2.v = strata2.v))
}

#Testing
#as.vector(distinct(dat ,dat$`Common Name` ))
species_choice = species("African Pipit")

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

############################# Quantiling function ##############################

Quantile_function = function(quantiles,specie){
  
  binning.values = c()                                  #initial empty binning dataset
  examples = species(specie)$study_level_data
  logits = sort(examples$logit_conf)
  
  
  for (i in 1:length(quantiles)){
    binning.values = c(binning.values,quantile(logits,quantiles[i]))
  }
  return(Quantiles = binning.values)
}

#Testing

Quantile_function(quantiles.3, "African Pipit")


################################################################################
#                   annotations  and bin weighs                                #
################################################################################
annotations = function(quantiles, specie){
  
  species.data = species(specie)
  Bin.bounds = Quantile_function(quantiles,specie)
  bounds = c(min(species.data$study_level_data$logit_conf),Bin.bounds, max(species.data$study_level_data$logit_conf))
  
  #annotations
  strata1.annotated = data.frame()
  strata2.annotated = data.frame()
  ground.truths = matrix(0,nrow = 1, ncol = 2)
  
  #bin weights
  strata1.bins.weights = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 1)
  strata2.bins.weights = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 1)
  study.level.bin.weights = matrix(NA,nrow = length(Bin.bounds)+1,ncol = 1)
  
  
  
  
  for (i in 2:length(bounds)){
    lower = bounds[i - 1]
    upper = bounds[i]
    
    #Strata 1 annotated data
    subset_rows.v = subset(species.data$strata1.v , logit_conf > lower & logit_conf <= upper)
    pi = nrow(filter(subset_rows.v,outcome == 1))/nrow(subset_rows.v)         #the proportion of positives
    
    subset_rows = subset(species.data$strata1, logit_conf > lower & logit_conf <= upper)
    x = sample(c(-1,0, 1), size = nrow(subset_rows),replace = TRUE, prob = c((1-pi), 0, pi))
    subset_rows$outcome = x
    strata1.annotated = rbind(strata1.annotated, subset_rows)
    
    strata1.bins.weights[i-1,1] = nrow(subset_rows)/nrow(species.data$strata1)  #bin weights
    
    #Strata 2 annotated data
    subset_rows.v = subset(species.data$strata2.v , logit_conf > lower & logit_conf <= upper)
    pi = nrow(filter(subset_rows.v,outcome == 1))/nrow(subset_rows.v)         #the proportion of positives
    
    subset_rows = subset(species.data$strata2, logit_conf > lower & logit_conf <= upper)
    x = sample(c(-1,0, 1), size = nrow(subset_rows),replace = TRUE, prob = c((1-pi), 0, pi))
    subset_rows$outcome = x
    strata2.annotated = rbind(strata2.annotated, subset_rows)
    
    strata2.bins.weights[i-1,1] = nrow(subset_rows)/nrow(species.data$strata2)  #bin weights
  }
  
  
  
  study.level.data.annotated = rbind(strata1.annotated, strata2.annotated)
  
  #study level bin weights
  for (i in 2:length(bounds)){
    lower = bounds[i - 1]
    upper = bounds[i]
    
    #Strata 1 annotated data
    subset_rows = subset(study.level.data.annotated , logit_conf > lower & logit_conf <= upper)
    study.level.bin.weights[i-1,1] = nrow(subset_rows)/nrow(study.level.data.annotated)
  }
  
  
  #ground truths
  strat1.ground = nrow(filter(strata1.annotated,outcome == 1))/nrow(strata1.annotated) 
  strat2.ground = nrow(filter(strata2.annotated,outcome == 1))/nrow(strata2.annotated) 
  study.ground = nrow(filter(study.level.data.annotated,outcome ==1))/nrow(study.level.data.annotated)
  
    return(list(strata1 = strata1.annotated,
                strata2 = strata2.annotated,
                study_level_data = study.level.data.annotated,
                strata1.ground = strat1.ground,
                strata2.ground = strat2.ground,
                study.ground = study.ground,
                strata1.bins.weights = strata1.bins.weights,
                strata2.bins.weights = strata2.bins.weights,
                study.level.bin.weights = study.level.bin.weights))
  
}

annotations(quantiles.3, "African Pipit")





################################################################################
#                            Validation dataset                                #
################################################################################
#-------------------------Getting the samples froms bins------------------------

validation.data = function(quantiles,specie,Size, dataset){

  Bin.bounds = Quantile_function(quantiles, specie)
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

validation.data(quantiles.3, "African Pipit", 25, annotate.data$study_level_data)


################################################################################
#                        Study level density estimation                        #
################################################################################

Study_level_call_density = function(quantiles,
                                    specie,
                                    Size){
  
#-------------------binning and weighting---------------------------------------
  Annotated.data = annotations(quantiles, specie)
  examples =  validation.data(quantiles,specie,Size,Annotated.data$study_level_data)
  small_constant = 0.1
  Bin.bounds = Quantile_function(quantiles, specie)   # as logit values
  species.data = species(specie)

  #bounds come from the unvalidated data
  bounds = as.numeric(c(min(examples$logit_conf),
             Bin.bounds,
             max(examples$logit_conf)))
  
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
  
  

#--------------------------Density  estimation----------------------------------
  
  # Create beta distributions for each bin.
  shape.parameters = matrix(NA,ncol=2,nrow=ncol(pos.neg.neu))
  
  for (i in 1:ncol(pos.neg.neu)){
    shape.parameters[i,] = c(pos.neg.neu[1,i] + small_constant,
                             pos.neg.neu[2,i] + small_constant)
  }
  
  #estimating the density
  density_matrix = matrix(NA,nrow = nrow(Annotated.data$study.level.bin.weights),ncol = 1)
  
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = ((shape.parameters[i,1]+ small_constant)/
      ((shape.parameters[i,1]+ small_constant)+(shape.parameters[i,2] + small_constant)))*Annotated.data$study.level.bin.weights[i,1]
  }
  
  density = sum(density_matrix)
#--------------------------function output -------------------------------------
  
  return(list(Density_Estimate = density,
              beta_parameters = shape.parameters,
              bin_densities = density_matrix,
              negative_counts = pos.neg.neu))
}


Study_level_call_density(quantiles.3,"African Pipit",25) 

########################## Study level bootstrap estimates #####################

study.level.bootstraps = function(number.samples, quantiles, specie,Size){

  samples = vector(length = number.samples)
  for(i in 1:number.samples){
  out = Study_level_call_density(quantiles, specie,Size)  
  samples[i] = out$Density_Estimate
  }
  
  return(samples)
}

study.level.bootstraps(200,quantiles.3, "African Pipit",25)

################################################################################
#                Strategy 1 strata level density estimation                    #
################################################################################

#sum(P(+|b)*P_s(b))
Strat1_call_density = function( quantiles, specie,Size){

  small_constant = 0.1
  species.data = species(specie)
  Bin.bounds = Quantile_function(quantiles, specie)   # as logit values
  Annotated.data = annotations(quantiles, specie)
  examples.strata1 =  validation.data(quantiles,specie,Size,Annotated.data$strata1)
  examples.strata2 =  validation.data(quantiles,specie,Size,Annotated.data$strata2)
  quantile_bounds = Quantile_function(quantiles, specie)
  out = Study_level_call_density(quantiles, specie,Size)

  bounds.strata1 = as.numeric(c(min(species.data$strata1$logit_conf),
                        Bin.bounds,
                        max(species.data$strata1$logit_conf)))
  
  bounds.strata2 = as.numeric(c(min(species.data$strata2$logit_conf),
                                Bin.bounds,
                                max(species.data$strata2$logit_conf)))
  
  
  #counts in each bin
  pos.neg.neu.strata1 = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu.strata1) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu.strata1) <- paste0("Bin", 1:(length(bounds.strata1) - 1))
  pos.neg.neu.strata1
  
  pos.neg.neu.strata2 = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu.strata2) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu.strata2) <- paste0("Bin", 1:(length(bounds.strata2) - 1))
  pos.neg.neu.strata2
  
  
  #allocating to pos.neg.neu.strata 1
  for(i in 2:length(bounds.strata1)){
    subset_rows = subset(examples.strata1,
                         logit_conf >= bounds.strata1[i - 1] &
                           logit_conf <= bounds.strata1[i])
    
    #sampling from the data
    pos.neg.neu.strata1["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu.strata1["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu.strata1["Unknown",i-1]  = sum(subset_rows$outcome==0)
  }
  
    #Strata2
  for(i in 2:length(bounds.strata2)){
    subset_rows = subset(examples.strata2,
                         logit_conf >= bounds.strata2[i - 1] &
                           logit_conf <= bounds.strata2[i])
    
    #sampling from the data
    pos.neg.neu.strata2["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu.strata2["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu.strata2["Unknown",i-1]  = sum(subset_rows$outcome==0)
  }
    

#----------------strata level density estimation--------------------------------
  density_matrix.strata1 = matrix(NA,nrow=nrow(Annotated.data$strata1.bins.weights), ncol =1)
  density_matrix.strata2 = matrix(NA,nrow=nrow(Annotated.data$strata2.bins.weights), ncol =1)

 
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix.strata1[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
    ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*Annotated.data$strata1.bins.weights[i,1]
    
    density_matrix.strata2[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
    ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*Annotated.data$strata2.bins.weights[i,1]
  }
  
  density.strata1 = sum(density_matrix.strata1)
  density.strata2 = sum(density_matrix.strata2)
  
#----------------------------------function output------------------------------
  return(list(Density_Estimate.strata1 = density.strata1,
              Density_Estimate.strata2 = density.strata2,
      
              bin_densities.strata1 = density_matrix.strata1,
              bin_densities.strata2 = density_matrix.strata2,
              
              allocations.strata1 = pos.neg.neu.strata1,
              allocations.strata2 = pos.neg.neu.strata2))
}

Strat1_call_density(quantiles.3,"African Pipit",25)
 
############################ repeated estimates ####################################
strategy1.bootstraps = function(numberSamples, quantiles, specie,Size){
  
  bootstrap.estimates.32 = vector(length = numberSamples)
  bootstrap.estimates.48 = vector(length = numberSamples)
  
  #includes a function call for the sampling rate
  for(i in 1:numberSamples){
    bootstrap.estimates.32[i] = Strat1_call_density(quantiles, specie, Size)$Density_Estimate.strata1
  }
  
  for(i in 1:numberSamples){
    bootstrap.estimates.48[i] = Strat1_call_density(quantiles, specie, Size)$Density_Estimate.strata2
  }
  
  return(list(Khz32.estimates = as.numeric(bootstrap.estimates.32),
              Khz48.estimates = as.numeric(bootstrap.estimates.48)))
}

 strategy1.bootstraps(20, quantiles.3, "African Pipit",25)

################################################################################
#                Strategy 2 strata level density estimation                    #
################################################################################
Strat2_call_density = function(quantiles, specie, Size) {
  
  
########################### Binning data and format ############################

  eps = 1e-12
  small_constant = 0.1
  Bin.bounds = Quantile_function(quantiles, specie)
  Annotated.data = annotations(quantiles, specie)
  examples.strata1 =  validation.data(quantiles,specie,Size,Annotated.data$strata1)
  examples.strata2 =  validation.data(quantiles,specie,Size,Annotated.data$strata2)
  out = Study_level_call_density(quantiles, specie,Size)
  species.data = species(specie)
  
  
  bounds.strata1 = as.numeric(c(min(species.data$strata1$logit_conf),
                                Bin.bounds,
                                max(species.data$strata1$logit_conf)))
  
  bounds.strata2 = as.numeric(c(min(species.data$strata2$logit_conf),
                                Bin.bounds,
                                max(species.data$strata2$logit_conf)))
  
  #counts in each bin
  pos.neg.neu.strata1 = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu.strata1) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu.strata1) <- paste0("Bin", 1:(length(bounds.strata1) - 1))
  pos.neg.neu.strata1
  
  pos.neg.neu.strata2 = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 3)
  rownames(pos.neg.neu.strata2) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu.strata2) <- paste0("Bin", 1:(length(bounds.strata2) - 1))
  pos.neg.neu.strata2
  
  #allocating to pos.neg.neu.strata 1
  for(i in 2:length(bounds.strata1)){
    subset_rows = subset(examples.strata1,
                         logit_conf >= bounds.strata1[i - 1] &
                           logit_conf <= bounds.strata1[i])
    
    #sampling from the data
    pos.neg.neu.strata1["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu.strata1["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu.strata1["Unknown",i-1]  = sum(subset_rows$outcome==0)
  }
    
  for(i in 2:length(bounds.strata2)){  
    #Strata2
    subset_rows = subset(examples.strata2,
                         logit_conf >= bounds.strata2[i - 1] &
                           logit_conf <= bounds.strata2[i])
    
    #sampling from the data
    pos.neg.neu.strata2["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu.strata2["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu.strata2["Unknown",i-1]  = sum(subset_rows$outcome==0)
  }
  
  
######################### site density estimation ##############################
  
  density_matrix.strata1 = matrix(NA,nrow=nrow(Annotated.data$strata1.bins.weights), ncol =1)
  density_matrix.strata2 = matrix(NA,nrow=nrow(Annotated.data$strata2.bins.weights), ncol =1)
  
  for(i in 1:nrow(out$beta_parameters)){
    density_matrix.strata1[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
        ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*Annotated.data$strata1.bins.weights[i,1]
    
    density_matrix.strata2[i,1] =  ((out$beta_parameters[i,1]+small_constant)/
        ((out$beta_parameters[i,1]+small_constant)+(out$beta_parameters[i,2]+small_constant)))*Annotated.data$strata2.bins.weights[i,1]
  }
  
  site_level_density.strata1 = sum(density_matrix.strata1)
  site_level_density.strata2 = sum(density_matrix.strata2)
  
  
 #-------- Study level: expected P(b|+) and P(b|-) -----------------------------
  alpha = out$beta_parameters[,1]                               # = k_pos + c
  beta = out$beta_parameters[,2]                                # = k_neg + c
  E_pos_given_b = alpha / (alpha + beta)                        # E[P(+|b)]
  P_b = as.vector(Annotated.data$study.level.bin.weights)       # study P(b)
  P_b = pmax(P_b, eps)
  P_b = P_b / sum(P_b)                          # Also to make sure in case < eps came true
  
  P_pos = sum(E_pos_given_b * P_b)
  P_b_pos = (E_pos_given_b * P_b) / P_pos              # P(+|b)P(b)/P(+)
  P_b_neg = ((1 - E_pos_given_b) * P_b) / (1 - P_pos)  # P(-|b)P(b)/(1-P(+))
  
  P_b_pos = pmax(P_b_pos, eps)
  P_b_pos = P_b_pos / sum(P_b_pos)                     # in case < eps was true
  P_b_neg = pmax(P_b_neg, eps)
  P_b_neg = P_b_neg / sum(P_b_neg) 

  
######################### grid iteration/search ################################
  
  Ps_b.strata1 = P_b_pos*site_level_density.strata1 + P_b_neg*(1-site_level_density.strata1)
  Ps_b.strata2 = P_b_pos*site_level_density.strata2 + P_b_neg*(1-site_level_density.strata2)
    
  q_seq = seq(0.001, 0.999, length.out = 2001)                # avoid exact 0/1
  KL_vals.strata1 = sapply(q_seq, function(q) {
    Q = q * P_b_pos + (1 - q) * P_b_neg
    Q = pmax(Q, eps)
    Q = Q / sum(Q)                                            # in case < eps was true
    sum(Ps_b.strata1 * (log(Ps_b.strata1) - log(Q)))
  })
  
  KL_vals.strata2 = sapply(q_seq, function(q) {
    Q = q * P_b_pos + (1 - q) * P_b_neg
    Q = pmax(Q, eps)
    Q = Q / sum(Q)                                            # in case < eps was true
    sum(Ps_b.strata2 * (log(Ps_b.strata2) - log(Q)))
  })
  
  
  strata1.density.estimate = q_seq[which.min(KL_vals.strata1)]
  strata2.density.estimate = q_seq[which.min(KL_vals.strata2)]
  
  
  return(list(Density_estimate.strata1 = strata1.density.estimate,
              Density_estimate.strata2 = strata2.density.estimate,
              KL_values.strata1 = KL_vals.strata1,
              KL_values.strata1 = KL_vals.strata1,
              qs = q_seq))
}


 Strat2_call_density(quantiles.3,"African Pipit", 25) 
 
############################ repeat estimates #####################################
 
strategy2.bootstraps = function(numberSamples, quantiles, specie, Size){
  bootstrap.estimates.32 = vector(length = numberSamples)
  bootstrap.estimates.48 = vector(length = numberSamples)
  
  #includes a function call for the sampling rate
  for(i in 1:numberSamples){
    strategy2.output = Strat2_call_density(quantiles,specie,Size)
    bootstrap.estimates.32[i] = strategy2.output$Density_estimate.strata2
    bootstrap.estimates.48[i] = strategy2.output$Density_estimate.strata1
  
  }
  
  return(list(Khz32.estimates = as.numeric(bootstrap.estimates.32),
              Khz48.estimates = as.numeric(bootstrap.estimates.48)))
  
}


strategy2.bootstraps(20,quantiles.3, "Baglafecht Weaver",25)
 
################################################################################
#                Strategy 3 strategy level density estimation                  #
################################################################################
#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(strategy1.re_estimates, strategy2.re_estimates){
  
  strategy1.32.density = mean(strategy1.re_estimates$Khz32.estimates)
  strategy2.32.density = mean(strategy2.re_estimates$Khz32.estimates)
  
  strategy1.48.density = mean(strategy1.re_estimates$Khz48.estimates)
  strategy2.48.density = mean(strategy2.re_estimates$Khz48.estimates)
  
  density_estimate.32 = (strategy1.32.density * strategy2.32.density)^(1/2)
  density_estimate.48 = (strategy1.48.density * strategy2.48.density)^(1/2)
  
  return (list(density.32Khz = density_estimate.32,
               density.48Khz = density_estimate.48))
}

out1 = strategy1.bootstraps(20, quantiles.3, "African Pipit",25)
out2 = strategy2.bootstraps(20, quantiles.3, "African Pipit",25)
Strat3_call_density(out1,out2)

################################################################################
#                Estimates and tables                                          #
################################################################################
Estimates = function(){
  list.of.species = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                      "African Gray Flycatcher", "Red-eyed Dove")
  
  estimates.groundtruths  = matrix(NA, nrow = 2, ncol = 5)
  rownames(estimates.groundtruths) = c("Estimate","Ground_Truth")
  colnames(estimates.groundtruths) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                       "African Gray Flycatcher", "Red-eyed Dove")
  
  estimates.groundtruths.1  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.1) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.1) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-eyed Dove")
  
  estimates.groundtruths.2  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.2) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.2) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-eyed Dove")
  
  estimates.groundtruths.3  = matrix(NA, nrow = 5, ncol = 4)
  colnames(estimates.groundtruths.3) = c("Estimate-32","Ground_Truth-32","Estimate-48","Ground_Truth-48")
  rownames(estimates.groundtruths.3) = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                                         "African Gray Flycatcher", "Red-eyed Dove")
  
  #variables
  quantiles = quantiles.3
  Size = 25
  numberSamples = 5
  
  
  for( i in 1:ncol(estimates.groundtruths)){
    
    #Study_level--------------------------
    
    estimates.groundtruths[1,i] = Study_level_call_density(quantiles, list.of.species[i], Size)$Density_Estimate
    estimates.groundtruths[2,i] = annotations(quantiles, list.of.species[i])$study.ground 
    
    #Strategy 1 --------------------------
    out1 = strategy1.bootstraps(numberSamples, quantiles, list.of.species[i],Size)
    estimates.groundtruths.1[i,1] = mean(out1$Khz32.estimates) 
    estimates.groundtruths.1[i,2] = annotations(quantiles, list.of.species[i])$strata1.ground
    
    out1.1 = strategy1.bootstraps(numberSamples, quantiles, list.of.species[i],Size)
    estimates.groundtruths.1[i,3] = mean(out1.1$Khz48.estimates)  
    estimates.groundtruths.1[i,4] = annotations(quantiles, list.of.species[i])$strata2.ground
    
    #Strategy 2 ---------------------------
    out2 = strategy2.bootstraps(numberSamples, quantiles, list.of.species[i],Size)
    estimates.groundtruths.2[i,1] = mean(out2$Khz32.estimates)  
    estimates.groundtruths.2[i,2] = annotations(quantiles, list.of.species[i])$strata1.ground
    
    out2.2 = strategy2.bootstraps(numberSamples, quantiles, list.of.species[i],Size)
    estimates.groundtruths.2[i,3] = mean(out2.2$Khz48.estimates) 
    estimates.groundtruths.2[i,4] = annotations(quantiles, list.of.species[i])$strata2.ground
    
    #Strategy 3 --------------------------
    out3 = Strat3_call_density(out1,out2)
    estimates.groundtruths.3[i,1] = out3$density.32Khz  
    estimates.groundtruths.3[i,2] = annotations(quantiles, list.of.species[i])$strata1.ground
    
    out3.3 = out3 = Strat3_call_density(out1.1,out2.2)
    estimates.groundtruths.3[i,3] = out3.3$density.48Khz  
    estimates.groundtruths.3[i,4] = annotations(quantiles, list.of.species[i])$strata2.ground
  }
    
  
  #View as table
  estimates.groundtruths.table = kable(estimates.groundtruths)
  estimates.groundtruths.table.1 = kable(estimates.groundtruths.1)
  estimates.groundtruths.table.2 = kable(estimates.groundtruths.2)
  estimates.groundtruths.table.3 = kable(estimates.groundtruths.3)
  
  
  return(list(estimates.groundtruths.table,
              estimates.groundtruths.table.1,
              estimates.groundtruths.table.2,
              estimates.groundtruths.table.3))
  
  
}  

Study.level.ground = kable(Estimates(), format = "latex")
writeLines(Study.level.ground, "_output/study.level.3bins.25.samples.tex")


debug(Estimates)
Estimates()

################################################################################
#                                Bin number play                               #
################################################################################
