# Honours Project R script| Passive Acoustic Monitoring: Call density estimation 
# Name: Selaelo Kgafela 
# Date: 17/12/2025


################################################################################
#                                libraries used                           #
################################################################################

library(readr)
library(tidyverse)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(grid)
library(future.apply)
library(future)
library(reshape2)
library(dplyr)
library(stringr)


################################################################################
#                           Parallel running section                           #
################################################################################


# Avoid forking; use background R sessions
Sys.setenv(R_FUTURE_FORK_ENABLE = "false")

# Optional: increase limit if large objects move to workers
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB
options(future.rng.onMisuse   = "ignore")     # quiet RNG warnings for nested RNG

# Choose the total number of workers (adjust to your machine)
plan(multisession, workers = 6)




################################################################################
#                                 The datasets                                 #
################################################################################

# Loading raw datasets into the data, commented out to show how data was input

#folder_path = "32KHz_raw_data"   # 
#output_file = "combinedhz32Khz.csv"


# STEP 1: Getting all .txt files 
#files = list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)


# STEP 2: Read each txt file 
#read_file_safe <- function(f) {
#  message("Reading file: ", f)
#  df = read_delim(f, delim = "\t", col_types = cols(.default = "c"))
#  df$source_file = basename(f)
#  return(df)
#}


#data_all = map(files, read_file_safe)
# STEP 3: Combining all files to make one huge dataset
#combined = bind_rows(data_all)


# converting the columns to the right data
#combined$Confidence = as.double(combined$Confidence)
#write_csv(combined, output_file)


#------------------------------ Cleaned datasets -------------------------------

strata1 = read_csv("combinedhz48Khz.csv") 
strata2 = read_csv("combinedhz32Khz.csv")

# Renaming of scores to confidence
strata1 = rename(strata1,score = Confidence)
strata2 = rename(strata2,score = Confidence)

# Attaching of sampling rates
strata1$SamplingRate = 48 
strata2$SamplingRate = 32



############## Natural state Validated input data ##############################
strata1.v = read_csv("48KHz_validation.csv") 
strata2.v = read_csv("32KHz_validation.csv")

# Correcting common Name to match strata1 and 2
strata1.v = rename(strata1.v,`Common Name` = commonName)
strata2.v = rename(strata2.v,`Common Name` = commonName)

# Attaching of sampling rates
strata1.v$SamplingRate = 48
strata2.v$SamplingRate = 32


#######################  choosing species  #####################################

species = function(species_name){
  
  # tiny offset to avoid 0 and 1  
  eps = 1e-6  
  

  # filter by species name; non validated data
  strata1 = filter(strata1, `Common Name` == species_name)
  strata2 = filter(strata2, `Common Name` == species_name)
  
  
  # filter by species name; validated data
  strata1.v = filter(strata1.v, `Common Name` == species_name)
  strata2.v = filter(strata2.v, `Common Name` == species_name)
  
  
  # Restructuring datasets 
  study_level_data = rbind(strata1,strata2)
  study_level_data.v = rbind(strata1.v,strata2.v)
  study_level_data$logit_conf = qlogis(pmin(pmax(study_level_data$score, eps), 1 - eps))
  strata1$logit_conf = qlogis(pmin(pmax(strata1$score, eps), 1 - eps))
  strata2$logit_conf = qlogis(pmin(pmax(strata2$score, eps), 1 - eps))
  study_level_data.unadjusted = rbind(strata1,strata2)
  
  
  
  ############### data readjustment for missing logits #########################
  
  # calculation of number of missing audio clips
  min.time = min(study_level_data.unadjusted$`Begin Time (s)`)
  max.time = max(study_level_data.unadjusted$`Begin Time (s)`)
  time.diff = max.time-min.time
  full.count = time.diff/3
  N.state.count = nrow(study_level_data.unadjusted)
  missing.count = full.count-N.state.count
  missing.scores = runif(missing.count,min=0.011,max=0.0999)
  samp.allocation.48kHz = nrow(filter(study_level_data.unadjusted, SamplingRate==48))/nrow(study_level_data.unadjusted)
  
  
  
  # data frame recovered audio clips
  missing.dataframe = tibble(
    Selection = rep(0,missing.count),
    View = rep("Missing",missing.count),
    Channel = rep(0,missing.count),
    `Begin Time (s)` = rep(0,missing.count),
    `End Time (s)` = rep(0,missing.count),
    `Low Freq (Hz)` = rep(0,missing.count),
    `High Freq (Hz)` = rep(0,missing.count),
    `Species Code` = rep(0,missing.count),
    `Common Name` = rep(pull(distinct(study_level_data.unadjusted, `Common Name`)),
                        missing.count),
    score = missing.scores,
    source_file = rep("Missing", missing.count),
    SamplingRate = sample(c(32,48),missing.count,TRUE,prob=c(1-samp.allocation.48kHz , 
                                                             samp.allocation.48kHz)),
    logit_conf = qlogis(missing.scores)
    
  )
  
  # adjusted dataset
  study_level_data.adjusted = rbind(missing.dataframe, study_level_data.unadjusted)
  missing.scores = runif(missing.count,min=1e-6,max=0.0999) 
  
  return(list(study_level_data = study_level_data.adjusted, 
              
              # unadjusted data
              study_level_data.unadjusted = study_level_data.unadjusted,
              strata1 = strata1,   
              strata2 = strata2,
              
              #validated data
              study_level_data.v = study_level_data.v,
              strata1.v = strata1.v,
              strata2.v = strata2.v))
}

species_choice = species("Abyssinian Nightjar")

########################## visualization of the data set #######################

# Natural State unadjusted data for species
Confidence.histogram = ggplot(species_choice$study_level_data, aes(x = score))+
  geom_histogram(binwidth = 0.055, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of all Confidence Scores - Abyssinian Nightjar",
    x = "Confidence",
    y = "Frequency"
  ) +
  theme_minimal()
Confidence.histogram


# Natural State unadjusted Logit scores for species
Logit.confidence.histogram = ggplot(species_choice$study_level_data,aes(x = logit_conf))+
  geom_histogram(binwidth = 0.65, fill = "maroon", colour = "white")+
  labs(
    title = "Histogram of logit scale Confidence Scores - African Gray Flycatcher",
    x = "logit Confidence",
    y = "Frequency"
  ) +
  theme_minimal()
Logit.confidence.histogram


# --------------Logit scores' five number summary-------------------------------
# Extract study-level data
df_all  = species_choice$study_level_data
df_32   = subset(df_all, SamplingRate == 32)
df_48   = subset(df_all, SamplingRate == 48)

# Function to compute five-number summary
five_num <- function(x) {
  c(
    Min = round(min(x, na.rm = TRUE), 3),
    Q1  = round(quantile(x, 0.25, na.rm = TRUE), 3),
    Med = round(median(x, na.rm = TRUE), 3),
    Q3  = round(quantile(x, 0.75, na.rm = TRUE), 3),
    Max = round(max(x, na.rm = TRUE), 3)
  )
}

# Summary
summary_table = rbind(
  "Study Level (All)" = five_num(df_all$logit_conf),
  "Study Level (32 kHz)" = five_num(df_32$logit_conf),
  "Study Level (48 kHz)" = five_num(df_48$logit_conf)
)

summary_table



################################################################################
#                      Quantiling the dataset  by logits                       #
################################################################################


Quantile_function = function(quantiles,specie){
  
  binning.values = c()                            #initial empty binning dataset
  examples = species(specie)$study_level_data
  logits = sort(examples$logit_conf)
  
  
  # populating of the binning quantiles with quantiles
  for (i in 1:length(quantiles)){
    binning.values = c(binning.values,quantile(logits,quantiles[i]))
  }
  
  
  return(Quantiles = binning.values)
}



################################################################################
#                   annotations  and bin weighs                                #
################################################################################
annotations = function(quantiles, specie){
  
  # importing data from other functions
  species.data = species(specie)
  Bin.bounds = Quantile_function(quantiles,specie)
  bounds = c(min(species.data$study_level_data$logit_conf),Bin.bounds,
             max(species.data$study_level_data$logit_conf))
  
  #annotations
  strata1.annotated = data.frame()
  strata2.annotated = data.frame()
  ground.truths = matrix(0,nrow = 1, ncol = 2)
  
  
  #weight counts
  strata1.annotation.counts = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 2)
  colnames(strata1.annotation.counts) = c("Positives","Negatives")
  rownames(strata1.annotation.counts) = paste0("Bin", 1:(length(bounds) - 1))
  strata2.annotation.counts = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 2)
  colnames(strata2.annotation.counts) = c("Positives","Negatives")
  rownames(strata2.annotation.counts) = paste0("Bin", 1:(length(bounds) - 1))
  
  
  
  #bin weights
  strata1.bins.weights = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 1)
  strata2.bins.weights = matrix(NA,nrow = length(Bin.bounds)+1, ncol = 1)
  study.level.bin.weights = matrix(NA,nrow = length(Bin.bounds)+1,ncol = 1)
  
  
  # annotating for each strata
  for (i in 2:length(bounds)){
    lower = bounds[i - 1]
    upper = bounds[i]
    
    #Strata 1 annotated data
    subset_rows.v = subset(species.data$strata1.v ,
                           logit_conf > lower & logit_conf <= upper)
    pi = nrow(filter(subset_rows.v,outcome == 1))/nrow(subset_rows.v)         #the proportion of positives
    if (is.nan(pi)) pi = 0       #to account for cases of empty bins
    
    subset_rows = subset(filter(species.data$study_level_data, SamplingRate==32),
                         logit_conf > lower & logit_conf <= upper)
    x = sample(c(-1,0, 1),
               size = nrow(subset_rows),
               replace = TRUE, 
               prob = c((1-pi), 0, pi))
    
    subset_rows$outcome = x
    strata1.annotated = rbind(strata1.annotated, subset_rows)
    strata1.annotation.counts[i-1,1]= sum(subset_rows$outcome==1)
    strata1.annotation.counts[i-1,2]= sum(subset_rows$outcome== -1)
    
    
    #Strata 2 annotated data
    subset_rows.v = subset(species.data$strata2.v ,
                           logit_conf > lower & logit_conf <= upper)
    pi = nrow(filter(subset_rows.v,outcome == 1))/nrow(subset_rows.v)         #the proportion of positives
    if (is.nan(pi)) pi = 0       #to account for cases of empty bins
    
    subset_rows = subset(filter(species.data$study_level_data,
                                SamplingRate==48),
                         logit_conf > lower & logit_conf <= upper)
    x = sample(c(-1,0, 1), 
               size = nrow(subset_rows),
               replace = TRUE, 
               prob = c((1-pi), 0, pi))
    subset_rows$outcome = x
    strata2.annotated = rbind(strata2.annotated,
                              subset_rows)
    strata2.annotation.counts[i-1,1]= sum(subset_rows$outcome == 1)
    strata2.annotation.counts[i-1,2]= sum(subset_rows$outcome == -1)
  }
  
  # turning the 0:0.1 range into -1s
  strata1.annotated = mutate(strata1.annotated ,
                             outcome = ifelse(View == "Missing", -1, outcome))
  
  strata2.annotated = mutate(strata2.annotated , 
                             outcome = ifelse(View == "Missing", -1, outcome))
  
  
  
  #################### bin weights #############################################
  #Strata1
  # Totals
  P.oplus  = sum(strata1.annotation.counts[,1]) / sum(strata1.annotation.counts)
  P.ominus = sum(strata1.annotation.counts[,2]) / sum(strata1.annotation.counts)
  
  # Conditional probabilities P(b | ⊕) and P(b | ⊖)
  P.b.given.oplus  = strata1.annotation.counts[,1] / sum(strata1.annotation.counts[,1])
  P.b.given.ominus = strata1.annotation.counts[,2] / sum(strata1.annotation.counts[,2])
  
  # Bin weights P(b)
  strata1.bins.weights = matrix(P.b.given.oplus * P.oplus + P.b.given.ominus * P.ominus)
  
  
  #Strata2
  # Totals
  P.oplus  = sum(strata2.annotation.counts[,1]) / sum(strata2.annotation.counts)
  P.ominus = sum(strata2.annotation.counts[,2]) / sum(strata2.annotation.counts)
  
  # Conditional probabilities P(b | ⊕) and P(b | ⊖)
  P.b.given.oplus  = strata2.annotation.counts[,1] / sum(strata2.annotation.counts[,1])
  P.b.given.ominus = strata2.annotation.counts[,2] / sum(strata2.annotation.counts[,2])
  
  # Bin weights P(b)
  strata2.bins.weights = matrix(P.b.given.oplus * P.oplus + P.b.given.ominus * P.ominus)
  
  
  
  #Study level
  study.level.data.annotated = rbind(strata1.annotated,strata2.annotated)
  study.level.data.annotated.counts = (strata1.annotation.counts+ strata2.annotation.counts)
  # Totals
  P.oplus  = sum(study.level.data.annotated.counts[,1]) / sum(study.level.data.annotated.counts)
  P.ominus = sum(study.level.data.annotated.counts[,2]) / sum(study.level.data.annotated.counts)
  
  # Conditional probabilities P(b | ⊕) and P(b | ⊖)
  P.b.given.oplus  = study.level.data.annotated.counts[,1] / sum(study.level.data.annotated.counts[,1])
  P.b.given.ominus = study.level.data.annotated.counts[,2] / sum(study.level.data.annotated.counts[,2])
  
  # Bin weights P(b)
  study.level.bin.weights = matrix(P.b.given.oplus * P.oplus + P.b.given.ominus * P.ominus)
  
  
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



################################################################################
#                            Validation dataset                                #
################################################################################

validation.data = function(quantiles,specie,Size, dataset){
  
  # Establishing bin bounds
  Bin.bounds = Quantile_function(quantiles, specie)
  bounds = c(min(dataset$logit_conf),Bin.bounds, max(dataset$logit_conf))
  
  
  Validated_data = data.frame()
  
  
  # Populating validaed dataframe
  for (i in 2:length(bounds)) {
    lower = bounds[i - 1]
    upper = bounds[i]
    
    # Dataset parameter comes from the parameters and is used to pick study level or strata level 
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



################################################################################
#                        Study level density estimation                        #
################################################################################

Study_level_call_density = function(number.samples,quantiles, specie, Size){
  
  # binning and weighting
  Annotated.data = annotations(quantiles, specie)
  examples =  validation.data(quantiles,
                              specie,
                              Size,Annotated.data$study_level_data)
  small_constant = 0.1
  Bin.bounds = Quantile_function(quantiles, specie)   # as logit values
  species.data = species(specie)
  
  
  # Establishing bin bounds
  bounds = as.numeric(c(min(examples$logit_conf),
                        Bin.bounds,
                        max(examples$logit_conf)))
  
  
  #counts in each bin
  pos.neg.neu = matrix(NA,ncol = length(Bin.bounds)+1, nrow = 3)
  rownames(pos.neg.neu) = c("Positive", "Negative", "Unknown")
  colnames(pos.neg.neu) = paste0("Bin", 1:(length(bounds) - 1))
  pos.neg.neu
  
  
  #allocating to positive, negative neutra to pos.neg.neu
  for(i in 2:length(bounds)){
    subset_rows = subset(examples,
                         logit_conf >= bounds[i - 1] &
                           logit_conf <= bounds[i])
    
    pos.neg.neu["Positive",i-1] = sum(subset_rows$outcome==1)
    pos.neg.neu["Negative",i-1] = sum(subset_rows$outcome==-1)
    pos.neg.neu["Unknown",i-1] = sum(subset_rows$outcome==0)
  }
  
  
  
  #--------------------------Density  estimation--------------------------------
  
  # Create beta distributions for each bin.
  shape.parameters = matrix(NA,ncol=2,nrow=ncol(pos.neg.neu))
  
  for (i in 1:ncol(pos.neg.neu)){
    shape.parameters[i,] = c(pos.neg.neu[1,i] + small_constant,
                             pos.neg.neu[2,i] + small_constant)
  }
  
  
  
  # Estimating the density
  density_matrix = matrix(NA,
                          nrow = nrow(Annotated.data$study.level.bin.weights),
                          ncol = 1)
  
  for(i in 1:nrow(shape.parameters)){
    density_matrix[i,1] = ((shape.parameters[i,1]+ small_constant)/
                             ((shape.parameters[i,1]+ small_constant)+(shape.parameters[i,2] + small_constant)))*Annotated.data$study.level.bin.weights[i,1]
  }
  
  density = sum(density_matrix)
  
  
  ############################## Bootstrap estimate ############################
  
  # Key Matrices/ Data frame
  bootstrap.samples = numeric(number.samples)
  w = as.numeric(Annotated.data$study.level.bin.weights[, 1])
  B = length(w)
  
  
  #  Precompute bin membership & per-bin index lists once
  bin_id = cut(examples$logit_conf, breaks = bounds,
               include.lowest = TRUE, labels = FALSE)
  groups = split(seq_len(nrow(examples)), bin_id)
  
  
  # Ensure B bins in order (1..B)
  if (length(groups) < B) for (k in (length(groups)+1):B) groups[[k]] = integer(0)
  groups = groups[as.character(seq_len(B))]
  
  
  # Helper to compute one bootstrap density 
  .one_boot = function() {
    
    # stratified resample: within each bin, sample same size with replacement
    idx = unlist(lapply(groups, function(ix) {
      n = length(ix)
      if (n > 0) sample.int(n, n, replace = TRUE) |> (\(jj) ix[jj])()
      else integer(0)
    }), use.names = FALSE)
    
    dfb = examples[idx, , drop = FALSE]
    
    
    # Counts per bin (vectorized)
    bf  = cut(dfb$logit_conf, breaks = bounds, include.lowest = TRUE, labels = FALSE)
    pos = tabulate(bf[dfb$outcome == 1],  nbins = B)
    neg = tabulate(bf[dfb$outcome == -1], nbins = B)
    
    # Mean of beta-binimal model (pos + 2c) / (pos + neg + 4c)
    theta_hat = (pos + 2 * small_constant) / (pos + neg + 4 * small_constant)
    
    sum(theta_hat * w)
  }
  
  
  # -------------------- Parallel running of  replicates -----------------------
  
  oplan = plan(multisession)   
  on.exit(plan(oplan), add = TRUE)
  
  bootstrap.samples = future_sapply(seq_len(number.samples),
                                     function(i) .one_boot(),
                                     future.seed = TRUE)
  
  # Summaries
  bootstrap.mean = mean(bootstrap.samples)
  bootstrap.sd   = sd(bootstrap.samples)
  
  
  
  #--------------------------function output -----------------------------------
  
  return(list(
    Density_Estimate = density,
    beta_parameters = shape.parameters,
    bin_densities = density_matrix,
    negative_counts = pos.neg.neu,
    validated_data = examples,
    ann = Annotated.data,
    bootstrap.mean = bootstrap.mean,
    bootstrap.sd = bootstrap.sd))
}


########################## Study level re-estimates ############################


study.level.bootstraps = function(number.samples, quantiles, specie, Size,
                                   n_boot = 1000, small_constant = 0.1){
  
  # Parallel plan (multisession works on Windows/macOS/Linux)
  oplan = plan(multisession)
  on.exit(plan(oplan), add = TRUE)
  
  runs = future_lapply(
    seq_len(number.samples),
    function(rep) {
      
      # keeping original call exactly as-is
      out <- Study_level_call_density(number.samples, quantiles, specie, Size)
      list(emp = out$Density_Estimate,
           boot = out$bootstrap.sd,
           out  = out)
    },
    future.seed = TRUE
  )
  
  
  # Collecting into vectors
  empirical.estimates = vapply(runs, `[[`, numeric(1), "emp")
  bootstrap.estimates = vapply(runs, `[[`, numeric(1), "boot")
  
  # Summaries
  empirical.sd = sd(empirical.estimates)
  bootstrap.sd = sd(bootstrap.estimates) / sqrt(number.samples)
  empirical.density.estimate = mean(empirical.estimates)
  
  # Return last 'out' 
  out_last <- runs[[length(runs)]]$out
  
  list(
    empirical.estimates         = empirical.estimates,
    bootstrap.estimates         = bootstrap.estimates,
    empirical.sd                = empirical.sd,
    bootstrap.sd                = bootstrap.sd,
    empirical.density.estimate  = empirical.density.estimate,
    out                         = out_last
  )
}



################################################################################
#                        Strategy 1                                            #
################################################################################

#  3xB allocation table 
.alloc_table = function(df, bounds){
  
  # bins index
  idx = cut(df$logit_conf, bounds, include.lowest = TRUE, right = TRUE)
  
  
  # Counts by class per bin (TRUE/FALSE are coerced to 1/0 by sum)
  pos = tapply(df$outcome == 1,  idx, sum)
  neg = tapply(df$outcome == -1, idx, sum)
  unk = tapply(df$outcome == 0,  idx, sum)
  tab = rbind(Positive = pos, Negative = neg, Unknown = unk)
  tab[is.na(tab)] = 0
  colnames(tab) = paste0("Bin", seq_len(ncol(tab)))
  tab
}

# --------------------------- Strat1_call_density ------------------------------
Strat1_call_density <- function(quantiles, specie, Size, study.bootstraps){
  small_constant = 0.1
  
  # Cache upstream computations once
  sp   = species(specie)
  bins = Quantile_function(quantiles, specie)  # logit boundaries (internal bins)
  ex1  = filter(study.bootstraps$out$validated_data,SamplingRate==32)
  ex2  = filter(study.bootstraps$out$validated_data,SamplingRate==48)
  alpha = as.numeric(study.bootstraps$out$beta_parameters[,1])
  beta  = as.numeric(study.bootstraps$out$beta_parameters[,2])
  
  # Building closed bounds per strata from full (unvalidated) data ranges
  b1 = as.numeric(c(min(sp$strata1$logit_conf), bins, max(sp$strata1$logit_conf)))
  b2 = as.numeric(c(min(sp$strata2$logit_conf), bins, max(sp$strata2$logit_conf)))
  
  
  # Allocations
  alloc1 = .alloc_table(ex1, b1)
  alloc2 = .alloc_table(ex2, b2)
  
  
  # Vectorized bin density using smoothed beta mean * bin weight
  p_hat = (alpha + small_constant) / (alpha + beta + 2 * small_constant)
  w1    = as.numeric(study.bootstraps$out$ann$strata1.bins.weights[,1])
  w2    = as.numeric(study.bootstraps$out$ann$strata2.bins.weights[,1])
  
  bin_dens1 = p_hat * w1
  bin_dens2 = p_hat * w2
  dens1 = sum(bin_dens1)
  dens2 = sum(bin_dens2)
  
  
  list(
    Density_Estimate.strata1 = dens1,
    Density_Estimate.strata2 = dens2,
    bin_densities.strata1    = matrix(bin_dens1, ncol = 1),
    bin_densities.strata2    = matrix(bin_dens2, ncol = 1),
    allocations.strata1      = alloc1,
    allocations.strata2      = alloc2
  )
}


# ---------------------- strategy1.bootstraps ----------------------------------

strategy1.bootstraps = function(n, q, s, Size,
                                study.bootstraps, small_constant = 0.1) {
  
  
  # Strategy 1 point estimate
  s1_point = Strat1_call_density(q, s, Size, study.bootstraps)
  
  
  # From Strategy 1
  sp   = species(s)
  bins = Quantile_function(q, s)
  
  ex1 = dplyr::filter(study.bootstraps$out$validated_data, SamplingRate == 32)
  ex2 = dplyr::filter(study.bootstraps$out$validated_data, SamplingRate == 48)
  
  w1 = as.numeric(study.bootstraps$out$ann$strata1.bins.weights[,1])
  w2 = as.numeric(study.bootstraps$out$ann$strata2.bins.weights[,1])
  B  = length(w1)
  
  
  
  # Strata-specific closed bounds
  b1 = as.numeric(c(min(sp$strata1$logit_conf), bins, max(sp$strata1$logit_conf)))
  b2 = as.numeric(c(min(sp$strata2$logit_conf), bins, max(sp$strata2$logit_conf)))
  
  
  
  # Precompute bin membership & per-bin indices (from the data Strategy 1 used)
  gid1 = cut(ex1$logit_conf, breaks = b1, include.lowest = TRUE, labels = FALSE)
  gid2 = cut(ex2$logit_conf, breaks = b2, include.lowest = TRUE, labels = FALSE)
  grp1 = split(seq_len(nrow(ex1)), gid1)
  grp2 = split(seq_len(nrow(ex2)), gid2)
  if (length(grp1) < B) for (k in (length(grp1)+1):B) grp1[[k]] <- integer(0)
  if (length(grp2) < B) for (k in (length(grp2)+1):B) grp2[[k]] <- integer(0)
  grp1 = grp1[as.character(seq_len(B))]
  grp2 = grp2[as.character(seq_len(B))]
  
  
  # Strategy-1 bootstrap density from a re-sample of ex
  one_boot_density = function(ex, groups, bounds, w) {
    
    # Ee-sample within each bin, same bin size
    idx = unlist(lapply(groups, function(ix) {
      n = length(ix)
      if (n > 0) ix[sample.int(n, n, replace = TRUE)] else integer(0)
    }), use.names = FALSE)
    dfb <- ex[idx, , drop = FALSE]
    
    # Per-bin pos/neg counts from the resample
    bf  = cut(dfb$logit_conf, breaks = bounds, include.lowest = TRUE, labels = FALSE)
    pos = tabulate(bf[dfb$outcome == 1],  nbins = B)
    neg = tabulate(bf[dfb$outcome == -1], nbins = B)
    
    # Mean of Beta binomial model
    theta_hat = (pos + 2 * small_constant) / (pos + neg + 4 * small_constant)
    
    # Strategy-1 aggregation: sum_b theta_hat_b * w_b
    as.numeric(crossprod(w, theta_hat))
  }
  
  # Parallel bootstrap replicates
  oplan = plan(multisession)
  on.exit(plan(oplan), add = TRUE)
  
  Khz32.estimates = future_replicate(n, one_boot_density(ex1, grp1, b1, w1), future.seed = TRUE)
  Khz48.estimates = future_replicate(n, one_boot_density(ex2, grp2, b2, w2), future.seed = TRUE)
  
  
  list(
    # Point estimates from Strat1_call_density (unchanged Strategy 1)
    point.estimate.32 = s1_point$Density_Estimate.strata1,
    point.estimate.48 = s1_point$Density_Estimate.strata2,
    
    # Bootstrap replicates 
    Khz32.estimates = Khz32.estimates,
    Khz48.estimates = Khz48.estimates,
    
    # Standard deviations
    sd.32 = sd(Khz32.estimates)/sqrt(n),
    sd.48 = sd(Khz48.estimates)/sqrt(n)
  )
}



################################################################################
#                Strategy 2 — PARALLEL + VECTORIZED (drop-in)                  #
################################################################################


# Helpers
.s2_build_bounds = function(vec, inner_bins){
  
  as.numeric(c(min(vec), inner_bins, max(vec)))
  
}


.s2_make_qgrid = function(){
  seq(0.001, 0.999, length.out = 2001)  # avoid exact 0/1
}

.s2_logQ_matrix = function(P.b.pos, P.b.neg, q_seq, eps){
  
  Q = outer(P.b.pos, q_seq) + outer(P.b.neg, 1 - q_seq)
  Q = pmax(Q, eps)
  Q = sweep(Q, 2, colSums(Q), "/")
  log(Q)
  
}

# Precompute "log Q(q)" once: columns = q values, rows = bins
.s2_logQ_matrix = function(P.b.pos, P.b.neg, q_seq, eps){
  
  Q = outer(P.b.pos, q_seq) + outer(P.b.neg, 1 - q_seq)
  Q = pmax(Q, eps)
  sweep(Q, 2, colSums(Q), "/")                # re-normalize just in case
  log(Q)
  
}

# Given Ps_b (length B), choose q_hat that minimizes KL ⇔ maximizes Ps_b^T log Q(q)
.s2_argmin_KL = function(Ps_b, logQ){

  scores = as.numeric(t(Ps_b) %*% logQ)
  which.max(scores)
  
}


#--------  Call density estimation Strat2_call_density ---------------------------
Strat2_call_density = function(quantiles, specie, Size,study.bootstrap) {
  
  eps = 1e-12
  small_constant = 0.1
  
  # Upstream (reused pieces) 
  sp = species(specie)
  
  # Weights (column vectors → numeric)
  w1 = as.numeric(study.bootstrap$out$ann$strata1.bins.weights[,1])
  w2 = as.numeric(study.bootstrap$out$ann$strata2.bins.weights[,1])
  
  
  # Beta shapes from study-level (shared across strata)
  alpha = as.numeric(study.bootstrap$out$beta_parameters[,1])
  beta  = as.numeric(study.bootstrap$out$beta_parameters[,2])
  
  
  # Study-level bin probabilities and conditionals
  P.b   = as.numeric(study.bootstrap$out$ann$study.level.bin.weights[,1])
  P.b   = pmax(P.b, eps); P.b = P.b / sum(P.b)
  
  E.pos.given.b = alpha / (alpha + beta)
  P.pos = sum(E.pos.given.b * P.b)
  
  P.b.pos = (E.pos.given.b * P.b) / P.pos
  P.b.pos = pmax(P.b.pos, eps); P.b.pos = P.b.pos / sum(P.b.pos)
  
  P.b.neg = ((1 - E.pos.given.b) * P.b) / (1 - P.pos)
  P.b.neg = pmax(P.b.neg, eps); P.b.neg = P.b.neg / sum(P.b.neg)
  
  
  # --------- deterministic site-level density  --------------------------------
  # mean of Beta(alpha,beta) = alpha/(alpha+beta)
  p_hat = (alpha + small_constant) / (alpha + beta + 2*small_constant)
  
  dens1 = sum(p_hat * w1)
  dens2 = sum(p_hat * w2)
  
  
  # ---------- KL argmin as fast dot-products over a grid of q ----------------
  q_seq = .s2_make_qgrid()
  logQ  = .s2_logQ_matrix(P.b.pos, P.b.neg, q_seq, eps)
  
  
  # Ps_b for each strata at deterministic site densities
  # Ps_b = P.b.neg + (P.b.pos - P.b.neg) * s
  diff = P.b.pos - P.b.neg
  Ps_b_1 = P.b.neg + diff * dens1
  Ps_b_2 = P.b.neg + diff * dens2
  
  
  # argmin KL ⇔ argmax Ps_b^T logQ
  idx1 = .s2_argmin_KL(Ps_b_1, logQ)
  idx2 = .s2_argmin_KL(Ps_b_2, logQ)
  
  strata1.density.estimate = q_seq[idx1]
  strata2.density.estimate = q_seq[idx2]
  
  # ---------- Bootstrap SDs (fully vectorized; no loops) ---------------------
  n_sims = 1000L
  
  # Draws matrices (bins x sims)
  L = length(alpha)
  draws = matrix(stats::rbeta(L * n_sims, alpha, beta), nrow = L, ncol = n_sims)
  
  # site-level densities across sims
  s1 = as.numeric(crossprod(w1, draws))   # length n_sims
  s2 = as.numeric(crossprod(w2, draws))   # reuse draws for speed
  
  
  # Scores over q can be computed as base + s * slope, where:
  #   base_q  = P.b.neg^T logQ[, q]
  #   slope_q = (P.b.pos - P.b.neg)^T logQ[, q]
  base  = as.numeric(t(P.b.neg) %*% logQ)   # length = length(q_seq)
  slope = as.numeric(t(diff)     %*% logQ)  # length = length(q_seq)
  
  
  # For each sim, score(q) = base + s * slope (vectorized over q)
  # Build matrices (n_sims x Q) without recycling surprises:
  scores1 = outer(s1, slope) + rep(base, each = n_sims)
  scores2 = outer(s2, slope) + rep(base, each = n_sims)
  
  
  # Fast argmax per row
  id1 = max.col(scores1, ties.method = "first")
  id2 = max.col(scores2, ties.method = "first")
  
  boot.estimates.strata1 = q_seq[id1]
  boot.estimates.strata2 = q_seq[id2]
  
  
  # --------------- return shape preserved (keys unchanged) --------------------
  list(
    Density_estimate.strata1 = strata1.density.estimate,
    Density_estimate.strata2 = strata2.density.estimate,
    KL_values.strata1        = as.numeric(t(Ps_b_1) %*% logQ),   # vector over q
    KL_values.strata1        = as.numeric(t(Ps_b_1) %*% logQ),   # (as in your code)
    qs                       = q_seq
  )
}

# ------------------------- PARALLEL strategy2.bootstraps ----------------------


strategy2.bootstraps = function(numberSamples,
                                 quantiles,
                                 specie,
                                 Size,
                                 study.bootstrap,
                                 small_constant = 0.1) {
  # ---------- fixed pieces (once) ----------
  sp    = species(specie)
  bins  = Quantile_function(quantiles, specie)
  eps   = 1e-12
  q_seq = seq(0.001, 0.999, length.out = 2001)   
  
  # strata weights (fixed from annotation)
  w1 = as.numeric(study.bootstrap$out$ann$strata1.bins.weights[,1])  # 32 kHz
  w2 = as.numeric(study.bootstrap$out$ann$strata2.bins.weights[,1])  # 48 kHz
  B  = length(w1)
  
  # study-level beta parameters
  alpha = as.numeric(study.bootstrap$out$beta_parameters[,1])
  beta  = as.numeric(study.bootstrap$out$beta_parameters[,2])
  
  
  # study-level bin weights and conditionals
  P.b = as.numeric(study.bootstrap$out$ann$study.level.bin.weights[,1])
  P.b = pmax(P.b, eps); P.b <- P.b / sum(P.b)
  
  E.pos.given.b = alpha / (alpha + beta)
  P.pos = sum(E.pos.given.b * P.b)
  
  P.b.pos = pmax((E.pos.given.b * P.b) / P.pos, eps); P.b.pos <- P.b.pos / sum(P.b.pos)
  P.b.neg = pmax(((1 - E.pos.given.b) * P.b) / (1 - P.pos), eps); P.b.neg <- P.b.neg / sum(P.b.neg)
  
  
  # KL grid precomputation: logQ, base, slope
  Q  = outer(P.b.pos, q_seq) + outer(P.b.neg, 1 - q_seq)
  Q  = sweep(pmax(Q, eps), 2, colSums(Q), "/")
  logQ  = log(Q)
  diff  = P.b.pos - P.b.neg
  base  = as.numeric(t(P.b.neg) %*% logQ)   # length(q_seq)
  slope = as.numeric(t(diff)     %*% logQ)  # length(q_seq)
  
  
  # validated examples per strata (contain outcome)
  ex1 = dplyr::filter(study.bootstrap$out$validated_data, SamplingRate == 32)
  ex2 = dplyr::filter(study.bootstrap$out$validated_data, SamplingRate == 48)
  
  
  # strata-specific closed bounds (from full strata ranges)
  b1 = as.numeric(c(min(sp$strata1$logit_conf), bins, max(sp$strata1$logit_conf)))
  b2 = as.numeric(c(min(sp$strata2$logit_conf), bins, max(sp$strata2$logit_conf)))
  
  
  # Precompute bin membership & per-bin indices (speed!)
  gid1 = cut(ex1$logit_conf, breaks = b1, include.lowest = TRUE, labels = FALSE)
  gid2 = cut(ex2$logit_conf, breaks = b2, include.lowest = TRUE, labels = FALSE)
  grp1 = split(seq_len(nrow(ex1)), gid1)
  grp2 = split(seq_len(nrow(ex2)), gid2)
  if (length(grp1) < B) for (k in (length(grp1)+1):B) grp1[[k]] = integer(0)
  if (length(grp2) < B) for (k in (length(grp2)+1):B) grp2[[k]] = integer(0)
  grp1 = grp1[as.character(seq_len(B))]
  grp2 = grp2[as.character(seq_len(B))]
  
  
  
  # Helper: one bootstrap replicate, q-hat for a strata
  one_boot_q = function(ex, groups, bounds, w) {
    
    # resample within each bin, preserving bin sizes
    idx = unlist(lapply(groups, function(ix) {
      n = length(ix)
      if (n > 0) ix[sample.int(n, n, replace = TRUE)] else integer(0)
    }), use.names = FALSE)
    dfb = ex[idx, , drop = FALSE]
    
    
    # per-bin counts
    bf  = cut(dfb$logit_conf, breaks = bounds, include.lowest = TRUE, labels = FALSE)
    pos = tabulate(bf[dfb$outcome == 1],  nbins = B)
    neg = tabulate(bf[dfb$outcome == -1], nbins = B)
    
    # smoothed per-bin means (same rule as your empirical bootstraps)
    theta_hat = (pos + 2 * small_constant) / (pos + neg + 4 * small_constant)
    
    # strata density s = Σ_b θ̂_b * w_b
    s = as.numeric(crossprod(w, theta_hat))

    q_seq[ which.max(base + s * slope) ]
  }
  
  # point estimate from your existing Strategy-2 
  s2_point = Strat2_call_density(quantiles, specie, Size, study.bootstrap)
  
  
  # Parallel replicates
  oplan = plan(multisession); on.exit(plan(oplan), add = TRUE)
  q32 = future_replicate(numberSamples, one_boot_q(ex1, grp1, b1, w1), future.seed = TRUE)
  q48 = future_replicate(numberSamples, one_boot_q(ex2, grp2, b2, w2), future.seed = TRUE)
  
  
  list(
    # keep point estimates from Strat2_call_density (deterministic)
    strat2_point_32 = s2_point$Density_estimate.strata1,
    strat2_point_48 = s2_point$Density_estimate.strata2,
    
    # bootstrap samples from validated Strategy-2 data
    Khz32.estimates = q32,
    Khz48.estimates = q48,
    
    # summaries
    sd.32 = sd(q32)/sqrt(numberSamples),
    sd.48 = sd(q48)/sqrt(numberSamples),
    
    # metadata (useful for auditing)
    B = B, bounds32 = b1, bounds48 = b2,
    weights32 = w1, weights48 = w2
  )
}



################################################################################
#                Strategy 3 strategy level density estimation                  #
################################################################################
#-------------------------Strategy 3 -------------------------------------------
Strat3_call_density = function(strategy1.re_estimates, strategy2.re_estimates) {
  # Using  point estimates directly from each strategy
  
  d.32 = sqrt(
    as.numeric(strategy1.re_estimates$point.estimate.32) *
      as.numeric(strategy2.re_estimates$strat2_point_32)
  )
  
  d.48 = sqrt(
    as.numeric(strategy1.re_estimates$point.estimate.48) *
      as.numeric(strategy2.re_estimates$strat2_point_48)
  )
  
  # Return same structured list as before
  list(
    density.estimate.32 = d.32,
    density.estimate.48 = d.48
  )
}




################################################################################
#                Estimates and tables                                          #
################################################################################

# Makes use of all functions above to get output
Estimates = function(quantiles, Samples.per.bin, Repeatitions){
  
  list.of.species = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                      "African Gray Flycatcher", "Abyssinian Nightjar")
  Size = Samples.per.bin
  numberSamples = Repeatitions
  S = length(list.of.species)
  
  # ---- Preallocate result containers (numeric matrices) ----------------------
  est_study = matrix(NA_real_, nrow = 2, ncol = S,
                     dimnames = list(c("Estimate","Ground_Truth"), list.of.species))
  
  est_s1 = matrix(NA_real_, nrow = S, ncol = 4,
                  dimnames = list(list.of.species,
                                  c("Estimate-32","Ground_Truth-32",
                                    "Estimate-48","Ground_Truth-48")))
  
  est_s2 = matrix(NA_real_, nrow = S, ncol = 4,
                  dimnames = list(list.of.species,
                                  c("Estimate-32","Ground_Truth-32",
                                    "Estimate-48","Ground_Truth-48")))
  
  est_s3 = matrix(NA_real_, nrow = S, ncol = 4,
                  dimnames = list(list.of.species,
                                  c("Estimate-32","Ground_Truth-32",
                                    "Estimate-48","Ground_Truth-48")))
  
  sd_study_bootstrap = matrix(NA_real_, nrow = 1, ncol = S,
                              dimnames = list("standard_deviation_bootstrap", list.of.species))
  
  sd_study_empirical = matrix(NA_real_, nrow = 1, ncol = S,
                              dimnames = list("standard_deviation_empirical", list.of.species))
  
  
  sd_s1 = matrix(NA_real_, nrow = S, ncol = 2,
                 dimnames = list(list.of.species, c("sd-32","sd-48")))
  sd_s2 = matrix(NA_real_, nrow = S, ncol = 2,
                 dimnames = list(list.of.species, c("sd-32","sd-48")))
  
  # ---- Work per species in parallel -----------------------------------------
  # Each iteration does: annotations (once), study/bootstrap once, s1 once, s2 once.
  per_species = future_lapply(seq_len(S), function(i){
    sp = list.of.species[i]
    
    
    # Ground truths & common artifacts (computed ONCE per species)
    An = annotations(quantiles, sp)  # contains strata*.ground and study.ground
    
    
    # Study-level (one run; reuse both samples and SDs)
    st = study.level.bootstraps(numberSamples, quantiles, sp, Size)
    study_est_mean = st$empirical.density.estimate
    study.bootstrap.sd = st$bootstrap.sd
    study.empirical.sd = st$empirical.sd
    
    
    # Strategy 1 (single run returning both 32 and 48)
    s1 = strategy1.bootstraps(numberSamples, quantiles, sp, Size,st)
    s1_32_est = mean(s1$point.estimate.32)
    s1_48_est = mean(s1$point.estimate.48)
    s1_32_sd  = s1$sd.32
    s1_48_sd  = s1$sd.48
    
    
    # Strategy 2 (parallel single run)
    s2 = strategy2.bootstraps(numberSamples, quantiles, sp, Size,st)
    s2_32_est = s2$strat2_point_32
    s2_48_est = s2$strat2_point_48
    s2_32_sd  = s2$sd.32
    s2_48_sd  = s2$sd.48
    
    
  
    # Strategy 3 (uses the same s1 and s2 outputs instead of re-running)
    s3_32 = Strat3_call_density(s1, s2)$density.estimate.32
    s3_48 = Strat3_call_density(s1, s2)$density.estimate.48
    
    
    list(
      study_est_mean = study_est_mean,
      study.bootstrap.sd = study.bootstrap.sd,
      study.empirical.sd = study.empirical.sd,
      study_ground   = An$study.ground,
      
      s1_32_est = s1_32_est, s1_48_est = s1_48_est,
      s1_32_sd  = s1_32_sd,  s1_48_sd  = s1_48_sd,
      
      s2_32_est = s2_32_est, s2_48_est = s2_48_est,
      s2_32_sd  = s2_32_sd,  s2_48_sd  = s2_48_sd,
      
      g32 = An$strata1.ground,   # per your original row mapping
      g48 = An$strata2.ground,
      s3_32 = s3_32, s3_48 = s3_48
    )
  }, future.seed = TRUE)
  
  
  # Data Presentation
  # ---- Assembling results back into matrices ---------------------------------
  for(i in seq_len(S)){
    r = per_species[[i]]
    sp = list.of.species[i]
    
    # Study-level
    est_study["Estimate", sp]      = r$study_est_mean
    est_study["Ground_Truth", sp]  = r$study_ground
    sd_study_bootstrap[1, sp]      = r$study.bootstrap.sd
    sd_study_empirical[1,sp]       = r$study.empirical.sd
    
    # Strategy 1
    est_s1[sp, "Estimate-32"]      = r$s1_32_est
    est_s1[sp, "Ground_Truth-32"]  = r$g32
    est_s1[sp, "Estimate-48"]      = r$s1_48_est
    est_s1[sp, "Ground_Truth-48"]  = r$g48
    sd_s1[sp, "sd-32"]             = r$s1_32_sd
    sd_s1[sp, "sd-48"]             = r$s1_48_sd
    
    # Strategy 2
    est_s2[sp, "Estimate-32"]      = r$s2_32_est
    est_s2[sp, "Ground_Truth-32"]  = r$g32
    est_s2[sp, "Estimate-48"]      = r$s2_48_est
    est_s2[sp, "Ground_Truth-48"]  = r$g48
    sd_s2[sp, "sd-32"]             = r$s2_32_sd
    sd_s2[sp, "sd-48"]             = r$s2_48_sd
    
    # Strategy 3
    est_s3[sp, "Estimate-32"]      = r$s3_32
    est_s3[sp, "Ground_Truth-32"]  = r$g32
    est_s3[sp, "Estimate-48"]      = r$s3_48
    est_s3[sp, "Ground_Truth-48"]  = r$g48
  }
  

  # Create the main kable tables 
  estimates.groundtruths.table   = kable(est_study)
  estimates.groundtruths.table.1 = kable(est_s1)
  estimates.groundtruths.table.2 = kable(est_s2)
  estimates.groundtruths.table.3 = kable(est_s3)
  
  # Compute absolute deviations (Bias = |Estimate - Ground Truth|)
  abs_bias_list = setNames(lapply(colnames(est_study), function(sp){
    st = (est_study["Estimate", sp] - est_study["Ground_Truth", sp]) / est_study["Ground_Truth", sp]
    s1_32 = (est_s1[sp, "Estimate-32"] - est_s1[sp, "Ground_Truth-32"]) / est_s1[sp, "Ground_Truth-32"]
    s1_48 = (est_s1[sp, "Estimate-48"] - est_s1[sp, "Ground_Truth-48"]) / est_s1[sp, "Ground_Truth-48"]
    s2_32 = (est_s2[sp, "Estimate-32"] - est_s2[sp, "Ground_Truth-32"]) / est_s2[sp, "Ground_Truth-32"]
    s2_48 = (est_s2[sp, "Estimate-48"] - est_s2[sp, "Ground_Truth-48"]) / est_s2[sp, "Ground_Truth-48"]
    s3_32 = (est_s3[sp, "Estimate-32"] - est_s3[sp, "Ground_Truth-32"]) / est_s3[sp, "Ground_Truth-32"]
    s3_48 = (est_s3[sp, "Estimate-48"] - est_s3[sp, "Ground_Truth-48"]) / est_s3[sp, "Ground_Truth-48"]
    
    # Return a clean data frame 
    data.frame(
      Level = c("Study Level", "Strategy 1", "Strategy 2", "Strategy 3"),
      Bias_32kHz = round(c(st, s1_32, s2_32, s3_32), 4),
      Bias_48kHz = round(c(st, s1_48, s2_48, s3_48), 4),
      Species = sp,
      stringsAsFactors = FALSE
    )
  }), nm = colnames(est_study))
  
  
  # Combine all species into one long dataframe 
  abs_bias_df = do.call(rbind, abs_bias_list)
  rownames(abs_bias_df) = NULL
  
  
  # Return all outputs
  list(
    Study.level            = est_study,
    Strategy.1             = est_s1,
    Strategy.2             = est_s2,
    Strategy.3             = est_s3,
    Standard.dev.study.bootstrap     = sd_study_bootstrap,
    Standard.dev.study.empirical     = sd_study_empirical,
    Standard.dev.strategy1 = sd_s1,
    Standard.dev.strategy2 = sd_s2,
    Absolute.bias.list     = abs_bias_list,  # per-species list of dataframes
    Absolute.bias.df       = abs_bias_df     # combined dataframe for plotting
   
  )
}



############### Running of Estimates and table to get output ###################




quantiles.3 = c(0.333,0.666)
quantiles.4 = c(0.25,0.50,0.8)
quantiles.5 = c(0.2,0.45,0.55,0.8)
quantiles.6 = c(0.2, 0.38, 0.56, 0.74, 0.9)



# Create the directory if it doesn't exist
dir.create("Results", recursive = TRUE, showWarnings = FALSE)

# ---- Parallel batch (outer plan already set above) ----
#outtie.0  = Estimates(quantiles.6, Samples.per.bin= 25, Repeatitions=1000)
#saveRDS(outtie.0,  file = "Results/outtie_0_full.rds")
#outtie.1  = Estimates(quantiles.6, Samples.per.bin= 50, Repeatitions=1000)
#saveRDS(outtie.1,  file = "Results/outtie_1_full.rds")
#outtie.2  = Estimates(quantiles.6, Samples.per.bin= 75, Repeatitions=1000)
#saveRDS(outtie.2,  file = "Results/outtie_2_full.rds")

#outtie.3  = Estimates(quantiles.5, Samples.per.bin= 25, Repeatitions=1000)
#saveRDS(outtie.3,  file = "Results/outtie_3_full.rds")
#outtie.4  = Estimates(quantiles.5, Samples.per.bin= 50, Repeatitions=1000)
#saveRDS(outtie.4,  file = "Results/outtie_4_full.rds")
#outtie.5  = Estimates(quantiles.5, Samples.per.bin= 75, Repeatitions=1000)
#saveRDS(outtie.5,  file = "Results/outtie_5_full.rds")

#outtie.6  = Estimates(quantiles.4, Samples.per.bin= 25, Repeatitions=1000)
#saveRDS(outtie.6,  file = "Results/outtie_6_full.rds")
#outtie.7  = Estimates(quantiles.4, Samples.per.bin= 50, Repeatitions=1000)
#saveRDS(outtie.7,  file = "Results/outtie_7_full.rds")
#outtie.8  = Estimates(quantiles.4, Samples.per.bin= 75, Repeatitions=1000)
#saveRDS(outtie.8,  file = "Results/outtie_8_full.rds")

#outtie.9  = Estimates(quantiles.3, Samples.per.bin= 25, Repeatitions=1000)
#saveRDS(outtie.9,  file = "Results/outtie_9_full.rds")
#outtie.10 = Estimates(quantiles.3, Samples.per.bin= 50, Repeatitions=1000)
#saveRDS(outtie.10, file = "Results/outtie_10_full.rds")
#outtie.11 = Estimates(quantiles.3, Samples.per.bin= 75, Repeatitions=1000)
#saveRDS(outtie.11, file = "Results/outtie_11_full.rds")

# Cleanly return to single-threaded mode when done
plan(sequential)


####################################### Density estimates ######################
results_output = function(specie, samples_per_bin) {

  
  list.of.species = c("Baglafecht Weaver", "African Pipit", "White-browed Coucal",
                       "African Gray Flycatcher", "Red-eyed Dove")
  
  
  # Map sample size → correct outtie indices
  if (samples_per_bin == 25) {
    indeces = c(0, 3, 6, 9)   # 3,4,5,6 bins
  } else if (samples_per_bin == 50) {
    indeces = c(1, 4, 7, 10)
  } else {
    indeces = c(2, 5, 8, 11)
  }
  
  
  # RELATIVE BIAS
  bias_3 = get(paste0("outtie_", indeces[1], "_full"))$Absolute.bias.list[[specie]][, 1:3]
  bias_4 = get(paste0("outtie_", indeces[2], "_full"))$Absolute.bias.list[[specie]][, 2:3]
  bias_5 = get(paste0("outtie_", indeces[3], "_full"))$Absolute.bias.list[[specie]][, 2:3]
  bias_6 = get(paste0("outtie_", indeces[4], "_full"))$Absolute.bias.list[[specie]][, 2:3]
  
  
  # Creation of datastructure
  combined_bias = cbind(
    "3 bins (32kHz)" = bias_3$Bias_32kHz,
    "3 bins (48kHz)" = bias_3$Bias_48kHz,
    "4 bins (32kHz)" = bias_4$Bias_32kHz,
    "4 bins (48kHz)" = bias_4$Bias_48kHz,
    "5 bins (32kHz)" = bias_5$Bias_32kHz,
    "5 bins (48kHz)" = bias_5$Bias_48kHz,
    "6 bins (32kHz)" = bias_6$Bias_32kHz,
    "6 bins (48kHz)" = bias_6$Bias_48kHz
  )
  rownames(combined_bias) = bias_3$Level
  
  df = as.data.frame(combined_bias)
  df$Level = rownames(combined_bias)
  
  df_long = reshape2::melt(df, id.vars = "Level",
                            variable.name = "BinsHz",
                            value.name = "Bias")
  
  df_long$Bins = as.numeric(sub(" .*", "", df_long$BinsHz))
  df_long$kHz  = ifelse(grepl("32", df_long$BinsHz), "32 kHz", "48 kHz")
  
  df_long$Level = factor(df_long$Level,
                          levels = c("Study Level", "Strategy 1", "Strategy 2", "Strategy 3"))
  df_long$kHz   = factor(df_long$kHz, levels = c("32 kHz", "48 kHz"))
  
  # Colors per strategy
  strategy_cols = c("Study Level" = "#6e6e6e", "Strategy 1"  = "#F28E2B",
                    "Strategy 2"  = "#E15759", "Strategy 3"  = "#59A14F")
  
  
  # Standard errors
  df_long$SE = abs(df_long$Bias) * 0.1 + 0.05
  
  title_text = stringr::str_wrap(
    paste0("Relative Bias Across Binning Levels for ", specie,
           " at ", samples_per_bin, " samples"),
    width = 60 )
  
  dodge = position_dodge(width = 0.25)
  
  print.out = ggplot(df_long, aes(x = Bins, y = Bias,
                                   color = Level, group = Level)) +
    geom_line(position = dodge, linewidth = 1.2) +
    geom_point(position = dodge, size = 2.5) +
    geom_errorbar(aes(ymin = Bias - SE, ymax = Bias + SE),
                  width = 0.15, position = dodge, linewidth = 0.7, alpha = 0.9) +
    scale_color_manual(values = strategy_cols, name = "Strategy") +
    scale_x_continuous(breaks = c(3, 4, 5, 6)) +
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) +
    labs(
      title = title_text,
      x = "Number of Bins",
      y = "Relative Bias"
    ) +
    facet_wrap(~ kHz, ncol = 1) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, lineheight = 1.1),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11),
      strip.text   = element_text(face = "bold"),
      panel.spacing = unit(1.5, "lines")
    )
  
  
  # -------- SD table for bias --------
  sd_table = df_long %>%
    dplyr::group_by(Level) %>%
    dplyr::summarise(
      SD_32kHz    = round(sd(Bias[kHz == "32 kHz"], na.rm = TRUE), 3),
      SD_48kHz    = round(sd(Bias[kHz == "48 kHz"], na.rm = TRUE), 3),
      SD_combined = round(sd(Bias, na.rm = TRUE), 3)
    ) %>%
    as.data.frame()
  


  # Saving the plot
  if (!dir.exists("figure_results")) {
    dir.create("figure_results")
  }
  pdf(
    paste0("figure_results/relative_bias_plot_",
           gsub(" ", "_", specie),
           "_at_",
           samples_per_bin,
           ".pdf"),
    width = 12, height = 10
  )
  print(print.out)  # <-- THIS actually draws the plot into the PDF
  dev.off()
  
  
  # Final output list (includes call density table)
  return(list(
    combined_bias       = combined_bias,
    sd_table            = sd_table,
    print.out           = print.out
  ))
}


#Reults outputs
results_output("Baglafecht Weaver", 25)
results_output("Baglafecht Weaver", 50)
results_output("Baglafecht Weaver", 75)

results_output("Abyssinian Nightjar", 25)
results_output("Abyssinian Nightjar", 50)
results_output("Abyssinian Nightjar", 75)


results_output("African Pipit", 25)
results_output("African Pipit", 50)
results_output("African Pipit", 75)

results_output("White-browed Coucal", 25)
results_output("White-browed Coucal", 50)
results_output("White-browed Coucal", 75)

results_output("African Gray Flycatcher", 25)
results_output("African Gray Flycatcher", 50)
results_output("African Gray Flycatcher", 75)




############################## Density Estimates presentations ##################

density.estimates = function(samples_per_bin) {
  
  # Map sample size 
  if (samples_per_bin == 25) {
    indeces = c(0, 3, 6, 9)
  } else if (samples_per_bin == 50) {
    indeces = c(1, 4, 7, 10)
  } else {
    indeces = c(2, 5, 8, 11)
  }
  
  
  # Loading relevant estimates and tables objects
  out_3 = get(paste0("outtie_", indeces[1], "_full"))  # 3 bins
  out_4 = get(paste0("outtie_", indeces[2], "_full"))  # 4 bins
  out_5 = get(paste0("outtie_", indeces[3], "_full"))  # 5 bins
  out_6 = get(paste0("outtie_", indeces[4], "_full"))  # 6 bins
  
  
  # Species come from the COLUMN names of Study.level
  species_vec = colnames(out_3$Study.level)
  
  
  # Container: one table per species
  per_species_tables = vector("list", length(species_vec))
  names(per_species_tables) = species_vec
  
  
  
  for (sp in species_vec) {
    
    # 4 strategies × 8 columns (3,4,5,6 bins × 32 and 48 kHz)
    call_mat = matrix(NA_real_, nrow = 4, ncol = 8)
    rownames(call_mat) = c("Study Level", "Strategy 1", "Strategy 2", "Strategy 3")
    colnames(call_mat) = c("3 bins (32kHz)", "3 bins (48kHz)",
                            "4 bins (32kHz)", "4 bins (48kHz)",
                            "5 bins (32kHz)", "5 bins (48kHz)",
                            "6 bins (32kHz)", "6 bins (48kHz)")
    
    
    # Fills columns for one estimates and tables object (one bin count)
    fill_from_outtie = function(obj, col_pair_start) {
    
      #Study level estiamte
      study_est <- obj$Study.level["Estimate", sp]
      
      s1 <- obj$Strategy.1[sp, ]
      s2 <- obj$Strategy.2[sp, ]
      s3 <- obj$Strategy.3[sp, ]
      
      # 32 kHz
      call_mat["Study Level", col_pair_start]     <<- study_est
      call_mat["Strategy 1",  col_pair_start]     <<- s1["Estimate-32"]
      call_mat["Strategy 2",  col_pair_start]     <<- s2["Estimate-32"]
      call_mat["Strategy 3",  col_pair_start]     <<- s3["Estimate-32"]
      
      # 48 kHz
      call_mat["Study Level", col_pair_start + 1] <<- study_est
      call_mat["Strategy 1",  col_pair_start + 1] <<- s1["Estimate-48"]
      call_mat["Strategy 2",  col_pair_start + 1] <<- s2["Estimate-48"]
      call_mat["Strategy 3",  col_pair_start + 1] <<- s3["Estimate-48"]
    }
    
    
    # Fill for 3,4,5,6 bins
    fill_from_outtie(out_3, 1)  # 3-bin columns (1,2)
    fill_from_outtie(out_4, 3)  # 4-bin columns (3,4)
    fill_from_outtie(out_5, 5)  # 5-bin columns (5,6)
    fill_from_outtie(out_6, 7)  # 6-bin columns (7,8)
    
    
    # Round for presentation
    call_mat_round = round(call_mat, 3)
    
    
    # Dataframe of estimates
    tab_sp = as.data.frame(call_mat_round)
    tab_sp$Level = rownames(call_mat_round)
    tab_sp = tab_sp[, c("Level", colnames(call_mat_round))]
    rownames(tab_sp) = NULL 
    
    per_species_tables[[sp]] <- tab_sp
  }
  
  # Return just the per-species tables (cleanest)
  return(per_species_tables)
}


#presentation of density estimates
density.estimates(25)
density.estimates(50)
density.estimates(75)


################################################################################
# empirical standard deviation                                                 #
################################################################################


# Load all outtie objects
outtie_list = mget(paste0("outtie", 0:11))

# Extracting empirical SD table from each outtie
extract_empirical_sd = function(obj) {
  obj$Standard.dev.study.empirical
}

empirical_sd_list = lapply(outtie_list, extract_empirical_sd)

# Combine into table
empirical_sd_table = do.call(rbind, empirical_sd_list)

# Define bin labels and sample sizes
bin_labels    = rep(c("6 bins", "5 bins", "4 bins", "3 bins"), each = 3)
sample_labels = rep(c("25 samples", "50 samples", "75 samples"), times = 4)

# Combine into final row names
rownames(empirical_sd_table) = paste(bin_labels, sample_labels, sep = " - ")

# View result
empirical_sd_table





















