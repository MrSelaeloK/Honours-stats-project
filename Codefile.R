
#-------------------------------------------------------------------------------
library(readr)
library(tidyverse)

Data = read_csv("Abyssinian Nightjar validation results.csv")
Data = data.frame(Data)

#Data collection for use
#This is the loading up of the data
UtilsDataRSV::view_cols(Data)
Data|>select("commonName",
             "scientificName",
             "confidence",
             logit_conf,
             outcome)


#Truncation of the data to remove any entries with a confidence of less than 0.1 and allows for a focuses on higher quality recordings
Data = Data |> filter(confidence>0.1) #Retraing only the high quality records

#need to seperate the data into quantiles for binning
Data |> arrange(logit_conf)
Data = Data |> mutate(quantile_group = ntile(Data[,"logit_conf"],5))

#Creation of bins
#Actually it turns out that i binned the data wrong
Bin1 = Data |> filter(quantile_group == 1)
Bin2 = Data |> filter(quantile_group == 2)
Bin3 = Data |> filter(quantile_group == 3)
Bin4 = Data |> filter(quantile_group == 4)
Bin5 = Data |> filter(quantile_group == 5)

#note that the data is not separated into bins by species name

#getting the distribution of binned data
distribution.data = Data |> filter(outcome != 0) |> select(logit_conf)

logit_distribution = ggplot(distribution.data, aes(x = logit_conf, y = ..density..))+
  geom_histogram(binwidth = 0.4, fill = "steelblue", colour = "white")+
  geom_density(colour= "red")+
  labs(title = "Histogram + Density Curve of logit_conf", x = "logit_conf", y = "Density")

#----Perch approach-------------------------------------------------------------

#Chunk 1 - validation example
#time_stamp is not necessary 
#score is the confidence score
#is_position will turn 1 if target species, -1 if not and 0 if unknown
#bin number the validation was assigned to
#Proportion of total observation in this bin


create_validation_example <- function(filename, 
                                      timestamp_offset,
                                      score, is_pos,
                                      bin,
                                      bin_weight)
{
  list(
    filename = filename,
    timestamp_offset = timestamp_offset, #i do not think that i need the time_stamp
    score = score, #score is confidence score
    is_pos = is_pos,
    bin = bin,
    bin_weight = bin_weight
  )
}


#Chunk 2 - takes in a list item and then sections it to get different parts. 
#Allows for it to be added the csv file
#filename is the file name
#timestamp_offset, I dont think will be important
#score is the confidence score associated with the example
#is_pos will tell if its a positive_id = 1, negative_id = -1 and 0 for unsure
#bin is the bin number
#bin_weight is the bin weight within the data space 

validation_to_row <- function(example) {
  c(
    filename = example$filename,
    timestamp_offset = example$timestamp_offset,
    score = example$score,
    is_pos = example$is_pos,
    bin = example$bin,
    bin_weight = example$bin_weight
  )
}



#Chunck 3 - 
#df reads the csv file
#examples refers to validation examples: 
#it takes a list from the data frame and then iterates through the csv file to name the different elements of the csv file.
#Returns a list of validation examples
load_validation_log <- function(filepath) {
  df <- read.csv(filepath)
  examples <- lapply(1:nrow(df), function(i) {
    row <- df[i, ]
    create_validation_example(
      filename = as.character(row$filename),
      timestamp_offset = as.numeric(row$timestamp_offset),
      score = as.numeric(row$score),
      is_pos = as.integer(row$is_pos),
      bin = as.integer(row$bin),
      bin_weight = as.numeric(row$bin_weight)
    )
  })
  return(examples)
}


#Chunk 4 - Call density estimation
#infix operator, if a is not missing return a else return b. It is handling missing keys
`%||%` <- function(a, b) if (!is.null(a)) a else b 

#bin_pos will be turns bin_pos (positive) into a list
#bin_neg will turns bin_neg in each example into a list
#bin weights are also turned into a list
estimate_call_density <- function(examples, num_beta_samples = 10000, beta_prior = 0.1) {
  bin_pos <- list()
  bin_neg <- list()
  bin_weights <- list()
  
  # Collect positive/negative counts and weights
  #b is the bin id
  #bin_weights stores bin_weight at index key index b
  #The if: if bin_pos is a positive id, bin_pos at key index b is assigned to 1 if value not missing
  #So each bin will have a count of how many items reside within it
  for (ex in examples) {
    b <- as.character(ex$bin) 
    bin_weights[[b]] <- ex$bin_weight
    if (ex$is_pos == 1) {
      bin_pos[[b]] <- (bin_pos[[b]] %||% 0) + 1
    } else if (ex$is_pos == -1) {
      bin_neg[[b]] <- (bin_neg[[b]] %||% 0) + 1
    }
  }
  
  # Maximum Likelihood Estimate of density
  #density_ev, counter initially set to 0 
  #p is the positive count
  #n is the negative count
  #w is the bin weight at the name b
  #there is a small error added here
  #reference to page 6 top left, The one about the law of total probability
  density_ev <- 0
  for (b in names(bin_weights)) {
    p <- bin_pos[[b]] %||% 0
    n <- bin_neg[[b]] %||% 0
    w <- bin_weights[[b]]
    density_ev <- density_ev + w * (p / (p + n + 1e-6))
  }
  
  # Bayesian sampling using beta distributions
  # num_beta_samples passed as a parameter up top.
  # outer for loop: iterates through the num_beta_samples
  # total is a counter set to 0
  # Inner for loop: for each bin, p is set to positive, n is set to false classification and both are set to 0 if      there is nothing counted.
  # will need an explanation of the total.
  q_betas <- numeric(num_beta_samples)
  for (i in 1:num_beta_samples) {
    total <- 0
    for (b in names(bin_weights)) {
      p <- bin_pos[[b]] %||% 0
      n <- bin_neg[[b]] %||% 0
      w <- bin_weights[[b]]
      total <- total + w * rbeta(1, p + beta_prior, n + beta_prior)
    }
    q_betas[i] <- total
  }
  
  return(list(density_ev = density_ev, samples = q_betas))
}

# Chunk 5 Estimating ROC AUC 
# function takes in the validation examples
estimate_roc_auc <- function(examples) {
  bin_pos <- list()
  bin_neg <- list()
  bin_weights <- list()
  
  # Collect bin-level counts and weights, described in the same way as up top.
  for (ex in examples) {
    b <- as.character(ex$bin)
    bin_weights[[b]] <- ex$bin_weight
    if (ex$is_pos == 1) {
      bin_pos[[b]] <- (bin_pos[[b]] %||% 0) + 1
    } else if (ex$is_pos == -1) {
      bin_neg[[b]] <- (bin_neg[[b]] %||% 0) + 1
    }
  }
  
  # Bin-level positive and negative probabilities, also described in the same way as up top. Then for the denominator, a separate variable was made.  
  p_pos_bin <- list()
  p_neg_bin <- list()
  for (b in names(bin_weights)) {
    p <- bin_pos[[b]] %||% 0
    n <- bin_neg[[b]] %||% 0
    denom <- p + n + 1e-6
    p_pos_bin[[b]] <- p / denom
    p_neg_bin[[b]] <- n / denom
  }
  
  # Global P(+) estimate, density estimates. (P(+) is the probability of a positive ID)
  # The function in the sapply gives the contribution to the positive density
  density_ev <- sum(sapply(names(bin_weights), function(b) {
    bin_weights[[b]] * p_pos_bin[[b]]
  }))
  
  # P(bin | +) and P(bin | -) 
  # Calculates the probability that a positive id comes from a particular bin
  # The same is being done for a negative id
  p_bin_pos <- list()
  p_bin_neg <- list()
  for (b in names(bin_weights)) {
    p_bin_pos[[b]] <- bin_weights[[b]] * p_pos_bin[[b]] / density_ev
    p_bin_neg[[b]] <- bin_weights[[b]] * p_neg_bin[[b]] / (1.0 - density_ev)
  }
  
  # Sort bins
  # Bins are sorted in decreasing order
  bins <- sort(as.integer(names(bin_weights)))
  
  # Off-diagonal contributions
  #seq_along bins will generate a sequence corresponding to the bins
  #for mismatching bins, the frequentist probability of a positive validation is matched with the probability of a    frequentist negative in another bin. Hence the term off diagonal.
  roc_auc <- 0
  for (i in seq_along(bins)) {
    for (j in (i + 1):length(bins)) {
      if (j > length(bins)) next
      b_i <- as.character(bins[i])
      b_j <- as.character(bins[j])
      roc_auc <- roc_auc + (p_bin_pos[[b_j]] %||% 0) * (p_bin_neg[[b_i]] %||% 0)
    }
  }
  
  # In-bin contributions (diagonal terms)
  # looks at positive scores and negative scores and then returns them and if the score is missing, it is removed from the score list. 
  # The next line is satisfied if there are no positive/negative scores.
  #if they are there
  for (b in bins) {
    b <- as.character(b)
    pos_scores <- sapply(examples, function(e) if (e$bin == as.integer(b) && e$is_pos == 1) e$score else NA_real_)
    neg_scores <- sapply(examples, function(e) if (e$bin == as.integer(b) && e$is_pos == -1) e$score else NA_real_)
    pos_scores <- na.omit(pos_scores)
    neg_scores <- na.omit(neg_scores)
    
    if (length(pos_scores) == 0 || length(neg_scores) == 0) next
    
    # Count how many times a positive score is greater than a negative
    # Hits is a pairwise check that +ve > -ve. 
    # bin_roc_auc is the probability that the pos_score are greater than the negative ones
    # roc_auc is how they are then added in the end
    hits <- sum(outer(pos_scores, neg_scores, FUN = "-") > 0)
    bin_roc_auc <- hits / (length(pos_scores) * length(neg_scores))
    roc_auc <- roc_auc + bin_roc_auc * (p_bin_pos[[b]] %||% 0) * (p_bin_neg[[b]] %||% 0)
  }
  
  return(roc_auc)
}


#Writing validation in log R - This part of the code works to eliminate duplicates. When a line is written to it with the same time stamp and it is duplicated, there is a replacement that happens.
write_validation_log <- function(validation_examples, output_path, target_class) {
  filename <- paste0("validation_", target_class, ".csv")
  validation_log_filepath <- file.path(output_path, filename)
  
  # Deduplicate by (filename, timestamp_offset)
  #examples_map creates a new environment
  examples_map <- new.env(hash = TRUE, parent = emptyenv())
  
  #This stores the keys in the hash map for each item.
  if (file.exists(validation_log_filepath)) {
    existing <- load_validation_log(validation_log_filepath)
    for (ex in existing) {
      key <- paste0(ex$filename, ":", ex$timestamp_offset)
      examples_map[[key]] <- ex
    }
  }
  
  for (ex in validation_examples) {
    key <- paste0(ex$filename, ":", ex$timestamp_offset)
    examples_map[[key]] <- ex
  }
  
  final_examples <- as.list(examples_map)
  
  # Convert to data frame and write
  rows <- lapply(final_examples, validation_to_row)
  df <- do.call(rbind, rows)
  write.csv(df, file = validation_log_filepath, row.names = FALSE)
  
  return(validation_log_filepath)
}

#Prune random results in R- just for the size constraining and keeping only of samples per bin
# if a bin has 

prune_random_results <- function(results,
                                 all_scores,
                                 quantile_bounds,
                                 samples_per_bin) {
  
  value_bounds <- quantile(all_scores, quantile_bounds)
  num_bins <- length(quantile_bounds) - 1
  binned <- vector("list", num_bins)
  
  # Assign each result to a bin
  for (res in results) {
    bin_idx <- max(which(res$score < value_bounds)) - 1
    bin_idx <- max(0, bin_idx)
    binned[[bin_idx + 1]] <- c(binned[[bin_idx + 1]], list(res))
  }
  
  # Sample from each bin
  sampled_bins <- lapply(binned, function(bin) {
    if (length(bin) > samples_per_bin) {
      sample(bin, samples_per_bin)
    } else {
      bin
    }
  })
  
  combined_results <- unlist(sampled_bins, recursive = FALSE)
  combined_results <- sample(combined_results)  # Shuffle
  return(combined_results)
}

#This is the part i need to explain
#convert_combined_results() in R
convert_combined_results <- function(combined_results,
                                     target_class,
                                     quantile_bounds,
                                     value_bounds) {
  examples <- list()
  for (r in combined_results) {
    ex <- from_search_result(r, target_class, quantile_bounds, value_bounds)
    if (!is.null(ex)) {
      examples <- c(examples, list(ex))
    }
  }
  return(examples)
}

#get_random_sample_size() in R
get_random_sample_size <- function(quantile_bounds, samples_per_bin) {
  bin_weights <- diff(quantile_bounds)
  rarest_weight <- min(bin_weights)
  return(ceiling(2 * samples_per_bin / rarest_weight))
}


# ----Making use of the code on data--------------------------------------------





