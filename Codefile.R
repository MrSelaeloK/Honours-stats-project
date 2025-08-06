
#-------------------------------------------------------------------------------
library(readr)
library(tidyverse)

Data = read_csv("Abyssinian Nightjar validation results.csv")
Data = data.frame(Data)
UtilsDataRSV::view_cols(Data)

#----Perch approach-------------------------------------------------------------
# Data binning
# Will have to code for the binning of the data here
# One to one mapping between the confidence and the logit so i will map it out that way
# Data$confidence
bin_column = matrix(NA,nrow = nrow(Data), ncol = 1)
bins = function(unbinned_data){
  for(i in 1 : nrow(unbinned_data)){
    if(Data$confidence[i]<=0.5){
      bin_column[i]=1
      next
    }
    if(Data$confidence[i]>0.5 && Data$confidence[i]<=0.75){
      bin_column[i]=2
      next
    }
    if(Data$confidence[i]>0.75 && Data$confidence[i]<=0.875){
      bin_column[i]=3
      next
    }
    else{
      bin_column[i]=4
    }
  }
  return(bin_column)
}
bins(Data)

Data1 = mutate(Data,bins(Data))    #Binned data
names(Data1)[11] <- "binNumber"    #Renaming of an item in the bin
UtilsDataRSV::view_cols(Data1)



#Bin weights
binWeight = function(binned_data){
  
  numberOfBins = count(distinct(Data1, binNumber))
  bin1Weight = count(filter(Data1,binNumber==1))/nrow(Data1)
  bin2Weight = count(filter(Data1,binNumber==2))/nrow(Data1)
  bin3Weight = count(filter(Data1,binNumber==3))/nrow(Data1)
  bin4Weight = count(filter(Data1,binNumber==4))/nrow(Data1)
  
  weights = matrix(NA, nrow=nrow(binned_data), ncol = 1)
  for(i in 1 : nrow(binned_data)){
    if(binned_data$binNumber[i]==1){
      weights[i]=bin1Weight
      next
    }
    if(binned_data$binNumber[i]==2){
      weights[i]=bin2Weight
      next
    }
    if(binned_data$binNumber[i]==3){
      weights[i]=bin3Weight
      next
    }
    else{
      weights[i]=bin4Weight
    }
  }
  
  weights = as.numeric(weights)
  Data2 = mutate(binned_data,weights = weights)
  return(Data2)
}

Data2 = binWeight(Data1)
UtilsDataRSV::view_cols(Data2)

#so the data that will be used for density estimation will be called Data2



# Chunk 4 - Call density estimation
# infix operator, if a is not missing return a else return b. It is handling missing keys
`%||%` <- function(a, b) if (!is.null(a)) a else b 

# bin_pos will be turns bin_pos (positive) into a list
# bin_neg will turns bin_neg in each example into a list
# bin weights are also turned into a list
estimate_call_density <- function(examples, num_beta_samples = 10000 , beta_prior = 0.1) {
  
  # infix operator, if a is not missing return a else return b. It is handling missing keys
  `%||%` <- function(a, b) if (!is.null(a)) a else b 
  
  bin_pos <- list()
  bin_neg <- list()
  bin_weights <- list()
  #bin_unknown = list()
  
  # Collect positive/negative counts and weights
  #b is the bin id
  #bin_weights stores bin_weight at index key index b
  #The if: if bin_pos is a positive id, bin_pos at key index b is assigned to 1 if value not missing
  #So each bin will have a count of how many items reside within it
examples = Data2

#Collect positive/negative counts and weights
for (i in 1:nrow(examples)) {
  b <- as.character(examples$binNumber[i])  # convert to character for list indexing
  w <- examples$weights[i]
  
  bin_weights[[b]] <- w
  
  if (examples$outcome[i] == 1) {
    bin_pos[[b]] <- (bin_pos[[b]] %||% 0) + 1
  } else if (examples$outcome[i] == -1) {
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

call.density = estimate_call_density(Data2,10000,0.1)

beta(59,0.001)


#---------call density estimate re-estimation-----------------------------------
call_density = function(examples,num_beta_samples=1000,beta_prior=0.1){
  #examples = Data2
  
  #infinix operator
  `%||%` <- function(a, b) if (!is.null(a)) a else b 
  
  #Positive/negative bin counts
  bin_positive = matrix(NA,nrow = 4, ncol = 2)
  bin_negative = matrix(NA,nrow = 4, ncol = 2)
  bin_number  = as.matrix(distinct(examples,binNumber))
  bin_weights = matrix(NA, nrow = nrow(bin_number),ncol=1)
  
  #Bin weights
  for(i in 1:nrow(bin_number)){
    bin_weights[i,1]=nrow(filter(examples,binNumber==i))/nrow(examples)
  }
  
  
  #Putting in the counts of each bin (positive/negative)
  # i do not know
  for(i in 1:nrow(bin_number)){
    examples2 = filter(examples,binNumber==i)
    bin_positive[i,2] = as.numeric(count(filter(examples2,outcome==1)))
    bin_negative[i,2] = as.numeric(count(filter(examples2,outcome==-1)))
  }
  bin_positive[,1]=bin_number ; bin_negative[,1] = bin_number
  

#Density  
  #will have to use positive beta weights
  # Create beta distributions.
  shape.parameters = matrix(NA,ncol=2,nrow=nrow(bin_number))
  
  for (b in bin_number){
    shape.parameters[b,] = c(
      bin_positive[b,2] + beta_prior, bin_negative[b,2] + beta_prior
    )
  }
  
#MLE of positive rate in each bin
  density_eval = 0
    for(j in 1:nrow(bin_number)){
    density_eval = density_eval = bin_weights[j] * bin_positive[j,2] /
      (bin_positive[j,2] + bin_negative[j,2] + 1e-6)
  }
    
#probability of a positive outcome. Supposed to be the probability of bin b given outcome
  q_betas=c()
  NoBetaSamples = 10000
  for(i in 1:NoBetaSamples){
    for(j in 1:nrow(bin_number)){
      rvalue = 0
      rvalue = rvalue + rbeta(1,shape.parameters[j,1],shape.parameters[j,2])
    }
    q_betas=c(q_betas,rvalue)
  }
  return(list(Density_Estimate = density_eval,
              Samples = q_betas))
}
  
out = call_density(Data2,1000,0.1)  
  
out$Density_Estimate
out$Samples # still wrapping my heat around what the samples are supposed to be.
  












#For low sample sizes 
Booststrapping = function(){
  
}




