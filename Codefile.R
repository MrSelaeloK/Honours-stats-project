
#-------------------------------------------------------------------------------
library(readr)
library(tidyverse)

Data = read_csv("Abyssinian Nightjar validation results.csv")
Data = data.frame(Data)
UtilsDataRSV::view_cols(Data)

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

#debug(generate_fake_validation_data)
fdata = generate_fake_validation_data(n = 100,
                                      bins= 5)

#---------------- Binning the data----------------------------------------------
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



#--------Bin weights: already redifined it so this may not be necessary--------- 
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



# -------Chunk 2 - Call density estimation-------------------------------------
#the algorithm will work as long it is supplied with data that is not binned and weighted but is validated binning and weighting happen within the function.
call_density = function(examples,num_beta_samples=1000,beta_prior=0.1){
  #examples = Data2
  examples = unlist(examples)
  df <- as.data.frame(do.call(rbind, fdata), stringsAsFactors = FALSE) #This will be the one that i check my bin and bin_weight generation against. 
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
  density = 0 
  for(i in 1:nrow(shape.parameters)){
    density = density + rbeta(1,shape.parameters[i,1]+beta_prior,
                             shape.parameters[i,2]+beta_prior)*bin_weights[i,1]
  }

    
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
              Samples = q_betas))
}

#Implementation of call density  
out = call_density(fdata,1000,0.1)  
out$Density_Estimate
out$Samples # still wrapping my heat around what the samples are supposed to be.
hist(out$Samples)






