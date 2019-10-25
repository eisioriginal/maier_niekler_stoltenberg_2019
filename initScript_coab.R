#####################################################
## run topic models #################################
#####################################################

source("script_dm.r")

#--------------------------------------------------
# new model
#--------------------------------------------------
makeSimpleTripletMatrix <- function(sparseDTMMatrix){
  tmp <- Matrix::summary(sparseDTMMatrix)
  
  return(slam::simple_triplet_matrix(i=tmp[,1], j=tmp[,2], v=tmp[,3],
                                     dimnames=dimnames(sparseDTMMatrix),nrow = nrow(sparseDTMMatrix), ncol = ncol(sparseDTMMatrix)))
}

#load dtm
#Download from server first
load(file="dtm_coab.Rdata")

# path to badtokens
blacklistBadTokens <- readLines("./backlist_badtokens.txt", encoding = "UTF-8")

intersect(colnames(dtm),tmca.util::preprocess(blacklistBadTokens))

#tokensToIgnore <- tmca.util::preprocess(blacklistBadTokens)
tokensToIgnore <-blacklistBadTokens %>%
  tokens(remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE)  %>%  tokens_tolower() 
# path to storage Directory
storageDirectory <- paste0("./res_corpora/models_coab/CorpusSampling")
sampleName <- "online_food_us"

# Komplettierung Fixed Sample/Fixed Ini
# Set LDA parameters
# Number of topics
K <- 50
# Iterations should be at least 1000 for final evaluation
iterations <- 1000
nCoocs <- 5000
alphaPriorsToTest <- c(0.5)
#alphaPriorsToTest <- c(0.5)
# Vector of values for term-topic distribution priors to test
etaPriorsToTest <- c(1 / K)
runs <- 1
fraction_sizes <- c(0.01, 0.05, 0.1, 0.2, 0.5)
#fraction_sizes <- c(0.5)

library(tictoc)
results <- list()

# Sample Models 5 runs, Random Initialization, for 1%,5%,10%, 20%, 50% sample size
for (frac in fraction_sizes){
  tic(paste(frac))
  for (i in 1:runs){
    run_number <- i
    print(paste0("run: ", run_number, ", frac: ", frac, " IniMeth: rnd"))
    modelEvaluationData <- runLDAModelEvalution(dtm, run_number, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                                blacklist = tokensToIgnore,  
                                                initTopicAssignments = F, 
                                                clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = nCoocs, fraction = frac, fixedSample = TRUE)
  }
  results[[length(results)+1]] <- toc()
}

#Reference Model 5 runs, Random Initialization

for (i in 1:runs){
  run_number <- i
  print(paste0("run: ", run_number, " Reference Model", " IniMeth: rnd"))
  
  tic(paste("1"))
  
  modelEvaluationData <- runLDAModelEvalution(dtm, run_number, K, iterations, alphaPriorsToTest, etaPriorsToTest, 
                                              blacklist = tokensToIgnore, 
                                              initTopicAssignments = F, 
                                              clusterMethod = "PAM", cooccurrenceMeasure = "LL", nCoocs = nCoocs, fraction = NULL, fixedSample = TRUE)
  results[[length(results)+1]] <- toc()
}