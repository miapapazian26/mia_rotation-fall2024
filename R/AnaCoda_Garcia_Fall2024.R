#### Fall 2024 - Rotation w/ Dr. Gilchrist ###

#Objective 1: Learn R Package 'AnaCoDa'

#set working directory
setwd("C:/Users/nyanc/OneDrive/Documents/All Files - personal/School_UT/Fall2024/Rotation_Gilchrist")

#download r package - dowloaded .tar.gz from CRAN, had to also install VGAM first
#install.packages("~/All Files - personal/School_UT/Fall2024/Rotation_Gilchrist/AnaCoDa_0.1.4.4.tar.gz", 
#                 repos = NULL, type = "source")

    #When installing, it mentions this packages most likely requires manual configuration via script

#Load libraries
library(tidyverse)
library(AnaCoDa)
library(seqinr)
#install.packages("bayesplot")
library(bayesplot)
library(data.table)

#Define function
chomp<-function(dir){  
  
  genome <- read_delim(dir, delim = "\n", col_names = FALSE, col_types = "c")
  genome.list <- genome %>% na.omit() %>% data.frame()
  
  i=0                                             #Initialize empty objects for faster processing time
  flag=0
  y<-data.frame()
  gene.length<-list()
  sequence.all<-list()
  alist<-list()
  #species.name<-species.list                      #this will be user supplied, must match how it looks in FASTA
  b<-''
  c<-''
  
  cat(paste("Now importing:",dir, "\n"))        #basic progress checker
  for(i in genome.list[ ,1]){                 #for each line of inputted genomeFASTA
    if(str_detect(i,"^>")){                     #if row starts with >
      if(flag==1){                              #flag check to concatenate sequence lines
        sequence<-str_c(alist,collapse = "")    #collapse
        aasize<-nchar(sequence)                 #get amino acid size
        gene.length<-rbind.data.frame(gene.length,aasize)  #store size
        sequence.all<-rbind.data.frame(sequence.all,sequence)  #store sequence
        alist<-list()                           #empty temp variable
        flag=0                                  #reset flag
      }
      flag=0                                    #reset
      pattern.chr<-paste0("Chromosome:\\S+")#make species name pattern
      pattern.gn<-paste0("description:.*")             #make gene name pattern
      pattern.type<-paste0("gene_biotype:\\S+")            #make >tr or >sp pattern 
      a<-str_extract(i,pattern = pattern.chr) #find and store species name (user supplies species name)
      b<-str_extract(i,pattern = pattern.gn)    #find and store gene name
      c<-str_extract(i,pattern = pattern.type)  #find and store pattern type (>tr or >sp, should be only two options)
      abci<-c(a,b,c,i)                             #make above identifiers into easy to bind object
      y<-rbind(y,abci)                           #rbind to fill rows of new df
    }
    else{
      flag=1                                    #Raise the flag!
      b<-str_trim(i,side = c("right"))          #Cut off new line character
      alist[[i]]<-b                             #Storing all sequence lines of current protein
    }
  }
  if(flag==1){                                  #this catches the last protein's sequence VVV
    sequence<-str_c(alist,collapse = "")        #
    aasize<-nchar(sequence)                     #
    gene.length<-rbind.data.frame(gene.length,aasize) #
    sequence.all<-rbind.data.frame(sequence.all,sequence) #
    alist<-list()                               #
    flag=0                                      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  }
  y[,5]<-gene.length                             #add sizes to df
  y[,6]<-sequence.all                             #add sequences to df
  colnames(y)<-c("Location","Description","Type","Header","Length","Sequence") #add colnames to object 'y'
  y$Description <- str_replace_all(y$Description, "description:", "")  #remove OS= from species column
  y$Location <- str_replace_all(y$Location, "Chromosome:", "")  #remove GN= from gene name column
  y$Type <- str_replace_all(y$Type, "gene_biotype:", "")
  y$Header <- str_replace_all(y$Header, ">", "")
  
  #Count the number of 'genes' that don't have descriptions
  perc <- (sum(is.na(y$Description)))/nrow(y)*100
  cat("Percentage of individual sequences w/o description:", perc)
  
  genome <<- y
}   #Load 'chomp' function to import .fasta files
loadROCParameterObject <- function(parameter, files){
  
  setBaseInfo <- function(parameter, files){
    for (i in 1:length(files)) {
      tempEnv <- new.env();
      load(file = files[i], envir = tempEnv)
      if (i == 1) {
        categories <- tempEnv$paramBase$categories
        categories.matrix <- do.call("rbind", tempEnv$paramBase$categories)
        numMixtures <- tempEnv$paramBase$numMix
        numMutationCategories <- tempEnv$paramBase$numMut
        numSelectionCategories <- tempEnv$paramBase$numSel
        mixtureAssignment <- tempEnv$paramBase$curMixAssignment
        lastIteration <- tempEnv$paramBase$lastIteration
        max <- tempEnv$paramBase$lastIteration + 1
        grouplist <- tempEnv$paramBase$grouplist
        
        stdDevSynthesisRateTraces <- vector("list", length = numSelectionCategories)
        for (j in 1:numSelectionCategories) {
          stdDevSynthesisRateTraces[[j]] <- tempEnv$paramBase$stdDevSynthesisRateTraces[[j]][1:max]
        }
        stdDevSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$stdDevSynthesisRateAcceptRatTrace
        synthesisRateTrace <- vector("list", length = numSelectionCategories)
        for (j in 1:numSelectionCategories) {
          for (k in 1:length(tempEnv$paramBase$synthRateTrace[[j]])){
            synthesisRateTrace[[j]][[k]] <- tempEnv$paramBase$synthRateTrace[[j]][[k]][1:max]
          }
        }
        synthesisRateAcceptanceRateTrace <- tempEnv$paramBase$synthAcceptRatTrace
        mixtureAssignmentTrace <- vector("list", length = length(tempEnv$paramBase$mixAssignTrace))
        for (j in 1:length(tempEnv$paramBase$mixAssignTrace)){
          mixtureAssignmentTrace[[j]] <- tempEnv$paramBase$mixAssignTrace[[j]][1:max]
        }
        mixtureProbabilitiesTrace <- c()
        for (j in 1:numMixtures) {
          mixtureProbabilitiesTrace[[j]] <- tempEnv$paramBase$mixProbTrace[[j]][1:max]
        }
        codonSpecificAcceptanceRateTrace <- tempEnv$paramBase$codonSpecificAcceptRatTrace
        
        ### ERROR HERE ###
        withPhi <- tempEnv$paramBase$withPhi
        if (withPhi){
          phiGroups <- length(tempEnv$paramBase$synthesisOffsetTrace) #add $paramBase
          synthesisOffsetTrace <- c()
          for (j in 1:phiGroups) {
            synthesisOffsetTrace[[j]] <- tempEnv$paramBase$synthesisOffsetTrace[[j]][1:max]
          }
          
          
          synthesisOffsetAcceptanceRateTrace <- tempEnv$paramBase$synthesisOffsetAcceptRatTrace
          
          
          observedSynthesisNoiseTrace <- c()
          for (j in 1:phiGroups) {
            observedSynthesisNoiseTrace[[j]] <- tempEnv$paramBase$observedSynthesisNoiseTrace[[j]][1:max]
          }
          #need number of phi groups, not the number of mixtures apparently.
        }else {
          synthesisOffsetTrace <- c()
          synthesisOffsetAcceptanceRateTrace <- c()
          observedSynthesisNoiseTrace <- c()
        }
      } else {
        if (sum(categories.matrix != do.call("rbind", tempEnv$paramBase$categories)) != 0){
          stop("categories is not the same between all files")
        }#end of error check
        
        if (numMixtures != tempEnv$paramBase$numMix){
          stop("The number of mixtures is not the same between files")
        }
        
        if (numMutationCategories != tempEnv$paramBase$numMut){
          stop("The number of mutation categories is not the same between files")
        }
        
        if (numSelectionCategories != tempEnv$paramBase$numSel){
          stop("The number of selection categories is not the same between files")
        }
        
        if (length(mixtureAssignment) != length(tempEnv$paramBase$curMixAssignment)){
          stop("The length of the mixture assignment is not the same between files.
               Make sure the same genome is used on each run.")
        }
        
        if(length(grouplist) != length(tempEnv$paramBase$grouplist)){
          stop("Number of Amino Acids/Codons is not the same between files.")
        }
        if (withPhi != tempEnv$paramBase$withPhi){
          stop("Runs do not match in concern in with.phi")
        }
        
        curSynthesisOffsetTrace <- tempEnv$paramBase$synthesisOffsetTrace
        curSynthesisOffsetAcceptanceRateTrace <- tempEnv$paramBase$synthesisOffsetAcceptRatTrace
        curObservedSynthesisNoiseTrace <- tempEnv$paramBase$observedSynthesisNoiseTrace
        
        if (withPhi){
          combineTwoDimensionalTrace(synthesisOffsetTrace, curSynthesisOffsetTrace, max)
          size <- length(curSynthesisOffsetAcceptanceRateTrace)
          combineTwoDimensionalTrace(synthesisOffsetAcceptanceRateTrace, curSynthesisOffsetAcceptanceRateTrace, size)
          combineTwoDimensionalTrace(observedSynthesisNoiseTrace, curObservedSynthesisNoiseTrace, max)
        }
        
        
        curStdDevSynthesisRateTraces <- tempEnv$paramBase$stdDevSynthesisRateTraces
        curStdDevSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$stdDevSynthesisRateAcceptRatTrace
        curSynthesisRateTrace <- tempEnv$paramBase$synthRateTrace
        curSynthesisRateAcceptanceRateTrace <- tempEnv$paramBase$synthAcceptRatTrace
        curMixtureAssignmentTrace <- tempEnv$paramBase$mixAssignTrace
        curMixtureProbabilitiesTrace <- tempEnv$paramBase$mixProbTrace
        curCodonSpecificAcceptanceRateTrace <- tempEnv$paramBase$codonSpecificAcceptRatTrace
        
        lastIteration <- lastIteration + tempEnv$paramBase$lastIteration
        
        
        #assuming all checks have passed, time to concatenate traces
        max <- tempEnv$paramBase$lastIteration + 1
        combineTwoDimensionalTrace(stdDevSynthesisRateTraces, curStdDevSynthesisRateTraces, max)
        
        size <- length(curStdDevSynthesisRateAcceptanceRateTrace)
        stdDevSynthesisRateAcceptanceRateTrace <- c(stdDevSynthesisRateAcceptanceRateTrace,
                                                    curStdDevSynthesisRateAcceptanceRateTrace[2:size])
        
        
        combineThreeDimensionalTrace(synthesisRateTrace, curSynthesisRateTrace, max)
        size <- length(curSynthesisRateAcceptanceRateTrace)
        combineThreeDimensionalTrace(synthesisRateAcceptanceRateTrace, curSynthesisRateAcceptanceRateTrace, size)
        
        combineTwoDimensionalTrace(mixtureAssignmentTrace, curMixtureAssignmentTrace, max)
        combineTwoDimensionalTrace(mixtureProbabilitiesTrace, curMixtureProbabilitiesTrace, max)
        size <- length(curCodonSpecificAcceptanceRateTrace)
        combineTwoDimensionalTrace(codonSpecificAcceptanceRateTrace, curCodonSpecificAcceptanceRateTrace, size)
      }
    }
    
    parameter$setCategories(categories)
    parameter$setCategoriesForTrace()
    parameter$numMixtures <- numMixtures
    parameter$numMutationCategories <- numMutationCategories
    parameter$numSelectionCategories <- numSelectionCategories
    parameter$setMixtureAssignment(tempEnv$paramBase$curMixAssignment) #want the last in the file sequence
    parameter$setLastIteration(lastIteration)
    parameter$setGroupList(grouplist)
    
    trace <- parameter$getTraceObject()
    trace$setStdDevSynthesisRateTraces(stdDevSynthesisRateTraces)
    trace$setStdDevSynthesisRateAcceptanceRateTrace(stdDevSynthesisRateAcceptanceRateTrace)
    trace$setSynthesisRateTrace(synthesisRateTrace)
    trace$setSynthesisRateAcceptanceRateTrace(synthesisRateAcceptanceRateTrace)
    trace$setSynthesisOffsetTrace(synthesisOffsetTrace)
    trace$setSynthesisOffsetAcceptanceRateTrace(synthesisOffsetAcceptanceRateTrace)
    trace$setObservedSynthesisNoiseTrace(observedSynthesisNoiseTrace)
    trace$setMixtureAssignmentTrace(mixtureAssignmentTrace)
    trace$setMixtureProbabilitiesTrace(mixtureProbabilitiesTrace)
    trace$setCodonSpecificAcceptanceRateTrace(codonSpecificAcceptanceRateTrace)
    
    parameter$setTraceObject(trace)
    return(parameter)
  } #changed a single line
  
  parameter <- setBaseInfo(parameter, files)
  for (i in 1:length(files)){
    tempEnv <- new.env();
    load(file = files[i], envir = tempEnv)
    
    numMutationCategories <- tempEnv$paramBase$numMut
    numSelectionCategories <- tempEnv$paramBase$numSel
    max <- tempEnv$paramBase$lastIteration + 1
    
    if (i == 1){
      
      
      codonSpecificParameterTraceMut <- vector("list", length=numMutationCategories)
      for (j in 1:numMutationCategories) {
        codonSpecificParameterTraceMut[[j]] <- vector("list", length=length(tempEnv$mutationTrace[[j]]))
        for (k in 1:length(tempEnv$mutationTrace[[j]])){
          codonSpecificParameterTraceMut[[j]][[k]] <- tempEnv$mutationTrace[[j]][[k]][1:max]
        }
      }
      
      codonSpecificParameterTraceSel <- vector("list", length=numSelectionCategories)
      for (j in 1:numSelectionCategories) {
        codonSpecificParameterTraceSel[[j]] <- vector("list", length=length(tempEnv$selectionTrace[[j]]))
        for (k in 1:length(tempEnv$selectionTrace[[j]])){
          codonSpecificParameterTraceSel[[j]][[k]] <- tempEnv$selectionTrace[[j]][[k]][1:max]
        }
      }
    }else{
      
      curCodonSpecificParameterTraceMut <- tempEnv$mutationTrace
      curCodonSpecificParameterTraceSel <- tempEnv$selectionTrace
      combineThreeDimensionalTrace(codonSpecificParameterTraceMut, curCodonSpecificParameterTraceMut, max)
      combineThreeDimensionalTrace(codonSpecificParameterTraceSel, curCodonSpecificParameterTraceSel, max)
    }#end of if-else
  }#end of for loop (files)
  
  trace <- parameter$getTraceObject()
  
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceMut, 0)
  trace$setCodonSpecificParameterTrace(codonSpecificParameterTraceSel, 1)
  
  parameter$currentMutationParameter <- tempEnv$currentMutation
  parameter$currentSelectionParameter <- tempEnv$currentSelection
  parameter$proposedMutationParameter <- tempEnv$proposedMutation
  parameter$proposedSelectionParameter <- tempEnv$proposedSelection
  parameter$setTraceObject(trace)
  return(parameter)
}
load_file <- function(file){
  as.data.frame(fread(file, stringsAsFactors = FALSE))
}  

# #Following R script provided in documentation VVV
# 
# #Initialize Genome File of Interest (CDS) - Ex. E.coli (Escherichia coli str. K-12 substr. W3110 (GCA_000010245))
# #Downloaded 9/8/2424 (CDS, AA): 
# #     https://bacteria.ensembl.org/Escherichia_coli_str_k_12_substr_w3110_gca_000010245/Info/Index
# genome <- initializeGenomeObject(file = "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa")
# 
# #Initialize Parameter Object (Beware of typo in original code)
# parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1, 
#                                        gene.assignment = rep(1, length(genome)))
# 
# 
# model <- initializeModelObject(parameter = parameter, model = "ROC")
# 
# 
# mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)
# 
# 
# runMCMC(mcmc = mcmc, genome = genome, model = model)


#Storing Output
#...not sure what the output should look like, so had to stop the script for now...

#Restarting the script, letting it run for a bit. 

# Adjustment range: < 0.225 or > 0.325 
# AA	Acc.Rat
# A:	0.16
# C:	0.31
# D:	0.274
# E:	0.318
# F:	0.296
# G:	0.116
# H:	0.302
# I:	0.214
# K:	0.384
# L:	0.122
# N:	0.292
# P:	0.188
# Q:	0.308
# R:	0.152
# S:	0.164
# T:	0.12
# V:	0.256
# Y:	0.348
# Z:	0.292
# Acceptance rate for synthesis rate:
#   Target range: 0.225-0.325 
# Adjustment range: < 0.225 or > 0.325 
# acceptance rates below lower target of 0.225: 163
# acceptance rate above upper target of 0.325: 180
# ##################################################
# Geweke Score after 25000 iterations: 0.280067
# ##################################################
# Stopping run based on convergence after 25000 iterations

# Acceptance rates for Codon Specific Parameters
# Target range: 0.175-0.375 
# Adjustment range: < 0.225 or > 0.325 
# AA	Acc.Rat
# A:	0.236
# C:	0.276
# D:	0.322
# E:	0.342
# F:	0.25
# G:	0.172
# H:	0.308
# I:	0.216
# K:	0.338
# L:	0.134
# N:	0.242
# P:	0.156
# Q:	0.288
# R:	0.132
# S:	0.154
# T:	0.2
# V:	0.184
# Y:	0.302
# Z:	0.296
# Acceptance rate for synthesis rate:
#   Target range: 0.225-0.325 
# Adjustment range: < 0.225 or > 0.325 
# acceptance rates below lower target of 0.225: 141
# acceptance rate above upper target of 0.325: 150
# ##################################################
# Geweke Score after 50000 iterations: -27.5509
# ##################################################
# leaving MCMC loop

# #The script finished. I think it took around an hour
# 
# #Saving data using C++ methods (which won't work with R methods)
# writeParameterObject(parameter = parameter, file = "parameter_out.Rda")
# writeMCMCObject(mcmc = mcmc, file = "mcmc_out.Rda")
# 
# # #Use this code to load it back in
# # parameter <- loadParameterObject(file = "parameter_out.Rda")
# # mcmc <- loadMCMCObject(file = "mcmc_out.Rda")
# 
# #Now to view results
# csp_mat <- getCSPEstimates(parameter = parameter)   #, CSP="Mutation", mixture = 1, samples = 1000)
# 
# head(csp_mat)
# 
# #posterior estimates for gene specific parameters (all genes)
# phi_mat <- getExpressionEstimates(parameter = parameter, 
#                                   gene.index = 1:length(genome),
#                                   samples = 1000)
# head(phi_mat)
# 
# #posterior estimates for gene specific parameters (sample 100 random)
# phi_mat.random <- getExpressionEstimates(parameter = parameter, 
#                                   gene.index = sample(1:length(genome), 100), 
#                                   samples = 1000)
# head(phi_mat.random)
# 
# #Calculate selection coefficient s for each codon and each gene
# selection.coefficients <- getSelectionCoefficients(genome = genome, 
#                                                    parameter = parameter, 
#                                                    samples = 1000)
# head(selection.coefficients)
# 
# #Now to compare our values to the weights from CAI
# cai.weights <- getCAIweights(referenceGenome = genome)
# head(cai.weights)
# 
#     # GCA       GCC       GCG       GCT       TGC       TGT 
#     # 0.5964128 0.7599486 1.0000000 0.4506385 1.0000000 0.8002539 
# 
# nc.per.aa <- getNcAA(genome = genome)
# head(nc.per.aa)
# 
# #Comparing distributino of selection coefficients to CAI values estimated from ref gene set
# selection.coefficients <- getSelectionCoefficients(genome = genome, 
#                                                    parameter = parameter, 
#                                                    samples = 1000)
# s <- exp(selection.coefficients)
# cai.weights <- getCAIweights(referenceGenome = genome)
# 
# codon.names <- colnames(s)
# h <- hist(s[, 1], plot = F)
# plot(NULL, NULL, axes = F, xlim = c(0,1), ylim = range(c(0,h$counts)), 
#      xlab = "s", ylab = "Frequency", main = codon.names[1], cex.lab = 1.2)
# lines(x = h$breaks, y = c(0,h$counts), type = "S", lwd=2)
# abline(v = cai.weights[1], lwd=2, lty=2)
# axis(1, lwd = 3, cex.axis = 1.2)
# axis(2, lwd = 3, cex.axis = 1.2)
# 
# #Now for some diagnostics
# trace <- getTrace(parameter)
# plot.diag <- plot(x = trace, what = "Mutation", mixture = 1) #make figure panel larger to acc omdate graph
# 
# #visualize the results of the model fit
# #Ex. use the last 500 samples from mixture 1 for posterior estimate.
# plot(x = model, genome = genome, samples = 500, mixture = 1)
# 
# #This will compare different 'mixtures' i.e. gene sets
# plot(parameter, what = "Selection", samples = 500)

### Start here when loading 'post_initial_run.Rdata'
#Now that we can go from genome file to plots, next step is partitioning the E.coli genome into
# groups to compare against (i.e. define multiple mixtures) (Ex. lagging strand genes vs. leading strand genes)

#Idea:
# When downloading the genome (CDS) from Ensembl, each gene SHOULD have either a 1 or -1 assigned to them
# 1: forward chromosomal strand 
# -1: reverse chromosomal strand
#
# It should be simple enough to partition the genes into two groups based on forward or reverse strand.

#Method: Use modified chomp function to import the .fasta and pull chromosome identifier from .fasta header
#Idea: Rbind every gene makes it slower, initialize all lists to fit data, then populate
# #Testing an option to speed up chomp
# genome <- read_delim(dir, delim = "\n", col_names = FALSE, col_types = "c")
# 
# #This counts the number of fasta headers (starts with >) which should equate to gene number
# temp <- genome %>%
#   filter(str_detect(X1, "^>")) %>% nrow()

# ##### Start Here! #####
# 
# #Define .fasta file in working directory
# dir <- "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"
# 
# #Define function
# chomp<-function(dir){  
#   
#   genome <- read_delim(dir, delim = "\n", col_names = FALSE, col_types = "c")
#   genome.list <- genome %>% na.omit() %>% data.frame()
#   
#   i=0                                             #Initialize empty objects for faster processing time
#   flag=0
#   y<-data.frame()
#   gene.length<-list()
#   sequence.all<-list()
#   alist<-list()
#   #species.name<-species.list                      #this will be user supplied, must match how it looks in FASTA
#   b<-''
#   c<-''
#   
#   cat(paste("Now importing:",dir, "\n"))        #basic progress checker
#   for(i in genome.list[ ,1]){                 #for each line of inputted genomeFASTA
#     if(str_detect(i,"^>")){                     #if row starts with >
#       if(flag==1){                              #flag check to concatenate sequence lines
#         sequence<-str_c(alist,collapse = "")    #collapse
#         aasize<-nchar(sequence)                 #get amino acid size
#         gene.length<-rbind.data.frame(gene.length,aasize)  #store size
#         sequence.all<-rbind.data.frame(sequence.all,sequence)  #store sequence
#         alist<-list()                           #empty temp variable
#         flag=0                                  #reset flag
#       }
#       flag=0                                    #reset
#       pattern.chr<-paste0("Chromosome:\\S+")#make species name pattern
#       pattern.gn<-paste0("description:.*")             #make gene name pattern
#       pattern.type<-paste0("gene_biotype:\\S+")            #make >tr or >sp pattern 
#       a<-str_extract(i,pattern = pattern.chr) #find and store species name (user supplies species name)
#       b<-str_extract(i,pattern = pattern.gn)    #find and store gene name
#       c<-str_extract(i,pattern = pattern.type)  #find and store pattern type (>tr or >sp, should be only two options)
#       abci<-c(a,b,c,i)                             #make above identifiers into easy to bind object
#       y<-rbind(y,abci)                           #rbind to fill rows of new df
#     }
#     else{
#       flag=1                                    #Raise the flag!
#       b<-str_trim(i,side = c("right"))          #Cut off new line character
#       alist[[i]]<-b                             #Storing all sequence lines of current protein
#     }
#   }
#   if(flag==1){                                  #this catches the last protein's sequence VVV
#     sequence<-str_c(alist,collapse = "")        #
#     aasize<-nchar(sequence)                     #
#     gene.length<-rbind.data.frame(gene.length,aasize) #
#     sequence.all<-rbind.data.frame(sequence.all,sequence) #
#     alist<-list()                               #
#     flag=0                                      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   }
#   y[,5]<-gene.length                             #add sizes to df
#   y[,6]<-sequence.all                             #add sequences to df
#   colnames(y)<-c("Location","Description","Type","Header","Length","Sequence") #add colnames to object 'y'
#   y$Description <- str_replace_all(y$Description, "description:", "")  #remove OS= from species column
#   y$Location <- str_replace_all(y$Location, "Chromosome:", "")  #remove GN= from gene name column
#   y$Type <- str_replace_all(y$Type, "gene_biotype:", "")
#   y$Header <- str_replace_all(y$Header, ">", "")
#   
#   #Count the number of 'genes' that don't have descriptions
#   perc <- (sum(is.na(y$Description)))/nrow(y)*100
#   cat("Percentage of individual sequences w/o description:", perc)
#   
#   genome <<- y
# }   #Load 'chomp' function to import .fasta files
# 
# #Run function on .fasta file to import
# chomp(dir)
# 
# #Take generated genome file and separate into groups based on chromosome strand assignment (1 or -1)
# reverse <- genome %>% filter(str_detect(Location, ":-1$")) %>% mutate(pos = "R", assignment = 2)
# forward <- genome %>% filter(str_detect(Location, ":1$")) %>% mutate(pos = "F", assignment = 1)
# genome.pos <- rbind(reverse, forward)   #This changed the order, which messed up the gene.assignments!
# genome.pos <- genome %>% left_join(genome.pos) %>% mutate(id = str_extract(Header, "^[^\\s]+"))
# 
# #Determine gene distribution
# cat("Genes on forward strand:", nrow(forward), "-->", (nrow(forward)/nrow(genome))*100, "%", 
#     "\nGenes on reverse strand:", nrow(reverse), "-->", (nrow(reverse)/nrow(genome))*100, "%")
# 
#     #Ex. E.coli K-12:
#     #    Genes on forward strand: 2149 --> 49.71085 % 
#     #    Genes on reverse strand: 2174 --> 50.28915 %

#Output new edited .fasta file
#out.fasta <- genome.pos %>% mutate(Header = str_replace_all(Header, " ", "*"))
#write.fasta(sequences = as.list(out.fasta$Sequence),names = out.fasta$Header,file.out = "Ecoli.K12.edit.fasta", open = "w")

#Using package 'seqinr', output two individual .fasta files (one for forward, one for reverse)
#write.fasta(sequences = as.list(reverse$Sequence),names = reverse$Header,file.out = "reverse.Ecoli.K-12.fasta",open = "w")
#write.fasta(sequences = as.list(forward$Sequence),names = forward$Header,file.out = "forward.Ecoli.K-12.fasta",open = "w")

#We now have two individual .fasta files for both the forward and reverse strands of the genome!
#Now to return to AnaCoDa and import both files at two separate 'mixtures'

# genomes <- initializeGenomeObject(file = "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa")
# #genomes <- initializeGenomeObject(file = "Ecoli.K12.edit.fasta")
# 
# #Initialize Parameter Object (Beware of typo in original code)
# parameter <- initializeParameterObject(genome = genomes, sphi = c(0.5, 2), num.mixtures = 2, 
#                                        gene.assignment = genome.pos$assignment, model = "ROC")
# 
# #set mixture.definition = "selectionShared"
# 
# #parameter$fixDEta() #fixing our selection parameter (i.e. assuming neutral selection for strandss)
# 
# model <- initializeModelObject(parameter = parameter, model = "ROC")
# 
# mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)

### Already ran - load in saves files
#runMCMC(mcmc = mcmc, genome = genomes, model = model)
 
# Acceptance rates for Codon Specific Parameters
# Target range: 0.175-0.375 
# Adjustment range: < 0.225 or > 0.325 
# AA	Acc.Rat
# A:	0.162
# C:	0.248
# D:	0.26
# E:	0.204
# F:	0.264
# G:	0.086
# H:	0.288
# I:	0.144
# K:	0.264
# L:	0.072
# N:	0.278
# P:	0.086
# Q:	0.252
# R:	0.054
# S:	0.114
# T:	0.16
# V:	0.134
# Y:	0.262
# Z:	0.264
# Acceptance rate for synthesis rate:
#   Target range: 0.225-0.325 
# Adjustment range: < 0.225 or > 0.325 
# acceptance rates below lower target of 0.225: 383
# acceptance rate above upper target of 0.325: 277

### Running MCMC again, this time with two mixtures
# Acceptance rates for Codon Specific Parameters
# Target range: 0.175-0.375 
# Adjustment range: < 0.225 or > 0.325 
# AA	Acc.Rat
# A:	0.116
# C:	0.242
# D:	0.234
# E:	0.27
# F:	0.282
# G:	0.072
# H:	0.202
# I:	0.124
# K:	0.242
# L:	0.108
# N:	0.24
# P:	0.094
# Q:	0.27
# R:	0.108
# S:	0.15
# T:	0.13
# V:	0.12
# Y:	0.252
# Z:	0.262
# Acceptance rate for synthesis rate:
#   Target range: 0.225-0.325 
# Adjustment range: < 0.225 or > 0.325 
# acceptance rates below lower target of 0.225: 376
# acceptance rate above upper target of 0.325: 360
# ##################################################
# Geweke Score after 25000 iterations: 0.0858739
# ##################################################
# Stopping run based on convergence after 25000 iterations

### Final Status after MCMC end
# Status at thinned sample (iteration): 5000 (50000)
# current logPosterior: -1.31531e+06 
# current logLikelihood: -2.66602e+06
# current stdDevSynthesisRate estimate for selection category 0: 0.443495
# current stdDevSynthesisRate estimate for selection category 1: 0.410679
# current stdDevSynthesisRate proposal width: 0.0347892
# current Mixture element probability for element 0: 0.438402
# current Mixture element probability for element 1: 0.561598
# Acceptance rates for Codon Specific Parameters
# Target range: 0.175-0.375 
# Adjustment range: < 0.225 or > 0.325 
# AA	Acc.Rat
# A:	0.09
# C:	0.216
# D:	0.232
# E:	0.304
# F:	0.238
# G:	0.122
# H:	0.28
# I:	0.07
# K:	0.28
# L:	0.116
# N:	0.262
# P:	0.118
# Q:	0.24
# R:	0.082
# S:	0.056
# T:	0.102
# V:	0.11
# Y:	0.246
# Z:	0.278
# Acceptance rate for synthesis rate:
#   Target range: 0.225-0.325 
# Adjustment range: < 0.225 or > 0.325 
# acceptance rates below lower target of 0.225: 366
# acceptance rate above upper target of 0.325: 292
# ##################################################
# Geweke Score after 50000 iterations: -63.0981
# ##################################################

# #Saving data using C++ methods (which won't work with R methods)
# writeParameterObject(parameter = parameter, file = "test_p.Rda")
# writeMCMCObject(mcmc = mcmc, file = "test_mcmc.Rda")
# 
# #Use this code to load it back in
# parameter <- loadParameterObject(files = "test_p.Rda")
# mcmc <- loadMCMCObject(file = "test_mcmc.Rda")

# #Saving data using C++ methods (which won't work with R methods)
#writeParameterObject(parameter = parameter, file = "parameter_out_Ecoli_split_fixDEta_09172024.Rda")
#writeMCMCObject(mcmc = mcmc, file = "mcmc_out_Ecoli_split_fixDEta_09172024.Rda")

# #Use this code to load it back in
# #parameter <- loadParameterObject(files = "parameter_out_Ecoli_split_09112024.Rda")
# mcmc <- loadMCMCObject(file = "mcmc_out_Ecoli_split_fixDEta_09172024.Rda")
# 
# ### This code was used to import parameter file an alternative way
#   #i <- 1
#   files <- "parameter_out_Ecoli_split_fixDEta_09172024.Rda"
#   
#   #load("parameter_out_Ecoli_split_09112024.Rda")
#   #This function is currently not working - FIXED
# loadROCParameterObject()

#   parameter <- loadROCParameterObject(parameter, files)
# ###
# 
# #Now to view results
# csp_mat_1 <- getCSPEstimates(parameter = parameter, mixture = 1)   #, CSP="Mutation", mixture = 1, samples = 1000)
# csp_mat_2 <- getCSPEstimates(parameter = parameter, mixture = 2)
# 
# head(csp_mat_1)
# head(csp_mat_2)
# 
# #posterior estimates for gene specific parameters (all genes)
# phi_mat <- getExpressionEstimates(parameter = parameter, 
#                                   gene.index = 1:length(genomes),
#                                   samples = 3000)
# phi.df <- phi_mat %>% data.frame() %>% `rownames<-`(genome.pos$id)
# 
# # Histogram of posterior distribution
# ggplot(phi_mat, aes(x = Mean)) +
#   geom_histogram(fill = "blue", color = "black") +
#   labs(title = "Histogram of Values", x = "Phi", y = "Frequency")
# 
# # Density plot of posterior distribution
# plot(density(phi_mat), main = "Density Plot of Phi", xlab = "Phi")
# 
# # #posterior estimates for gene specific parameters (sample 100 random)
# # phi_mat.random <- getExpressionEstimates(parameter = parameter, 
# #                                          gene.index = sample(1:length(genome), 100), 
# #                                          samples = 1000)
# # head(phi_mat.random)
# 
# #Calculate selection coefficient s for each codon and each gene
# selection.coefficients <- getSelectionCoefficients(genome = genomes, 
#                                                    parameter = parameter, 
#                                                    samples = 3000)
# head(selection.coefficients)
# 
# #Now to compare our values to the weights from CAI
# cai.weights <- getCAIweights(referenceGenome = genomes)
# head(cai.weights)
# 
# # GCA       GCC       GCG       GCT       TGC       TGT 
# # 0.5964128 0.7599486 1.0000000 0.4506385 1.0000000 0.8002539 
# 
# nc.per.aa <- getNcAA(genome = genomes)
# head(nc.per.aa)
# 
# #Comparing distribution of selection coefficients to CAI values estimated from ref gene set
# selection.coefficients <- getSelectionCoefficients(genome = genomes, 
#                                                    parameter = parameter, 
#                                                    samples = 3000)
# s <- exp(selection.coefficients)
# cai.weights <- getCAIweights(referenceGenome = genomes)
# 
# codon.names <- colnames(s)
# h <- hist(s[, 1], plot = F)
# plot(NULL, NULL, axes = F, xlim = c(0,1), ylim = range(c(0,h$counts)), 
#      xlab = "s", ylab = "Frequency", main = codon.names[1], cex.lab = 1.2)
# lines(x = h$breaks, y = c(0,h$counts), type = "S", lwd=2)
# abline(v = cai.weights[1], lwd=2, lty=2)
# axis(1, lwd = 3, cex.axis = 1.2)
# axis(2, lwd = 3, cex.axis = 1.2)
# 
# #Now for some diagnostics (Trace plots - how the parameter values change each iteration)
# trace <- getTrace(parameter)
# #trace.s <- parameter$getTraceObject()
# plot.diag.1m <- plot(x = trace, what = "Mutation", mixture = 1) #make figure panel larger to accommodate graph
# plot.diag.2m <- plot(x = trace, what = "Mutation", mixture = 2)
# 
# plot.diag.1s <- plot(x = trace, what = "Selection", mixture = 1) #make figure panel larger to accommodate graph
# plot.diag.2s <- plot(x = trace, what = "Selection", mixture = 2)
# #save plots 15x15 cairo
# 
# #visualize the results of the model fit
# #Ex. use the last 500 samples from mixture 1 for posterior estimate.
# plot(x = model, genome = genomes, samples = 3000, mixture = 1)
# plot(x = model, genome = genomes, samples = 3000, mixture = 2)
# #save plots 10x10 cairo
# 
# #This will compare different 'mixtures' i.e. gene sets
# plot(parameter, what = "Selection", samples = 3000)
# plot(parameter, what = "Mutation", samples = 3000)
# #save plot 8x8 cairo

### NOTE ###
#
#  When exporting a figure to inkscape, use cairo.pdf to preserve Greek symbols (phi)
#
###      ###

### Next step: send figures to Mike; fit gene naming scheme VVV
#(default is to take the first string of characters before white space)
#Idea: for genes with missing description, give identifier, then concatenate (sep = ".")

#Objective: Use package 'Bayesplot' for visualization
#From AnaCoDa, we obtain codon-specific parameters (csp_mat) and gene-specific parameters (phi_mat)
#We also obtain selection coefficients for each codon and each gene.

#Current options for what = ""
#Selection, Mutation, MixtureProbability, Sphi, Mphi, ExpectedPhi, AcceptanceCSP

# ### Example: Store trace information as data.frame
# temp <- as.data.frame(trace$getMixtureProbabilitiesTrace())
# 
# phi_trace <- as.data.frame(trace$getSynthesisRateTrace())
# #This code gets a df where each row is an iteration and columns are ??? but number=n*2 (4323 + 4323)
# 
# #Gets trace of synthesis rates
# temp <- trace$getSynthesisRateTrace()
# #gets gene names
# temp <- getNames(genomes)
# 
# #This captures the trace of each gene's parameters! (x/5001)
# df <- trace$getSynthesisRateTrace()
# df_1 <- df[[1]]
# names(df_1) <- getNames(genomes)
# df.1 <- df_1 %>% data.frame() %>% slice(-(1:2000))
# 
# #This captures the actual estimate of each gene's parameters! (x/4323)
# df <- parameter$getSynthesisRate()
# df_1 <- df[[1]]
# names(df_1) <- getNames(genomes)
# df.1 <- df_1 %>% data.frame()

#Now, we need to capture the mutation and selection parameters for each codon
#df <- parameter$

#Trace is important for diagnostics, parameter is the actual estimate

# #Example Bayesplot method
# library("bayesplot")
# library("rstanarm")
# library("ggplot2")
# 
# fit <- stan_glm(mpg ~ ., data = mtcars)
# posterior <- as.matrix(fit)
# 
# plot_title <- ggtitle("Posterior distributions",
#                       "with medians and 80% intervals")
# 
# #This plots traces of the parameter
# mcmc_areas(posterior,
#            pars = c("cyl", "drat", "am", "wt"),
#            prob = 0.8) + plot_title
# 
# # #Now trying it out on our data - not sure what the ideal estimate is to plot
# # mcmc_areas(
# #   phi.df,
# #   pars = c("Mean", "Mean.log10", "Std.Dev", "log10.Std.Dev"),
# #   prob = 0.8
# # ) + plot_title
# 
# #Plotting the trace
# mcmc_areas(
#   phi.df,
#   pars = c("Mean", "Mean.log10", "Std.Dev", "log10.Std.Dev"),
#   prob = 0.8
# ) + plot_title
# 
# #Could try to construct a matrix with each codon's estimates...
# csp.df.1 <- csp_mat_1 %>% data.frame()
# csp.df.2 <- csp_mat_2 %>% data.frame()
# 
# hist(csp.df.1$Mutation.Mean)
# 
# #Take object of traces and feed into bayesplot
# names(mutationTrace) <- c("mut", "sel")
# 
# #Get list of codons
# names.aa <- aminoAcids()
# names.aa <- setdiff(names.aa, c("W", "X", "M"))
# #remove W, X, and M
# AA.df <- names.aa %>% data.frame(aa = names.aa) %>% select(-.)
#   
# test <- mutate(AA.df, Codon = lapply(AA.df$aa, AAToCodon)) #need to remove last codon from each aa (should end up with 40 total)
# Codon.AA <- test %>% unnest(cols = c("Codon"))
  
#40 entries?
#excluded stops, methionine, tryptophan, and we remove another 20 (ref codon)
# we split serine so 21 amino acids   #SERINE IS SPLIT INTO TWO
#59 potential codons, I'm only seeing 40
#reference codon is last codon for an amino acids

#Useful to have a function that extracts traces in our required format
# main goal is to be looking at the traces of the fits and using other tools
# that other people have developed

################################# Attempt two  ###############################################

#Define .fasta file in working directory
dir <- "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"

#Run function on .fasta file to import
chomp(dir)

#Take generated genome file and separate into groups based on chromosome strand assignment (1 or -1)
reverse <- genome %>% filter(str_detect(Location, ":-1$")) %>% mutate(pos = "R", assignment = 2)
forward <- genome %>% filter(str_detect(Location, ":1$")) %>% mutate(pos = "F", assignment = 1)
genome.pos <- rbind(reverse, forward)
genome.pos <- genome %>% left_join(genome.pos) %>% mutate(id = str_extract(Header, "^[^\\s]+"))

#Determine gene distribution
cat("Genes on forward strand:", nrow(forward), "-->", (nrow(forward)/nrow(genome))*100, "%", 
    "\nGenes on reverse strand:", nrow(reverse), "-->", (nrow(reverse)/nrow(genome))*100, "%")

#Initialize Genome Object
genomes <- initializeGenomeObject(
  file = "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"
  )

#Initialize Parameter Object
parameter <- initializeParameterObject(
  genome = genomes, 
  sphi = c(0.5, 2), 
  num.mixtures = 2, 
  gene.assignment = genome.pos$assignment, 
  model = "ROC",
  mixture.definition = "selectionShared") #share DEta for both mixtures

#Initialize Model Object
model <- initializeModelObject(parameter = parameter, model = "ROC")

#Initialize MCMC
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)

#Set restart settings
setRestartSettings(mcmc = mcmc, filename = "restart_out_Ecoli_split_sharedSelection_09182024", samples = 1000)

#Run MCMC
#runMCMC(mcmc = mcmc, genome = genomes, model = model)

#Saving data using C++ methods (which won't work with R methods)
#writeParameterObject(parameter = parameter, file = "parameter_out_Ecoli_split_sharedSelection_09182024.Rda")
#writeMCMCObject(mcmc = mcmc, file = "mcmc_out_Ecoli_split_sharedSelection_09182024.Rda")

#Load in saved files
mcmc <- loadMCMCObject(file = "mcmc_out_Ecoli_split_sharedSelection_09182024.Rda")

### This code was used to import parameter file an alternative way
#i <- 1
files <- "parameter_out_Ecoli_split_sharedSelection_09182024.Rda"

#load("parameter_out_Ecoli_split_09112024.Rda")
#This function is currently not working - FIXED

parameter <- loadROCParameterObject(parameter, files)
###

#Now for some diagnostics (Trace plots - how the parameter values change each iteration)
trace <- getTrace(parameter)

#plot.diag.1m <- plot(x = trace, what = "Mutation", mixture = 1)
#plot.diag.2m <- plot(x = trace, what = "Mutation", mixture = 2)

#plot.diag.1s <- plot(x = trace, what = "Selection", mixture = 1)
#plot.diag.2s <- plot(x = trace, what = "Selection", mixture = 2)

#visualize the results of the model fit
#plot(x = model, genome = genomes, samples = 3000, mixture = 1)
#plot(x = model, genome = genomes, samples = 3000, mixture = 2)

#This will compare different 'mixtures' i.e. gene sets
#plot(parameter, what = "Selection", samples = 3000)
#plot(parameter, what = "Mutation", samples = 3000)

#Now to grab codon specific parameters
mutationTrace <- trace$getCodonSpecificParameterTrace(0)

#Take object of traces and feed into bayesplot
names(mutationTrace) <- c("mix1", "mix2")

#Get list of codons
names.aa <- aminoAcids()
names.aa <- setdiff(names.aa, c("W", "X", "M")) #remove W, X, and M

#This should be the order conserved in mutationTrace[[i]]
AA.df <- names.aa %>%                 #pipe amino acids
  data.frame(aa = names.aa) %>%       #make a dataframe
  select(-.) %>%                      #remove duplicated column
  mutate(Codon = lapply(aa, AAToCodon)) %>% #make new column with associated codons for each amino acid
  mutate(Codon = map(Codon, ~ .x[-length(.x)])) %>% #remove the reference codon (last alphabetical codon for each AA)
  unnest(cols = c(Codon))                            #unnest column data

#This assumes that the codons are sorted by amino acid, then codon sequence (i.e. A: GCA, GCC, GCG, etc.)
names(mutationTrace[[1]]) <- AA.df$Codon
names(mutationTrace[[2]]) <- AA.df$Codon

#To verify this data matches our plotted trace data, plot a few codons individually to compare
#plot(mutationTrace$mix1$AGA, main = "Mixture 1: AGA", ylab = expression(Delta * M), xlab = "Samples")
#plot(mutationTrace$mix1$AGC, main = "Mixture 1: AGC", ylab = expression(Delta * M), xlab = "Samples")
#plot(mutationTrace$mix1$TTC, main = "Mixture 1: TTC", ylab = expression(Delta * M), xlab = "Samples")

#plot(mutationTrace$mix2$AGA, main = "Mixture 2: AGA", ylab = expression(Delta * M), xlab = "Samples")
#plot(mutationTrace$mix2$AGC, main = "Mixture 2: AGC", ylab = expression(Delta * M), xlab = "Samples")
#plot(mutationTrace$mix2$TTC, main = "Mixture 2: TTC", ylab = expression(Delta * M), xlab = "Samples")

#Making trace data frames and discard burn-in
mutTrace.1 <- mutationTrace$mix1 %>% data.frame() %>% slice(-1:-2000)
mutTrace.2 <- mutationTrace$mix2 %>% data.frame() %>% slice(-1:-2000)

#Moving on to Bayesplot
posterior.1 <- as.matrix(mutTrace.1)                       #fit trace into posterior matrix
posterior.2 <- as.matrix(mutTrace.2)

mcmc_areas(posterior.1,                                  #supply the posterior matrix
           pars = c("GCA", "GCC", "GCG"),                #supply the desired traces to plot
           prob = 0.8) +                                 #supply desired interval
  ggtitle("Posterior distributions: Codons of Alanine [Mixture 1]",  #supply title
          "with medians and 80% intervals") +
  coord_cartesian(xlim = c(-1.5, 0.5))                   #adjust x-axis

mcmc_areas(posterior.2,                                  #supply the posterior matrix
           pars = c("GCA", "GCC", "GCG"),                #supply the desired traces to plot
           prob = 0.8) +                                 #supply desired interval
  ggtitle("Posterior distributions: Codons of Alanine [Mixture 2]",  #supply title
          "with medians and 80% intervals") +
  coord_cartesian(xlim = c(-1.5, 0.5))                   #adjust x-axis

#Directly compare across mixtures
alanine <- data.frame(
  GCA.1 = mutTrace.1$GCA,
  GCA.2 = mutTrace.2$GCA,
  GCC.1 = mutTrace.1$GCC,
  GCC.2 = mutTrace.2$GCC,
  GCG.1 = mutTrace.1$GCG,
  GCG.2 = mutTrace.2$GCG
)

posterior <- as.matrix(alanine)

mcmc_areas(posterior, prob = 0.8) + 
  ggtitle("Alanine Posterior Distributions","with medians and 80% intervals")

#Testing a visualization method
ggplot() +
  geom_point(aes(x = mutTrace.1$GCA, y = mutTrace.2$GCA)) +
  geom_smooth(aes(x = mutTrace.1$GCA, y = mutTrace.2$GCA), method = "lm", color = "blue", se = TRUE) # Correlation line

regression <- lm(mutTrace.2$GCA ~ mutTrace.1$GCA) 
summary(regression)

#Calculate the correlation between each 'pair' of codons (mixture1:mixture2)
correlations <- sapply(names(mutTrace.1), function(col) {
  cor(mutTrace.1[[col]], mutTrace.2[[col]])
})

codon.corr <- data.frame(
  codon = names(correlations),
  corr = correlations,
  r.sq = correlations^2
)

hist(codon.corr$corr)
hist(codon.corr$r.sq)

print(codon.corr)

# Run regression analysis and store results
regression_results <- sapply(names(mutTrace.1), function(col) {
  model <- lm(mutTrace.2[[col]] ~ mutTrace.1[[col]])
  summary(model)$coefficients[2, ]  # Get the coefficient and its statistics
})

# Create a data frame from the results
regression_df <- data.frame(
  Codon = names(mutTrace.1),
  Estimate = regression_results[1, ],
  Std.Error = regression_results[2, ],
  t.value = regression_results[3, ],
  p.value = regression_results[4, ]
)

print(regression_df)

# #Perform One-way ANOVA
# temp <- mutTrace.1 %>% 
#   pivot_longer(cols = 1:ncol(.), names_to = "codon", values_to = "estimate")
# 
# anova_result <- aov(estimate ~ codon, data = temp)  
# 
# summary(anova_result)
# 
# tukey_result <- TukeyHSD(anova_result)
# 
# tukey_df <- as.data.frame(tukey_result$codon)
# 
# tukey_df_sig <- filter(tukey_df, `p adj` < 0.01 )

#Now, import ecoli expression data [Bhatia et al. 2022]
ecoli.data <- load_file("ecoli_expression_data.csv")
ecoli.data[ecoli.data == 0] <- NA
#ecoli.data[ecoli.data < 1] <- NA

e <- ecoli.data %>% select(-Locus_tag) %>% distinct() %>%
  mutate_at(.vars = 2:ncol(.), log) %>% column_to_rownames(var = "Gene_name")

obs.mean.phi <- mean(colMeans(e, na.rm = TRUE))

#Now, back to our trace
phi_trace <- trace$getSynthesisRateTrace()[[1]]
names(phi_trace) <- getNames(genomes)
phi.trace <- phi_trace %>% data.frame() %>% slice(-1:-2000)
est.mean.phi <- mean(colMeans(phi.trace, na.rm = TRUE))

########## Rerunning the simulation with new sphi parameters ##########

#Define .fasta file in working directory
dir <- "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"

#Run function on .fasta file to import
chomp(dir)

#Take generated genome file and separate into groups based on chromosome strand assignment (1 or -1)
reverse <- genome %>% filter(str_detect(Location, ":-1$")) %>% mutate(pos = "R", assignment = 2)
forward <- genome %>% filter(str_detect(Location, ":1$")) %>% mutate(pos = "F", assignment = 1)
genome.pos <- rbind(reverse, forward)
genome.pos <- genome %>% left_join(genome.pos) %>% mutate(id = str_extract(Header, "^[^\\s]+"))

#Determine gene distribution
cat("Genes on forward strand:", nrow(forward), "-->", (nrow(forward)/nrow(genome))*100, "%", 
    "\nGenes on reverse strand:", nrow(reverse), "-->", (nrow(reverse)/nrow(genome))*100, "%")

#Initialize Genome Object
genomes <- initializeGenomeObject(
  file = "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"
)

#Initialize Parameter Object
parameter <- initializeParameterObject(
  genome = genomes, 
  sphi = c(1,1), 
  num.mixtures = 2, 
  gene.assignment = genome.pos$assignment, 
  model = "ROC",
  mixture.definition = "selectionShared") #share DEta for both mixtures

#parameter$fixSphi()

#Initialize Model Object
model <- initializeModelObject(parameter = parameter, model = "ROC")

#Initialize MCMC
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)

#Set restart settings
setRestartSettings(mcmc = mcmc, filename = "restart_out_Ecoli_split_sharedSelection_sphi1_09302024", samples = 2500)

#Run MCMC
#runMCMC(mcmc = mcmc, genome = genomes, model = model)

#Saving data using C++ methods (which won't work with R methods)
#writeParameterObject(parameter = parameter, file = "parameter_out_Ecoli_split_sharedSelection_sphi1_09302024.Rda")
#writeMCMCObject(mcmc = mcmc, file = "mcmc_out_Ecoli_split_sharedSelection_sphi1_09302024.Rda")

########## Now let's rerun it starting where the previous run ended ############
#Define .fasta file in working directory
dir <- "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"

#Run function on .fasta file to import
chomp(dir)

#Take generated genome file and separate into groups based on chromosome strand assignment (1 or -1)
reverse <- genome %>% filter(str_detect(Location, ":-1$")) %>% mutate(pos = "R", assignment = 2)
forward <- genome %>% filter(str_detect(Location, ":1$")) %>% mutate(pos = "F", assignment = 1)
genome.pos <- rbind(reverse, forward)
genome.pos <- genome %>% left_join(genome.pos) %>% mutate(id = str_extract(Header, "^[^\\s]+"))

#Determine gene distribution
cat("Genes on forward strand:", nrow(forward), "-->", (nrow(forward)/nrow(genome))*100, "%", 
    "\nGenes on reverse strand:", nrow(reverse), "-->", (nrow(reverse)/nrow(genome))*100, "%")

#Initialize Genome Object
genomes <- initializeGenomeObject(
  file = "Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.cds.all.fa"
)

#Initialize Parameter Object
parameter <- initializeParameterObject(
  genome = genomes, 
  sphi = c(obs.mean.phi,obs.mean.phi), 
  num.mixtures = 2, 
  gene.assignment = genome.pos$assignment, 
  model = "ROC",
  mixture.definition = "selectionShared")#,
  #init.with.restart.file = "restart_out_Ecoli_split_sharedSelection_sphi1_09302024.rst_5000")

#Initialize Model Object
model <- initializeModelObject(parameter = parameter, model = "ROC")

#Initialize MCMC
mcmc <- initializeMCMCObject(samples = 5000, thinning = 10, adaptive.width=50)

#Run MCMC
#runMCMC(mcmc = mcmc, genome = genomes, model = model)

#Saving data using C++ methods (which won't work with R methods)
#writeParameterObject(parameter = parameter, file = "parameter_out_Ecoli_split_sharedSelection_sphi1_RESTART_5000_09302024.Rda")
#writeMCMCObject(mcmc = mcmc, file = "mcmc_out_Ecoli_split_sharedSelection_sphi1+RESTART_5000_09302024.Rda")

#Saving data
#writeParameterObject(parameter = parameter, file = "parameter_out_Ecoli_split_sharedSelection_sphiObsMeanPhi_10012024.Rda")
#writeMCMCObject(mcmc = mcmc, file = "mcmc_out_Ecoli_split_sharedSelection_sphiObsMeanPhi_10012024.Rda")

#This is the range of mean phi values when sphi is given obs.mean.phi
range(colMeans(phi.trace, na.rm = TRUE))
#[1] 5.132268e-06 1.186245e+15

#Plot selection on y-axis and mutation on x-axis for each codon
selectionTrace <- trace$getCodonSpecificParameterTrace(1)
sel_Trace <- selectionTrace[[1]]
names(sel_Trace) <- AA.df$Codon

#Making trace data frames and discard burn-in
selTrace <- sel_Trace %>% data.frame() %>% slice(-1:-2000)

plot(x = mutTrace.1$GCA, y = selTrace$GCA)
plot(x = mutTrace.2$GCA, y = selTrace$GCA)

#plot
#summation phi y (log)
#sample x

# at sample(i), what is the sum of all phi value
View(rowMeans(phi.trace))

test <- data.frame(rowMeans(phi.trace))

temp <- test %>% data.frame() %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(log = log(rowMeans.phi.trace.)) %>% select(-rowMeans.phi.trace.) %>% 
  plot(title = "test")

post.prob <- mcmc$getLogPosteriorTrace()
log.likelihood <- mcmc$getLogLikelihoodTrace()

#median(each gene's phi) over sample (5,000)

#Now, back to our trace
phi_trace <- trace$getSynthesisRateTrace()[[1]]
names(phi_trace) <- getNames(genomes)
phi.trace <- phi_trace %>% data.frame() %>% slice(-1:-2000)
est.mean.phi <- mean(colMeans(phi.trace, na.rm = TRUE))

#test - plot genes with extrememly large phi values
max(colMeans(phi.trace))

temp <- phi.trace$ENSB.RBWPQul3cZ.YYSP

plot(phi.trace$ENSB.RBWPQul3cZ.YYSP[2000:2100])
plot(log(phi.trace$ENSB.RBWPQul3cZ.YYSP))


