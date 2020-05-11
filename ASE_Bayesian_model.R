#read in the data
parent.dat <- read.csv("parent35_9679.csv",header=T)
hybrid.dat <- read.csv("hybrid35_9679.csv",header=T)


# Make the data for the model
n.genes <- length(unique(hybrid.dat$Gene)) #number of genes (scalar)
n.obs.off <- dim(hybrid.dat)[1] #number of offspring reads (scalar)
n.obs.parent <- dim(parent.dat)[1] #number of parent reads (scalar)
which.gene <- hybrid.dat$Gene #vector that assigns each offspring observation to a particular gene
is.A <- parent.dat$N.is.SD # vector that indicates if each obs in the parent data is from pop A (0/1)
is.B <- parent.dat$N.is.SC # vector that indicates if each obs in the parent data is from pop A (0/1)
which.gene.parent <- parent.dat$Gene #vector that assigns each offspring observation to a particular gene
N.off.A <- hybrid.dat$N.off.SD #vector of A reads for each offspring observation
tot.reads.off <- hybrid.dat$N.off.total #vector of total reads (A & B) for each offspring observation
N.correct <- as.integer(parent.dat$N.Correct) #number of correct parent reads (eg. A read as A) for each observation
tot.reads.parent <- parent.dat$N.total #vector of total reads (A & B) for each parent observation

# Specify model in BUGS language
#testModel = adjust P for each gene and treatment individually
#no correction for treatment or direction
library(R2jags)
sink("control.jags")
cat("
    model {
    
    ######## DATA NEEDED:  
    # n.genes <- number of genes (scalar)
    # n.obs.off <- number of offspring reads (scalar)
    # n.obs.parent <- number of parent reads (scalar)
    # which.gene <- vector that assigns each offspring observation to a particular gene
    # is.A <- vector that indicates if each obs in the parent data is from pop A (0/1)
    # is.B <- vector that indicates if each obs in the parent data is from pop B (0/1)
    # which.gene.parent <- vector that assigns each offspring observation to a particular gene
    # N.off.A <- vector of A reads for each offspring observation
    # tot.reads.off <- vector of total reads (A & B) for each offspring observation
    # N.correct <- number of correct parent reads (eg. A read as A) for each observation
    # tot.reads.off <- vector of total reads (A & B) for each parent observation
    
    ######## PARAMETERS: 
    # GENE EXPRESSION PROPORTIONS (TRUTH)   
    # P          <- prop of RNA from A
    
    # OBSERVATION ERROR
    # P.obs.A    <- proportion of RNA reported from A that is actually from A
    # P.obs.B    <- proportion of RNA reported from B that is actually from B
    
    
    # Priors and constraints (all parameters treated as fixed and separate for each gene)
    for (i in 1:n.genes){
    P[i] ~ dbeta(1,1) # the base % of each gene in offspring attributable to parent A
    
    
    P.obs.A[i] ~ dbeta(1,1)
    P.obs.B[i] ~ dbeta(1,1)
    }
    
    # DERIVED QUANTITIES 
    
    # Get the true expected proprtion of genes coming from A for each observtion i    
    for (i in 1:n.obs.off){
    logit(P.true[i])<-log(P[which.gene[i]] / (1-P[which.gene[i]]))
    }
    
    # Get expected proportion of misreads for each observation of parent data
    for (i in 1:n.obs.parent){
    logit(P.obs.final[i]) <- log(P.obs.A[which.gene.parent[i]]/(1-P.obs.A[which.gene.parent[i]])) * is.A[i] +  # get the proportion of A mis-reads (classed as A but are actually B)
    log(P.obs.B[which.gene.parent[i]]/(1-P.obs.B[which.gene.parent[i]])) * is.B[i]      # get the proportion of B mis-reads (classed as B but are actually A)
    }
    
    
    # LIKELIHOODS
    
    ##### determine the likelihood of parent observations ################################################
    for (i in 1:n.obs.parent){
    N.correct[i] ~ dbin(P.obs.final[i],tot.reads.parent[i])
    }
    ######################################################################################################	
    
    # now determine the likelihood of offspring observations #############################################
    for (i in 1:n.obs.off){
    # the below two lines calculate the % correct A and B reads for each observation (adjusted for temp treatment of each observation)
    # in other words, how many A are actually A, how many B are actually B
    logit(P.correct.A[i]) <- log(P.obs.A[which.gene[i]]/(1-P.obs.A[which.gene[i]]))
    logit(P.correct.B[i]) <- log(P.obs.B[which.gene[i]]/(1-P.obs.B[which.gene[i]]))
    
    # Proportion of A reads in total
    P.final[i]	 = 	(P.true[i] * P.correct.A[i]) +  		      # % of total reads read as A that are actually A, plus
    ((1-P.true[i]) * (1-P.correct.B[i])) 	    # % of total reads read as B that are actually A
    
    N.off.A[i] ~ dbin(P.final[i], tot.reads.off[i]) #likelihood of the num of A reads given P.a and total reads
    
    }
    #####################################################################################################
    
    }
    ",fill = TRUE)
sink()

# Bundle data
jags.data <- list(n.genes=n.genes,n.obs.off=n.obs.off,n.obs.parent=n.obs.parent,which.gene=which.gene,
                  which.gene.parent=which.gene.parent,N.off.A=N.off.A,
                  tot.reads.off=tot.reads.off,N.correct=N.correct,tot.reads.parent=tot.reads.parent, is.A=is.A,is.B=is.B) 

# Parameters monitored
parameters <- c("P", "P.obs.A","P.obs.B")

# MCMC settings
ni <- 5000
nt <- 10
nb <- 1000
nc <- 3

# Call JAGS from R (run 9679 genes 1 treatment for ~)
gene.ran <- jags(jags.data, inits=NULL, parameters, "control.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

attach.jags(gene.ran)

# check model 
print(gene.ran)

# checking MCMC parameters
plot(gene.ran)

#summarizing data
#can extract results directly from gene.ran
#gene.ran$BUGSoutput$mean $sd $median
library(gtools)
# try plotting posterior distribution
library(ggplot2)
P.norm <- as.data.frame(P[,30])
ggolot()

library(rethinking) #get HDPI from posterior distribution
library(tRophicPosition) #get mode from posterior distribution
HS.ASE <- NULL
CT.ASE <- NULL
HPDI.lower <- NULL
HPDI.upper <- NULL
P.mode <- NULL
P.mean <- round(gene.ran$BUGSoutput$mean$P,4)
P.sd <- round(gene.ran$BUGSoutput$sd$P,4)
P.med <- round(gene.ran$BUGSoutput$median$P,4)

#making raw data file
#table for the P values for future graphing purpose
P.raw <- NULL
for (i in 1:n.genes){
  P.norm <- as.data.frame(P[,i])
  P.raw <- append(P.raw,P.norm)
}
P.raw <- as.data.frame(P.raw)
P.raw <- t(P.raw) #transpose
write.csv(P.raw,"P_rawData_CT.csv")


#making summary file
for (i in 1:n.genes){
  P.norm <- as.data.frame(P[,i])
  P.mode <- append(P.mode,getPosteriorMode(P.norm[,1])) #this step takes a long time
  sig.value <- as.data.frame(HPDI(P.norm,prob=0.95))
  HPDI.lower <- append(HPDI.lower,sig.value[1,1])
  HPDI.upper <- append(HPDI.upper,sig.value[2,1])
  if(sig.value[1,]>0.5) {CT.ASE <- append(CT.ASE,"SD bias")}
  else if (sig.value[2,] < 0.5) {CT.ASE <- append(CT.ASE,"SC bias")}
  else {CT.ASE <- append(CT.ASE,"not significant")}
}

summary <- cbind(1:n.genes,P.mean,P.sd,P.med,P.mode,HPDI.lower,HPDI.upper,CT.ASE,HS.ASE)
write.csv(summary,"ASE20_results_v2.csv")
