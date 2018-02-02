###########################################################################################################
# Program:      C://Users//MichaelDDiDomenico//Desktop//HUD//FAFSA//Phase III//Random Assignment//Programs//mc.simulation.v4.R
# Programmer:   Mike DiDomenico
# Purpose:      estimate power of our study given our assignment procedure and defined clusters
# Create Date:  4/20/17
# Input:        C://Users//MichaelDDiDomenico//Desktop//HUD//FAFSA//Phase III//Random Assignment//Data//clusters.csv
# Output:

# make sure to update pha.params as necessary to reflect current assumptions in the experimental design

##########################################################################################################
rm(list = ls())

datadir<-"C://Users//MichaelDDiDomenico//Desktop//HUD//FAFSA//Phase III//Random Assignment//Data//"
progdir<-"C://Users//MichaelDDiDomenico//Desktop//HUD//FAFSA//Phase III//Random Assignment//Programs//"

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
# note on adding variation: we are adding variation (u) because we think that there are important 
# differences at the unit of randomization (OES_cluster, think one or a few buildings in close geographic
# proximity) within each PHA. This is plausible because it is common for different high schools to have
# different FAFSA completion rates, and it is likely that OES_clusters overlap to some extent with 
# different high schools. We could use only our estimates for the PHA-specific FAFSA completion rates if
# we thought there was no intra-PHA variation in FAFSA completion rates, but this is unlikely.
# This will result in a conservative estimate of power, but one we are comfortable with both for the 
# purposes of communicating with PHAs and for helping HUD think about how to frame future grant announcements 
# to allow for strong evaluation.
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################



# call in assignment function 
library(hudfafsaiii)
# call in other relevant libraries calculating for robust standard errors
library(sandwich)
library(lmtest)

# parallel libraries
library(parallel)
library(foreach)
library(doSNOW)
library(doParallel)
library(doRNG)

# Specify clusters (how many the machine has)
#cl <- makeCluster(cores)
#registerDoSNOW(cl)
cls <- makeCluster(detectCores())

# set seed for reproducibility
registerDoRNG(seed=123)

# set.seed(3472) # 3472 is random integer chosen using random.org:
# Timestamp: 2017-03-31 17:16:42 UTC


##########################################################################################################
# parameters for intervention
##########################################################################################################

# Overall individual takeup rate ("rate") for students offered services:
rate <- 0.5 # Set at 0.5 as a baseline

# List of PHA states, their number of navigators (nav), and max workload (ml) for each navigator
state <- c("CA", "IL", "PA", "WA", "WI")
nav <- c(3, 3, 2, 3, 1)
ml  <- c(125,150,150,82,125)
pha.params <- cbind.data.frame(state,nav,ml)

# Load data. This is at the amp level in case we want to re-organize the assignment clusters.
clusters.with.amps <- read.csv(paste(datadir,"clusters.with.amps.csv",sep=""))

# Aggregate total number of 15-20 year olds by OES cluster
clusters<-aggregate(tot_age_eligible~OES_cluster+state,data=clusters.with.amps,FUN=sum)

dta <- merge(clusters,pha.params,by="state",all.x=T)

# Create new variables needed for treatment assignments
# Calculate total number of age-eligible to-be-served in each AMP
dta$tot.served <- dta$tot_age_eligible*rate
# Total workload of navigator(s):
dta$tot.load <- 0
# Create a treatment indicator
dta$treat <- 0

# Remove all clusterss w/ fewer than 10 age-eligible youth
dta <- dta[dta$tot_age_eligible >= 10, ]
# Remove all clusters w/ too many age-eligible youth to be serviced by a single navigator
# Not applying to HACLA after consulting with the PHA and HUD
dta <- subset(dta, state == "CA" | dta$tot_age_eligible <= as.numeric(dta$ml)/rate)

# Number of AMPs in population
(N <- dim(dta)[1])


##########################################################################################################
# parameters for data generating process which will be:
# y=b_0+b_1*treat+b.chicago*chicago+b.la*la+b.milwaukee*milwaukee+b.philly*philly+u
##########################################################################################################
# from NCAN, the completion rates by city are the following for the HS class of 2015 (we are still waiting
# to see if we can get PHA specific figures from HUD) 
# https://www.insidehighered.com/sites/default/server_files/files/NCAN%20city%20report.pdf:
# Chicago:      .61
# Los Angeles:  .60
# Milwaukee:    .39
# Philadelphia: .54
# Seattle:      .49


# create PHA specific intercepts (with Seattle as reference)
b.chicago <- .61-.49
b.la <- .60-.49
b.milwaukee <- .39-.49
b.philly <- .54-.49
# mean of baseline fafsa completion (for Seattle)
b0 <- .49

# add state dummies
dta$chicago<-(dta$state=="IL")
dta$la<-(dta$state=="CA")
dta$milwaukee<-(dta$state=="WI")
dta$philly<-(dta$state=="PA")
dta$seattle<-(dta$state=="WA")

# unclear to me why it is the case, but it seems like this function needs to be specified within the program
# and not say in a library or it will not work with the foreach loop


##########################################################################################################
# The k loop will create a power curve for the sequence of effect sizes defined in the effect.size 
# vector. The current range (.05-.10) was chosen to find the lower end of what we think we can detect
# with adequate power. We think it is plausible to see an effect of .10 but want to see if we are able to 
# detect smaller, but still policy relevant changes
##########################################################################################################

# Range of plausible treatment effects for creating a power curve
b1 <- .09

dta$u[dta$state=="IL"] <- 
   runif(n=sum(dta$state=="IL"),-abs(min(b.chicago-0,1-b1-b.chicago)),abs(min(b.chicago-0,1-b1-b.chicago)))
dta$u[dta$state=="CA"] <- 
   runif(n=sum(dta$state=="CA"),-abs(min(b.la-0,1-b1-b.la)),abs(min(b.la-0,1-b1-b.la)))
dta$u[dta$state=="WI"] <- 
   runif(n=sum(dta$state=="WI"),-abs(min(b.milwaukee-0,1-b1-b.milwaukee)),abs(min(b.milwaukee-0,1-b1-b.milwaukee)))
dta$u[dta$state=="PA"] <- 
   runif(n=sum(dta$state=="PA"),-abs(min(b.philly-0,1-b1-b.philly)),abs(min(b.philly-0,1-b1-b.philly)))
dta$u[dta$state=="WA"] <- 
   runif(n=sum(dta$state=="WA"),-abs(min(b0-0,1-b1-b0)),abs(min(b0-0,1-b1-b0)))



# right now this will generate all of the info necessary to calculate the power for a specific 
# effect size. we could add another layer to trace out the power curve.

# use clusterExport to send data to workers, but requires to assign data to global memory (<<-)
# use clusterApply to direct workers to perform tasks
# use Reduce to combine results 

    # need to send the data to the workers to store in memory on each processor cache, so storing as global
    
    # this will pass the split lists of iterations to each of the workers
    gen.ri.estimates<-function(ri.iterations){
        # split up list of iterations, so that we can send chunks to different processors
        rlist<-clusterSplit(cls,1:ri.iterations)
        # this will send the substates data frame to each of the worker processors
        clusterExport(cls,"substates.global")
        # send to the make.assignments function (can build out as needed)
        substates2<-clusterApply(cls,rlist,make.assignments,substates.global,pha.params)
        # return the results in a matrix n=N rows in substates, m=number iterations
        Reduce(rbind,substates2)
    }
    
    # this will loop through the MC simulation and spit out N=length(rlist) null estimates 
    make.assignments <- function(rlist,input, pha.params){
        all.est<-list()
        for (k in rlist){
            # Create a container dset called "substates" (just initialize to input)
            # need to resort the order each iteration
            N<-dim(input)[1]
            input$rand <- runif(n = N)
            input <- input[order(input$state, input$rand), ]
    
            substates2 <- input[1, ]
            # Want the dset to be empty, so empty it
            substates2 <- substates2[-1, ]
    
            for (s in c("IL", "PA", "CA", "WA", "WI")){
                # this preserves the ys from the 'true' run, but will change the treat indicator
                st <- input[input$state == s, ]
                # reset treatment indicator
                st$treat <- 0
                
                # Number of AMPs in PHA s:
                (st_n <- dim(st)[1])
                
                # Number of navigators in PHA s:
                (nn <- pha.params$nav[pha.params$state == s])
                
                # Calculate in loop: (1) running total caseload and (2) max drive time
                
                # We will begin the loop by examining the second AMP in the sort order and calculating
                # the distance (drive time) between AMP2 and AMP1. Then we will proceed through
                # the loop, calculating the maximum distance between AMPi and all previous AMPs.
                
                # Calculate the total load for AMP1 outside the loop:
                st$tot.load[1] <- st$tot.served[1]
                # Assign the first AMP to treatment outside the loop:
                st$treat[1] <- 1
                
                # Initialize index varible k:
                l <- 2
                
                # Begin while loop, calculating as long as the total workload is less than or
                # equal to the maximum load per navigator times the number of navigators:
                while (st$tot.load[l-1] <= pha.params$ml[pha.params$state==s]*nn){
                    # Calculate running total load
                    st$tot.load[l] <- sum(st$tot.served[1:l])
                    # Assign to treatment as long as the while loop is still going
                    st$treat[l] <- 1
                    # Increment l
                    l <- l + 1
                }
                # Build the "substates2" data by stacking finished "st" data frames on top of
                # each other.
                substates2 <- rbind(substates2, st)
            }
            tmp.est <- coefficients(lm(y ~ treat + chicago + la + milwaukee + philly, data=substates2))[2]
    
        all.est<-rbind(all.est,tmp.est)
        }
        all.est
    }


# this will do a one processor outer loop of generating multiple draws of the sample and pass
# the inner loop, the randomization inference loop, to multiple processors and return a matrix of
# estimates for b1, the randomization inference p values, and the HC2 robust p values. It should take
# approximately 1 hr for loops of 1000,1000 using 3 cores (leaving one out for other work, but maybe not 
# necessary)
gen.b.estimates <- function(mc.iterations,ri.iterations){
    ri.p.values<-rep(0,mc.iterations)
    robust.p.values<-rep(0,mc.iterations)
    b1.estimates<-rep(0,mc.iterations)
    for (j in 1:mc.iterations){

        # Generate a random number for each OES_cluster and sort on it 
        N <- dim(dta)[1]
        dta$rand <- runif(n = N)
        dta <- dta[order(dta$state, dta$rand), ]
    
        substates<-assignment (input=dta,pha.params=pha.params)
    
        # adjust the outcome (y) for each cluster based on the data generating parameters defined
        # above
        substates$y <- b0 + b1*substates$treat + b.chicago*substates$chicago + b.la*substates$la +
            + b.milwaukee*substates$milwaukee + b.philly*substates$philly + substates$u
    
        # make sure there are no out of range y values
        stopifnot(max(substates$y)<=1 | min(substates$y>=0))
    
        b.est <- coefficients(lm(y ~ treat + chicago + la + milwaukee + philly, data=substates))[2]
        b1.model <- lm(y ~ treat + chicago + la + milwaukee + philly, data=substates)
        # including hc2 robust p values as extra check
        robust.p.values[j] <- coeftest(b1.model, vcov = vcovHC(b1.model, type = "HC2"))[2,4]
        b1.estimates[j] <- b.est

        # need to send the data to the workers to store in memory on each processor cache, so storing as 
        # global
        substates.global<<-substates

    # this step sends to the multiple processors to complete the randomization inference loops
    system.time(foo<-gen.ri.estimates(ri.iterations))
    # need to make this a list to calculate values without having to use something like reduce
    foo<-unlist(foo)
    # this is a two-sided p-value
    p.value<-sum(abs(median(foo) - foo) > abs(median(foo) - b.est))/length(foo)
    ri.p.values[j]<-p.value
    }
    return(cbind(b1.estimates,ri.p.values,robust.p.values))
}
# system.time(foo<-gen.b.estimates(mc.iterations=50,ri.iterations=10))




# need a wrapper loop here for the power curve
gen.power.curve <- function(effect.list,mc.iterations,ri.iterations){
    effect.size <- effect.list
    mean.effect <- rep(0,length(effect.size))
    ri.power <- rep(0,length(effect.size))
    robust.power <- rep(0,length(effect.size))
    for (i in 1:length(effect.size)){
        # Range of plausible treatment effects for creating a power curve
        b1 <<- effect.size[i]
        dta$u[dta$state=="IL"] <- 
           runif(n=sum(dta$state=="IL"),-abs(min(b.chicago-0,1-b1-b.chicago)),abs(min(b.chicago-0,1-b1-b.chicago)))
        dta$u[dta$state=="CA"] <- 
           runif(n=sum(dta$state=="CA"),-abs(min(b.la-0,1-b1-b.la)),abs(min(b.la-0,1-b1-b.la)))
        dta$u[dta$state=="WI"] <- 
           runif(n=sum(dta$state=="WI"),-abs(min(b.milwaukee-0,1-b1-b.milwaukee)),abs(min(b.milwaukee-0,1-b1-b.milwaukee)))
        dta$u[dta$state=="PA"] <- 
           runif(n=sum(dta$state=="PA"),-abs(min(b.philly-0,1-b1-b.philly)),abs(min(b.philly-0,1-b1-b.philly)))
        dta$u[dta$state=="WA"] <- 
           runif(n=sum(dta$state=="WA"),-abs(min(b0-0,1-b1-b0)),abs(min(b0-0,1-b1-b0)))
        dta <<- dta
        system.time(goo<-gen.b.estimates(mc.iterations,ri.iterations))
        
        mean.effect[i] <- mean(goo[,1])
        ri.power[i] <- 1-sum(goo[,2]>.05)/dim(goo)[1]
        robust.power[i] <- 1-sum(goo[,3]>.05)/dim(goo)[1]
    }   
return(cbind(effect.size,mean.effect,ri.power,robust.power))
}
    
# something is happening with the rows at the top that define b1 and the us. i cannot comment it out without
# throwing an error when trying to run it through. need to figure out what essential beice is missing. Maybe
# solved by making dta global in gen.power.curve
system.time(boo<-gen.power.curve(effect.list=0.05,mc.iterations=1000,ri.iterations=1000))
 
    

# benefits using just 3 cores when going large at the bottom level where the parallel processing happens--
# i increased the number of runs 100 fold, but it increased the time of the run 50 fold
# > system.time(boo<-gen.power.curve(effect.list=seq(0.05, 0.1, 0.01),mc.iterations=10,ri.iterations=10))
#    user  system elapsed 
#    2.31    0.05    5.05 
# > system.time(boo<-gen.power.curve(effect.list=seq(0.05, 0.1, 0.01),mc.iterations=10,ri.iterations=1000))
#    user  system elapsed 
#    2.37    0.03  255.84 


##########################################################################################################
##########################################################################################################
# below are the standard runs to compare timing; could also compare to foreach loops
##########################################################################################################
##########################################################################################################

system.time(foo<-timertest(substates,pha.params))

timertest<-function(input,pha.params){
    all.est<-list()
    for (i in 1:1000){
        # Create a container dset called "substates" (just initialize to input)

        N<-dim(input)[1]
        input$rand <- runif(n = N)
        input <- input[order(input$state, input$rand), ]

        substates2 <- input[1,]
        # Want the dset to be empty, so empty it
        substates2 <- substates2[-1, ]

        for (s in c("IL", "PA", "CA", "WA", "WI")){
            # this preserves the ys from the 'true' run, but will change the treat indicator
            st <- input[input$state == s,]
            # reset treatment indicator
            st$treat <- 0
            
            # Number of AMPs in PHA s:
            (st_n <- dim(st)[1])
            
            # Number of navigators in PHA s:
            (nn <- pha.params$nav[pha.params$state == s])
            
            # Calculate in loop: (1) running total caseload and (2) max drive time
            
            # We will begin the loop by examining the second AMP in the sort order and calculating
            # the distance (drive time) between AMP2 and AMP1. Then we will proceed through
            # the loop, calculating the maximum distance between AMPi and all previous AMPs.
            
            # Calculate the total load for AMP1 outside the loop:
            st$tot.load[1] <- st$tot.served[1]
            # Assign the first AMP to treatment outside the loop:
            st$treat[1] <- 1
            
            # Initialize index varible i:
            i <- 2
            
            # Begin while loop, calculating as long as the total workload is less than or
            # equal to the maximum load per navigator times the number of navigators:
            while (st$tot.load[i-1] <= pha.params$ml[pha.params$state==s]*nn){
                # Calculate running total load
                st$tot.load[i] <- sum(st$tot.served[1:i])
                # Assign to treatment as long as the while loop is still going
                st$treat[i] <- 1
                # Increment i
                i <- i + 1
            }
            # Build the "substates2" data by stacking finished "st" data frames on top of
            # each other.
            substates2 <- rbind(substates2, st)
        }
        tmp.est <- coefficients(lm(y ~ treat + chicago + la + milwaukee + philly, data=substates2))[2]

    all.est<-rbind(all.est,tmp.est)
    }
    return(all.est)
 }



# run time:
# > system.time(boo<-gen.power.curve(effect.list=0.05,mc.iterations=1000,ri.iterations=1000))
#    user  system elapsed 
#   36.27    1.28 4203.06 
# >  
# > 4203/60/60
# [1] 1.1675 hours
