###########################################################################################################
# Program:      C://Users//MichaelDDiDomenico//Desktop//HUD//FAFSA//Phase III//Random Assignment//Programs//assignment.function.R
# Programmer:   Mike DiDomenico
# Purpose:      This loop will create the assignments for HUD FAFSA Phase III. The input can be renamed
#               in the call, and the output can be renamed as the object that stores the function output
# Create Date:  4/21/17
##########################################################################################################

assignment.function <- function(input){
    # Create a container dset called "substates" (just initialize to input)
    substates <- input[1, ]
    # Want the dset to be empty, so empty it
    substates <- substates[-1, ]
    # Now we have a dset with all the same variables as input, but empty; we'll fill it
    # up with data once the random sorting and allocation have allowed us to make treatment
    # assignments.
    
    for (s in c("IL", "PA", "CA", "WA", "WI")){
    	# Deal only with one state at a time:
    	st <- input[input$state == s, ]
    	
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
    	# Build the "substates" data by stacking finished "st" data frames on top of
    	# each other.
    	substates <- rbind(substates, st)
    }

   	return(substates)
}
