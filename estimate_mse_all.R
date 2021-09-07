estimate_mse_all <- function(C){

# PURPOSE:
#   This subroutine calculates the Mean Square Error (MSE) for the user
#   supplied error covariance matrix (C). It is assumed that the MSE are
#   being caluclated for a interval given by dt. For example, if C was
#   calculated using a 4 initalizations per day, then dt =4;
# REFERENCES: 
#
# CALLING SEQUENCE:
#    [ mset ] = estimate_mse_all(C);
#
# INPUTS:
#       C: The lagged error covariance matrix. 
# 
# OUTPUTS:
#       mset: The [nE,L] The Mean Square Error computed for the lagged error covariance matrix. 
#       nE == number of ensemble members and L == lead (or initalization frequency). 
#        

dindex <- dim(C);
nlead <- dindex[1] # Number of leads
mset <- array(,c(dindex[1],dindex[1]))
lindex =1;
for (tau  in 0:(nlead-1)){
    for (nE in 1:(nlead-tau)){
    
    mset[nE,lindex] = (1/(nE)^2)*sum(sum(C[(tau+1):(tau+nE),(tau+1):(tau+nE)]));

    }
    lindex = lindex+1;
}


result<-list(mse = mset)
return(result)

}
