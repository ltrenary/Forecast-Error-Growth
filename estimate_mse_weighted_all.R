estimate_mse_weighted_all <- function(C){

#  PURPOSE:
#   This subroutine calculates the Mean Square Error (MSE) for the user
#   supplied error covariance matrix (C). It is assumed that the MSE are
#   being caluclated for a interval given by dt. For example, if C was
#   calculated using a 4 initalizations per day, then dt =4;
# REFERENCES: 
#
# CALLING SEQUENCE:
#    [ mset ] = estimate_mse_weighted_all(C);
#
# INPUTS:
#       C: The lagged error covariance matrix. 
#
# OUTPUTS:
#       mset: The [nE x L] Mean Square Error computed for the lagged error covariance matrix. 
#       nE == number of ensemble members and L == lead (or initalization frequency). 
#		
#		wts: The [nE,nE,L]  array of the optimal weights.    
# 

dindex <- dim(C)
nlead <- dindex[1] # Number of leads
mset <- array(0,c(dindex[1],dindex[1]))
wt <- array(0,c(dindex[1],dindex[1],dindex[1]))

lindex =1;
for (tau  in 0:(nlead-1)){
    for (nE in 1:(nlead-tau)){

      js <- array(1,c(nE,1)) 
      sig <- C[(tau+1):(tau+nE),(tau+1):(tau+nE)];
      sig
      num <- solve(sig)%*%js
      denom <- t(js)%*%solve(sig)%*%js
      wt1 = num%*%(1/denom);
      
    mset[nE,lindex] = t(wt1)%*%sig%*%wt1;
    wt[1:nE,nE,lindex] = wt1;
    }
    lindex = lindex+1;
}

 result<-list(mse= mset, wts=wt)
 return(result)

}
