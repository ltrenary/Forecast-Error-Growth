compute_error_cov <- function(data,dfi){
	
##############################################################	
# PURPOSE:
# This subroutine calculates the error covariance matrix for a   
# lagged ensemble forecast of a user specified climate time     
# series. The user must input the error matrix for the forecast  
# time series. Specifcally, the error matrix should be written   
# such that each column of the matrix is for a given target day  
# and the rows enteries give error at different leads.
#  
# REFERENCES: 
#
# CALLING SEQUENCE:
#    C <- compute_error_cov(data);
#
# INPUTS:
#       data: A matrix of the forecast error for a given target day as a 
#             function of lead time and year. Assumed dimensions are 
#             [lead (L), target day (ntime), year (nyr)]. 
#       dfi:  Initalization frequency. This is only used for plotting purposes.
#
# OUTPUTS:
#       C$COVy: The [nlead,nlead,nyr] lagged error covariance matrix for input data
#       supplied by user. nlead == lead (or initalization frequency) and
#       nyr == the number of years for which the forecasts exist.
#        
#       C$COV: The lagged error covariance matrix [nlead X nlead] for input data
#       supplied by user averaged across all forecast years. As before, 
#       nlead == lead (or initalization frequency) and
#       nyr == the number of years for which the forecasts exist.
#       
#		C$iday and C$jday: The dimensions of the lagged error covariance matrix.
#################################################################	

dindex <- dim(data) 
nlead <- dindex[1] # Number of leads
nfday <- dindex[2] # Number of target days
nyr <- dindex[3]   # Number of years

COV <- array(, c(nlead, nlead))
COVy <- array(, c(nlead, nlead, nyr))
	
for (ny in 1:nyr)
{
	dum = data[, ,ny]
	COVy[, ,ny] = var(t(dum), na.rm =TRUE)	
}
 # FIND THE AVERAGE LAGGED ERROR COVARIANCE MATRIX ACROSS ALL YEARS
 COV <- apply(COVy, c(1,2), function(x) mean(x))

 iday = seq(dfi,nlead,dfi)
 jday = seq(dfi,nlead,dfi)

 result<-list(C= COV, Cy=COVy,iday=iday,jday=jday)
 return(result)
 
}
