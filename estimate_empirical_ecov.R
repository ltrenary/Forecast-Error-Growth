estimate_empirical_ecov <- function(C,dfi){
	
#.  PURPOSE:
#   This subroutine fits a set of functions to the error covariance matrix for a lagged
#   ensemble forecast. The error covariance matrix is computed using the
#   provided subroutine named 'compute_error_cov'. 
#
#  OFF-DIAGONAL FIT:
#   The off-diagonal elements of the lagged error covariance are fit as a
#   decaying exponetial functions, such that the amplitude (a), decay rate
#   (gamma), and intercept (b) vary in relation to the lead time (initaliztion frequency). 
#   Formally, the model of off-diagonal error growth is written as follows:
#              C(i,j) = a_min[i,j]*exp(-1*(gamma_min[i,j])*|i-j|)+b_min[i,j]
#   Where: 
#   i,j == are index values of the lagged error covariance matrix, contingent 
#.         on i not being equal to j;
#
#   C(i,j) == the lagged error covariance matrix.
#
#   a == Amplitude of the expoential function. This parameter is a function
#   of lead time. 
#
#   gamma == Decay of the exponential function.  This parameter is a
#   function of lead time.
#
#   b == Saturation value of the off-digaonal elements. This parameter is a
#   function of lead time.
#   
#   The parameters a,gamma, and b are linear functions of lead time (tau).
#       a_tau = (beta_a)*tau + alpha_a;
#       gamma_tau = (beta_gamma)*tau ;
#       b_tau = (beta_b)*tau + alpha_b;
#   Where the beta terms are the slope parameters and the alpha terms
#   represent the intercepts. The variable tau represents the lead time.
#   Note that intercept term for gamma is set to zero, to ensure that
#   beta_gamma is always positive. 
#
#
#  ALONG-DIAGONAL FIT:
#   The diagonal elements are modeled by a logistic equation, 
#.  fit to the residuals between the diagonal elements of the
#   error covariance matrix and the off diagaonal elements. 
#           C_residual(tau) = epsilon_0/(1 + e^-1*alpha(tau-tau_0))
#   Where:
#   epsilon_0 == the maximum error.
#   alpha == the initial error growth rate
#   tau_0 == the inflection point.
#   The residual error growt
#
# REFERENCES: 
#  Lorenz, E. N. (1982), Atmospheric predictaiblity experiments with a 
#.                       large numerical model, #Tellus, 34, 505â€“513.
#
# CALLING SEQUENCE:
#    C.FIT <- compute_error_cov(Ce,dt);
#
# INPUTS:
#       Ce: The lagged error covariance matrix. The matrix is [LxL] dimensions, 
#           where L == lead time or initalization frequency. 
#      
# OUTPUTS:
#       C.FIT$Ce: The [LxL] Empirically derived lagged error covariance matrix. 
#        
#       C.FIT$coef.odiag: The parameters a, gamma, and b from the nonlinear fit 
#                   of the off-diagonal elements of the lagged error covariance matrix. 
#       C.FIT$coef.diag:  The parameters epsilon, alpha, and tau0 recovered when the 
#                  diagonal of the lagged error covariance matrix is fit as a logistics equation.
# 
#		C$iday and C$jday: The dimensions of the lagged error covariance matrix.

	
	library(nlmrt)
	mindex <- dim(C)
	
##############################################################	 
#FIT THE NON-LINEAR MODEL TO THE OFF-DIAGONAL ELEMENTS	
##############################################################	

 for (n in 1:(mindex[1]-1)){
    Ctmp <-array(,c((mindex[1]-n),1))
    Ctmp[,1] <- C[n,(n+1):mindex[1]]
    
    clen <- length(Ctmp)
    t2 <-array(,c(clen,1))
    t2[,1] = n
    
    t1 <-array(,c(clen,1))
    t1[,1] = 1:clen
    
    if (n == 1)
    {Ct = Ctmp
    	 tt1 =t1
    	 tt2 =t2}
    
    if ( n > 1)
    {Ct =rbind(Ct, Ctmp)
    	tt1 =rbind(tt1,t1)
    	tt2 =rbind(tt2,t2)}
    
    }
    
    
   #SPECIFY THE MODEL
   regmod <- Ct ~ (p1*tt2+p2)*exp((p3*tt2*tt1))+(p4*tt2+p5)
   #SPECIFY THE INITAL GUESS VALUES FOR THE PARAMETERS. THE FIT CAN BE SENSTIVE TO    THESE CHOICES. 
   strt <- c(p1 =0.1, p2 =0.04, p3= -0.03,p4=0.03,p5=.2)

Cdata <- data.frame(Ct =Ct, tt1 =tt1, tt2=tt2)
Cfit <-nlxb(regmod,start=strt, trace=FALSE, data=Cdata)
 par <- coef(Cfit)
    
    
    #parameters
    a <- array(,c((mindex[1]),1))
    a[,1] <- par[1]*(1:(mindex[1]))+par[2]
   
    coef_odiag <- array(,c(2,3))   
    coef_odiag[1,1] = par[1]
    coef_odiag[2,1] =par[2]
    
    
    gamma <- array(,c((mindex[1]),1))
    gamma[,1] <- par[3]*(1:(mindex[1]))
     
    coef_odiag[1,2] = par[3]
      
    
    b <- array(,c((mindex[1]),1))
    b[,1] <- par[4]*(1:(mindex[1]))+par[5]
    coef_odiag[1,3] = par[4]
    coef_odiag[2,3] =par[5]
    
    C2 <- array(,c(mindex[1],mindex[1]))
    
for (ii in 1:(mindex[1]-1)){
      t = 1:(mindex[1]-ii);
     C2[ii,(ii+1):mindex[1]] = a[ii]*exp(gamma[ii]*t)+b[ii];
     C2[(ii+1):mindex[1],ii] = a[ii]*exp(gamma[ii]*t)+b[ii];
 } 
 
 
##############################################################	 
#FIT THE NON-LINEAR MODEL TO THE ALONG-DIAGONAL ELEMENTS	
##############################################################	
      d <- diag(C);
      d = d[1:(mindex[1]-1)];
      
      dex <- array(,c((mindex[1]),1))
      for (ii in 1:(mindex[1])) {
      dex[ii,1] = a[ii]*exp(gamma[ii]*0.)+b[ii];
      }
    
   
   tt = seq(dfi,(mindex[1]-1),dfi)
   Rt = d - dex[1:(mindex[1]-1)]
# SPECIFY THE MODEL   
regmod <- Rt ~ b1/(1.+exp(-b2*(tt-b3)))

#SPECIFY THE INITAL GUESS VALUES FOR THE PARAMETERS. THE FIT CAN BE SENSTIVE TO    THESE CHOICES.
strt <- c(b1 =0.4, b2 =15.7, b3= 16.)

edata <- data.frame(Rt =Rt, tt =tt)
Cdiag <-nlxb(regmod,start=strt, trace=FALSE, data=edata)
 par <- coef(Cdiag)

     rfit <- array(,c((mindex[1]),1))
    for (ii in dfi:dfi:(mindex[1])) {
    	rfit[ii,1] = par[1]/(1+exp(-par[2]*(ii-par[3]))); #evaluate the fit
    }
     
         
     
     dfit <- rfit+dex; #The fit is to the residual (residual = diagonal-extrapolated)---> diagonal = residual+ extrapolated
     
   for (ii in 1:(mindex[1])) {
   	C2[ii,ii] = dfit[ii];
   }
     
end

coef_diag <- array(,c(1,3)) 
coef_diag.alpha <- array(,c(1,1)) 

coef_diag[1,1] = par[1];
coef_diag[1,2]= par[2];
coef_diag[1,3] = par[3];
 
 iday = seq(dfi,mindex[1],dfi)
 jday = seq(dfi,mindex[1],dfi)

 
 result<-list(Ce= C2, coef.odiag = coef_odiag, coef.diag = coef_diag,iday=iday,jday=jday)
 return(result)
}
