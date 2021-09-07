plot_nmse <- function(nmse,dti,dtf){

# PURPOSE:
#   This subroutine plots the Normalized Mean Square Error (NMSE) 
#.  for a user specified time interval.  For example, if C was
#   calculated using a 4 initalizations per day, then dt =4;
# REFERENCES: 
#
# CALLING SEQUENCE:
#   plot_mse(nmse,dti,dtf);
#
# INPUTS:
#       C: The empirical lagged error covariance matrix extrapolated 
#.          from the initalization frequency 'dti' to 'dtf'. 
#      dti: The original initalization frequency of the data used in 
#           fitting the lagged error covariance matrix. 
#      dtf: The initalization frequency for which the lagged error 
#           covariance matrix is extrapolated to. 
#
#
# OUTPUTS:
#      plot of the NMSE, shown as function of lagged ensemble size (horizontal axis) 
#      and lead (colored curves - the number denotes the forecast lead in data units).


dindex <- dim(nmse)
nmses <- array(0,c(dindex[1],dindex[1]))
nli =1;

# This portion of the code isolates the MSE for the unit time based on the original data
for (nl in seq(from=dtf, to=dindex[1], by=dtf)){
 ntaui =1;
for (ntau  in seq(from=dtf, to=dindex[1], by=dtf)){
    nmses[ntaui,nli] = nmse[ntau,nl];
    ntaui = ntaui+1;
}
 nli = nli+1;
}

dev.new()
ntmp = nmses[1:15,1:28] #Specify a max. lagged ensemble size of 15 and lead times from 1-28
yrange = range(ntmp)
xrange = c(1,15)
n.all  = seq(1,28,3)
col.say = colorRampPalette(c('blue','red'))(length(n.all))
for ( n in n.all) {
	if ( n == 1) {
		plot(1:15,ntmp[,n],type='l',xlim=xrange,ylim=yrange,xlab='ensemble size',ylab='NMSE',col=col.say[which(n==n.all)])
		dd <- array(n,c(15,1))
		text(1:15,ntmp[,n],sprintf("%i", dd))
	} else {
		par(new=TRUE)
		plot(1:15,ntmp[,n],type='l',xlim=xrange,ylim=yrange,xlab='',ylab='',axes=FALSE,col=col.say[which(n==n.all)])
		dd <- array(n,c(15,1))
		text(1:15,ntmp[,n],sprintf("%i", dd))
	}
}
       
}
