function plot_mse(nmse,dti,dtf)

% PURPOSE:
%   This subroutine plots the Normalized Mean Square Error (NMSE) for a user specified time interval.  For example, if C was
%   calculated using a 4 initalizations per day, then dt =4;
% REFERENCES: 
%
% CALLING SEQUENCE:
%   plot_mse(nmse,dti,dtf);
%
% INPUTS:
%       C: The empirical lagged error covariance matrix extrapolated from the initalization frequency 'dti' to 'dtf'. 
%      dti: The original initalization frequency of the data used in fitting the lagged error covariance matrix. 
%      dtf: The initalization frequency for which the lagged error covariance matrix is extrapolated to. 
%
%
% OUTPUTS:
%      plot of the NMSE, shown as function of lagged ensemble size (horizontal axis) and lead (colored curves - the number denotes the forecast lead in data units).

[ii,jj] =size(nmse);


nli =1;
for nl = 1:dtf:ii
 ntaui =1;
for ntau =dtf:dtf:ii
    nmses(ntaui,nli) = nmse(ntau,nl);
    ntaui = ntaui+1;
end
 nli = nli+1;
end


figure
 clf
colors = jet(30);
 
for nlead = 1:3:30
    aro = num2str(nlead);

plot(1:15,nmses(1:15,nlead),'Color',colors(nlead,:),'linewidth',1.2)

hold on
text(1:15,nmses(1:15,nlead),aro,'color','k')
end
grid

xlabel('Number of days in the lagged ensemble')
ylabel('Normalized MSE')
lstr = ['Initialized  ',num2str(dtf),' times per day'];
title({'NMSE of CSFV2 hindcast of MJO','November-February 1999-2010', lstr})

ylim([0 1.6])

 xlim([1 15])
 plot([1 15], [1 1],'k','linewidth',2)
 







end
