function [Cn] = extrapolate_empirical_ecov(coef_odiag,coef_diag,dti,dtf,ntot)
% PURPOSE:
%   This subroutine extrapolates the empirical fit of the error covariance found for the initalization frequency 'dti' to 
%   a finer initalization frequency (dtf). 
%
%   The off-diagonal fit is assumed to be an exponentially decreasing
%   function of the form:
%              C(i,j) = a_min[i,j]*exp(-1*(gamma_min[i,j])*|i-j|)+b_min[i,j]
%   Where: 
%   i,j == are index values of the lagged error covariance matrix, contingent on i not being equal to j;
%
%   C(i,j) == the lagged error covariance matrix.
%
%   a == Amplitude of the expoential function. This parameter is a function
%   of lead time. 
%
%   gamma == Decay of the exponential function.  This parameter is a
%   function of lead time.
%
%   b == Saturation value of the off-digaonal elements. This parameter is a
%   function of lead time.
%   
%   The parameters a,gamma, and b are linear functions of lead time (tau).
%       a_tau = (beta_a)*tau + alpha_a;
%       gamma_tau = (beta_gamma)*tau ;
%       b_tau = (beta_b)*tau + alpha_b;
%   Where the beta terms are the slope parameters and the alpha terms
%   represent the intercepts. The variable tau represents the lead time.
%   Note that intercept term for gamma is set to zero, to ensure that
%   beta_gamma is always positive. 
%
%   The diagonal elements are modeled by a logistic funtion, fit to the residuals between the diagonal elements of the
%   error covariance matrix and the off diagaonal elements. 
%           C_residual(tau) = epsilon_0/(1 + e^-1*alpha(tau-tau_0))
%   Where:
%   epsilon_0 == the maximum error.
%   alpha == the initial error growth rate
%   tau_0 == the inflection point.

% REFERENCES: 
%  Lorenz, E. N. (1982), Atmospheric predictaiblity experiments with a large numerical model, Tellus, 34, 505â€“513.
%
% CALLING SEQUENCE:
%    [Cn] = extrapolate_empirical_ecov(coef_odiag,coef_diag,dti,dtf,ntot);
%
% INPUTS:
%       coef_odiag: The parameters a, gamma, and b from the nonlinear fit of the off-diagonal elements of the lagged error covariance matrix. 
%       coef_diag:  The parameters epsilon, alpha, and tau0 recovered when the diagonal of the lagged error covariance matrix is fit as a logistics equation. 
%       dti:  The initalization frequency for data used in estimating the empirical fit of the error covariance matrix. 
%             For example if using a daily initalization frequency, dt =1. If initalized every five days, dt = 5. 
%       dtf:  The new initalization frequency for which the emprical fit
%             will be extrapolated. 
%       ntot: The total number of leads times.

% OUTPUTS:
%       C: The [ntotxntot] Empirically derived lagged error covariance matrix extroplated to have an initalization frequency of dtf. 
%        

ntt = ntot/dtf;
Cn = NaN(ntt,ntt);
t = dtf/dti:dtf/dti:ntot;

%Specifying the new initalization frequency, relative to the original fit

%% FIND THE OFF-DIAGONAL COEFFICENTS FOR THE NEW INITALIZATION FREQUENCY
at = coef_odiag.a(1,1)*t+coef_odiag.a(2,1);
bt = coef_odiag.b(1,1)*t+coef_odiag.b(2,1);
gt = coef_odiag.gamma(1,1)*t;


for n=1:ntt-1
    tt =t(1:ntt-n);
    ft = at(n)*exp(gt(n)*tt)+bt(n);
    Cn(n+1:end,n) = ft';
    Cn(n,n+1:end) = ft';
end

%% FIND THE ALONG-DIAGONAL COEFFICENTS FOR THE NEW INITALIZATION FREQUENCY

d = coef_diag.epsilon./(1+exp(-coef_diag.alpha*(t-coef_diag.tau0)));
     dt = d+at+bt;
  
 for n=1:ntt   
     Cn(n,n) = dt(n);
 end


 
figure
clevels = 0:.1:2;


contourf(Cn,clevels)
xlabel('i [days]')
ylabel('j [days]')
title({'Extrapolated Cross-lead Error Covariance Matrix'})
caxis([0 2.])
colorbar

end
