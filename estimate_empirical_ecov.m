function [C, coef_odiag, coef_diag] = estimate_empirical_ecov(Ce)
% PURPOSE:
%   This subroutine fits a set of functions to the error covariance matrix for a lagged
%   ensemble forecast. The error covariance matrix is computed using the
%   provided subroutine named 'compute_error_cov'. 
%
%  OFF-DIAGONAL FIT:
%   The off-diagonal elements of the lagged error covariance are fit as a
%   decaying exponetial functions, such that the amplitude (a), decay rate
%   (gamma), and intercept (b) vary in relation to the lead time (initaliztion frequency). 
%   Formally, the model of off-diagonal error growth is written as follows:
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
%
%  ALONG-DIAGONAL FIT:
%   The diagonal elements are modeled by a logistic equation, fit to the residuals between the diagonal elements of the
%   error covariance matrix and the off diagaonal elements. 
%           C_residual(tau) = epsilon_0/(1 + e^-1*alpha(tau-tau_0))
%   Where:
%   epsilon_0 == the maximum error.
%   alpha == the initial error growth rate
%   tau_0 == the inflection point.
%   The residual error growt
%
% REFERENCES: 
%  Lorenz, E. N. (1982), Atmospheric predictaiblity experiments with a large numerical model, Tellus, 34, 505â€“513.
%
% CALLING SEQUENCE:
%    [C,coef_odiag,coef_diag] = compute_error_cov(Ce,dt);
%
% INPUTS:
%       Ce: The lagged error covariance matrix. The matrix is [LxL] dimensions, where L == lead time or initalization frequency. 
%      
% OUTPUTS:
%       C: The [LxL] Empirically derived lagged error covariance matrix. 
%        
%       coef_odiag: The parameters a, gamma, and b from the nonlinear fit of the off-diagonal elements of the lagged error covariance matrix. 
%       coef_diag:  The parameters epsilon, alpha, and tau0 recovered when the diagonal of the lagged error covariance matrix is fit as a logistics equation. 
%
%% FIT MODEL TO THE OFF-DIAGONAL ELEMENTS
[i,j] = size(Ce);
C = NaN(i,j);
for ii =1:i
    t1 = (1:i-ii);
   
    Cs = Ce(ii,ii+1:end);
    
    clen = length(Cs);
    t2 = nan(clen,1);
    t2(:) = ii;
    t2= t2';
    if (ii == 1)
        tt1 =t1;
        tt2 =t2;
        Ct =Cs;
    end
    if (ii >1)
        tt1 = [tt1 t1];
        tt2 = [tt2 t2];
        Ct = [Ct Cs];
    end
    
end
tt1 =tt1';
tt2 =tt2';
tt = [tt1 tt2]; %lead time that are used in the fit.
Ct = Ct'; 


lb = [];
ub = [];
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

B0 = [.01 .04 -.03 .03 .2]; %Specifing inital guess values for the parameters a, gamma, and b. Nonlinear fitting can be sensitive to the inital parameters guess, this may need to be changed. 
F = @(b,tdata)(b(1).*tdata(:,2)+b(2)).*exp((b(3).*tdata(:,2)).*tdata(:,1))+(b(4).*tdata(:,2)+b(5));
 Cfit = lsqcurvefit(F,B0,tt,Ct,lb,ub,options);
      
  a = Cfit(1)*(1:(i-1))+Cfit(2);
  coef_odiag.a(1,1) = Cfit(1);
  coef_odiag.a(2,1) = Cfit(2);
  
  gamma = Cfit(3)*(1:(i-1));
  coef_odiag.gamma(1,1) = Cfit(3);

  b = Cfit(4)*(1:(i-1))+Cfit(5);
  coef_odiag.b(1,1) = Cfit(4);
  coef_odiag.b(2,1) = Cfit(5);

 
 % EVALUATE THE OFF-DIAGONAL FIT 
 for ii=1:(i-1)
      t = (1:i-ii);
     C(ii,ii+1:end) = a(ii)*exp(gamma(ii)*t)+b(ii);
     C(ii+1:end,ii) = a(ii)*exp(gamma(ii)*t)+b(ii);
 end
 
%% FIT MODEL TO THE ALONG-DIAGONAL ELEMENTS 
% Find the extrapolated values from the fit to the off diagonal terms
      d = diag(Ce);
      d = d(1:end-1);
      for ii =1:(i-1)
      dex(ii,1) = a(ii)*exp(gamma(ii)*0.)+b(ii);
      end
     
% Find the residual between the extropolated fit and the actual diagonal. Fit a logistic curve to this residual.      
      res = d-dex;
      lb = [];
      ub = [];
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

    B0 = [.9 .04 20]; %Initial guess values for parameters epsilon, alpha, and tau_0
   
    t = (1:(i-1))';
    
  
    
    F = @(b,tdata)(b(1)./(1+exp(-b(2).*(tdata-b(3)))));
    R = lsqcurvefit(F,B0,t,res,lb,ub,options);

     
    for ii=1:(i-1)
     rfit(ii,1) = R(1)/(1+exp(-R(2)*(ii-R(3)))); %evaluate the fit
    end
      
     
     dfit = rfit+dex; %The fit is to the residual (residual = diagonal-extrapolated)---> diagonal = residual+ extrapolated
     
for ii =1:(i-1)
     C(ii,ii) = dfit(ii);
end

coef_diag.epsilon(1,1) = R(1);
coef_diag.alpha(1,1)= R(2);
coef_diag.tau0(1,1) = R(3);


%PLOTS OF LAGGED ERROR COVARIANCE
figure
clevels = 0:.1:2;

subplot(1,2,1)
contourf(Ce,clevels)
xlabel('i [days]')
ylabel('j [days]')
title({'Cross-lead Error Covariance Matrix'})
caxis([0 2.])
colorbar

subplot(1,2,2)
contourf(C,clevels)
xlabel('i [days]')
ylabel('j [days]')
title({'Emprical Cross-lead Error Covariance Matrix'})
caxis([0 2.])
colorbar


end


