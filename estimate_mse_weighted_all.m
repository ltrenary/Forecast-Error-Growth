function [ mset,wt ] = estimate_mse_weighted_all(C)

% PURPOSE:
%   This subroutine calculates the Mean Square Error (MSE) for the user
%   supplied error covariance matrix (C). It is assumed that the MSE are
%   being caluclated for a interval given by dt. For example, if C was
%   calculated using a 4 initalizations per day, then dt =4;
% REFERENCES: 
%
% CALLING SEQUENCE:
%    [ mset ] = estimate_mse_weighted_all(C);
%
% INPUTS:
%       C: The lagged error covariance matrix. 
%
% OUTPUTS:
%       mset: The [nE,L] The Mean Square Error computed for the lagged error covariance matrix. 
%       nE == number of ensemble members and L == lead (or initalization frequency). 
%        

[i,j] =size(C);

%dbstop 25
lindex =1;
for tau =0:i-1 
    for nE =1:(i-tau)
      js = ones(nE,1);
      sig = C(tau+1:tau+nE,tau+1:tau+nE);
      wt1 = (inv(sig)*js)./(js'*inv(sig)*js);
      
    mset(nE,lindex) = wt1'*sig*wt1;
    wt(1:nE,nE,lindex) = wt1;
    end
    lindex = lindex+1;
end


