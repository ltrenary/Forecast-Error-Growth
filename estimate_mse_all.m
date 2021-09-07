[ mset ] = estimate_mse_all(C)

% PURPOSE:
%   This subroutine calculates the Mean Square Error (MSE) for the user
%   supplied error covariance matrix (C). It is assumed that the MSE are
%   being caluclated for a interval given by dt. For example, if C was
%   calculated using a 4 initalizations per day, then dt =4;
% REFERENCES: 
%
% CALLING SEQUENCE:
%    [ mset ] = estimate_mse_all(C);
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
    
    mset(nE,lindex) = (1/(nE)^2)*sum(sum(C(tau+1:tau+nE,tau+1:tau+nE)));

    end
    lindex = lindex+1;
end






end
