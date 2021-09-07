function [cov,covt] = compute_error_cov(data)

% PURPOSE:
%   This subroutine calculates the error covariance matrix for a lagged
%   ensemble forecast of a user specified climate time series. The user must input
%   the error matrix for the forecast time series. Specifcally, the error matrix
%   should be written such that each column of the matrix  
%   is for a given target day and the rows enteries give error at different leads.
%  
% REFERENCES: 
%
% CALLING SEQUENCE:
%    [cov,covt] = compute_error_cov(data);
%
% INPUTS:
%       data: A matrix of the forecast error for a given target day as a function of lead time and year. 
%       Assumed dimensions are [lead (L), target day (ntime), year (nyr)]. 
%
% OUTPUTS:
%       cov: The [L,L,nyr] lagged error covariance matrix for input data
%       supplied by user. L == lead (or initalization frequency) and
%       nyr == the number of years for which the forecasts exist.
%        
%       covt: The lagged error covariance matrix for input data
%       supplied by user averaged across all forecast years. As before, 
%       L == lead (or initalization frequency) and
%       nyr == the number of years for which the forecasts exist.

[L,ntime,nyr] = size(data);


for ny=1:nyr
  
  % Calculate the lagged error covariance matrix allowing for NaN's in data.  
  % This calculation in matrix notation is equivalent to C =
  % data*data'/ntime.
  
A = squeeze(data(:,:,ny));

    

cov(:,:,ny) =nancov(A'); %Lagged error covariance matrix per year


end

    % Average the lagged error covariance matrix across all years
covt = mean(cov,3);

end
