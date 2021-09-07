# Forecast-Error-Growth
The enclosed code provide a general method for determining the optimal lagged ensemble in subseasonal forecasting. 
The mean square error of a lagged ensemble depends only on a quantity called the cross-lead error covariance matrix,
which can be estimated from a short hindcast dataset and parameterized in terms of analytic functions of time. 
The resulting parameterization allows the skill of forecasts to be evaluated for an arbitrary ensemble size and initialization frequency, 
without the need of additional hindcast experiments. On this page, we provide the relevant numerical codes related to this methodology and 
then apply it to forecasts of the Madden Julian Oscillation (MJO) from version 2 of the Climate Forecast System (CFSv2). 
The details of this work can be found here:
Trenary, L., DelSole, T., Tippett, M. K., and Pegion, K. (2017), A new method for determining the optimal lagged ensemble, 
J. Adv. Model. Earth Syst., 9, 291– 306, doi:10.1002/2016MS000838.

