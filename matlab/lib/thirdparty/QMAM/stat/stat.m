function [ ] = stat(variable)
%STAT prints some quick statistics about a variable in the command window. 
% 
%% Syntax and description
% 
% stat(variable) prints the following properties of a variable in the
% command window: 
% 
%   size: size of variable.
%   numel: number of elements in variable.
%   NaNs: number of not-a-number values in variable.
%   maximum: maximum value in variable across all dimensions.
%   minimum: minimum value in variable across all dimensions.
%   mean: mean value of variable across all dimensions.
%   median: median value of variable across all dimensions.
%   mode: most frequent value in variable across all dimensions.
%   std dev: standard deviation of variable across all dimensions. 
%   variance: variance of variable across all dimensions. 
%
%% Examples 
% 
% A = magic(5)
% A =
%     17    24     1     8    15
%     23     5     7    14    16
%      4     6    13    20    22
%     10    12    19    21     3
%     11    18    25     2     9
% 
% stat(A)
%  
%  Properties of A:
%  size     = 5  5
%  numel    = 25
%  NaNs     = 0
%  maximum  = 25
%  minimum  = 1
%  mean     = 13
%  median   = 13
%  mode     = 1
%  std dev  = 7.3598
%  variance = 54.1667
% 
% 
% % Also works with NaNs: 
% 
% A(2:4)=NaN
% A =
%     17    24     1     8    15
%    NaN     5     7    14    16
%    NaN     6    13    20    22
%    NaN    12    19    21     3
%     11    18    25     2     9
% 
% stat(A)
%  
%  Properties of A:
%  size     = 5  5
%  numel    = 25
%  NaNs     = 3
%  maximum  = 25
%  minimum  = 1
%  mean     = 13.0909
%  median   = 13.5
%  mode     = 1
%  std dev  = 7.2697
%  variance = 52.8485
% 
%% Author Info
% This function was written by Chad A. Greene (www.chadagreene.com) of the
% University of Texas Institute for Geophysics, October 2014. 
% 
% See also: who, whos, mean, median, mode, size, numel. 

assert(nargin==1,'stat function requires exactly one input variable.')

disp(' ') 
disp([' Properties of ',inputname(1),':'])
disp([' size     = ',num2str(size(variable))])
disp([' numel    = ',num2str(numel(variable))])
disp([' NaNs     = ',num2str(sum(isnan(variable(:))))])

disp([' maximum  = ',num2str(max(variable(isfinite(variable))))])
disp([' minimum  = ',num2str(min(variable(isfinite(variable))))])
disp([' mean     = ',num2str(mean(variable(isfinite(variable))))])
disp([' median   = ',num2str(median(variable(isfinite(variable))))])
disp([' mode     = ',num2str(mode(variable(isfinite(variable))))])
disp([' std dev  = ',num2str(std(variable(isfinite(variable))))])
disp([' variance = ',num2str(var(variable(isfinite(variable))))])
disp(' ') 



end

