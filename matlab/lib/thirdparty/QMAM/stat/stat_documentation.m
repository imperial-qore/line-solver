%% |stat| 
% This function prints some quick statistics about a variable in the
% command window.
% 
%% Syntax and Description
% 
% |stat(variable)| prints the following properties of a variable in the
% command window: 
% 
% * |size|: size of variable.
% * |numel|: number of elements in variable.
% * |NaNs|: number of not-a-number values in variable.
% * |maximum|: maximum value in variable across all dimensions.
% * |minimum|: minimum value in variable across all dimensions.
% * |mean|: mean value of variable across all dimensions.
% * |median|: median value of variable across all dimensions.
% * |mode|: most frequent value in variable across all dimensions.
% * |std dev|: standard deviation of variable across all dimensions. 
% * |variance|: variance of variable across all dimensions. 
%
%% Examples 

A = magic(5)

%% 
stat(A)

%%
% This function also works with |NaN| values: 

A(2:4)=NaN

%% 

stat(A)

%% Author Info
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the
% University of Texas Institute for Geophysics (<http://www.ig.utexas.edu/people/students/cgreene/ UTIG>), October 2014. 