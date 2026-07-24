function [ lambda1,lambda2,a ] = estimateH2( lambda,CX2 )

% This function generates a 2 phase Hyper-exponential distribution
% according to the given mean rate: lambda, and the given squared
% coefficient of variation: CX2

u1 = 1/lambda;
%CX2 is the squares coefficient of variation
a = (1+sqrt((CX2-1)/(CX2+1)))/2;
lambda1 = 2*a/u1;
lambda2 = 2*(1-a)/u1;
end

