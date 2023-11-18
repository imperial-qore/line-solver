function [I,H,logI] = laplaceapprox(h,x0)
% I = laplaceapprox(f,x0)
% approximates I=int f(x)dx by Laplace approximation at x0
% example:  I = laplaceapprox(@(x) prod(x),[0.5,0.5])
d = length(x0); % number of dimensions
tol = 1e-5;
H = num_hess(@(x) log(h(x)), x0, tol);
detnH=det(-H);
if detnH<0
    tol = 1e-4;
    H = num_hess(@(x) log(h(x)), x0, tol);
    detnH=det(-H);
end
if detnH<0
    tol = 1e-3;
    H = num_hess(@(x) log(h(x)), x0, tol);
    detnH=det(-H);
end
if detnH<0
    warning('laplaceapprox.m: det(-H)<0')
end
I = h(x0) * sqrt((2*pi)^d/detnH);
logI = log(h(x0)) + (d/2)*log(2*pi) - log(detnH);
end