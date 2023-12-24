function [APH] = aph_bernstein(f,order)
% APH = = aph_bernstein(f,order)
%
% Fits an APH using Bernstein's approximation (Horvath-Vicario, 2023)
% 
% Input:
% f: function handle for the density function f(x) to approximation
% order: approximation order (e.g., 100)
%
% Output:
% APH: acyclic phase type distribution of order n fitting f
%
% References:
% András Horváth, Enrico Vicario: Construction of Phase Type Distributions 
% by Bernstein Exponentials. EPEW 2023: 201-215.
%
% Example:
% lambda=5;
% k = 2;
% n = 100;
% f = @(x) lambda^k*x^(k-1)*exp(-lambda*x)/factorial(k-1); % erlang-k
% 
% APH = aph_bernstein(f,n)
%
% ExactMean = k/lambda
% ApproximateMean = map_mean(APH)
% xset = 0.01:0.1:4;
% plot(xset,arrayfun(f,xset))
% if n<=100 % map_pdf gets very slow otherwise
%     hold all
%     g=@(x) map_pdf(MAP,x);
%     plot(xset,arrayfun(g,xset))
%     legend('Exact','Bernstein')
% end

n = order;
% Bernstein approximation of order n
c = 0;
for i=1:n
    c = c + f(-log(i/n))/i;
end
T = diag(-[1:n]) + diag([1:(n-1)],1);
alpha = zeros(1,n);
for i=1:n
    alpha(i) = f(-log(i/n))/(i*c);
end
% Conversion to D0,D1 format
P = repmat(alpha,n,1);
APH = {T,-T*P};

% % Verification of results
