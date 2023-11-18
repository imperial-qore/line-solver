function [mmmpp] = m3pp2m_fit_count(a, bt1, bt2, binf, m3t2, t1, t2, ...
                                    ai, dvt3, t3)
% Fits a second-order Marked MMPP.
% a: arrival rate
% bt1: IDC at scale t1
% bt2: IDC at scale t2
% binf: IDC for t->inf
% m3t2: third central moment
% t1: first time scale
% t2: second time scale
% ai: i-th element is the rate of class i
% dvt3: i-th element is the delta of variance of class i and the variance
%       of all other classes combined, at resolution t3
% t3: third time scale

m = size(ai,1);

[mmpp] = mmpp2_fit_count(a, bt1, bt2, binf, m3t2, t1, t2);

if isempty(mmpp)
    % infeasible
    mmmpp = {};
    return;
end

if size(mmpp{1},1) == 1
   % degenerate case: marked poisson process
   mmmpp = cell(1,2+m);
   mmmpp{1} = mmpp{1};
   mmmpp{2} = mmpp{2};
   for i = 1:m
       mmmpp{2+i} = ai(i);
   end
   return;
end

mmmpp = cell(1,2+m);
mmmpp{1} = mmpp{1};
mmmpp{2} = mmpp{2};

l1 = mmpp{2}(1,1);
l2 = mmpp{2}(2,2);
r1 = mmpp{1}(1,2);
r2 = mmpp{1}(2,1);
q = zeros(2,m);
for i = 1:(m-1)
    a_1 = ai(i);
    dv_1 = dvt3(i);
    t = t3;
    q(1,i) = -(dv_1*r1^4 + dv_1*r2^4 - 2*a_1*r1^4*t - 2*a_1*r2^4*t + 4*dv_1*r1*r2^3 + 4*dv_1*r1^3*r2 + l1*r2^4*t + l2*r1^4*t + 6*dv_1*r1^2*r2^2 + 4*a_1*l1*r2^3*t - 4*a_1*l2*r2^3*t - 8*a_1*r1*r2^3*t - 8*a_1*r1^3*r2*t + 3*l1*r1*r2^3*t + l1*r1^3*r2*t + l2*r1*r2^3*t + 3*l2*r1^3*r2*t - 12*a_1*r1^2*r2^2*t + 3*l1*r1^2*r2^2*t + 2*l1^2*r1*r2^2*t + 2*l1^2*r1^2*r2*t + 3*l2*r1^2*r2^2*t + 2*l2^2*r1*r2^2*t + 2*l2^2*r1^2*r2*t - 8*a_1*l1*r2^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 8*a_1*l2*r2^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l1^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l2^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 8*a_1*l1*r1*r2^2*t + 4*a_1*l1*r1^2*r2*t - 8*a_1*l2*r1*r2^2*t - 4*a_1*l2*r1^2*r2*t - 4*l1*l2*r1*r2^2*t - 4*l1*l2*r1^2*r2*t - 8*a_1*l1*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 8*a_1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 8*l1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/(4*l1*r2*(r1 + r2)*(2*l1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*r1*t - l1*r2*t + l2*r1*t + l2*r2*t));
    q(2,i) = (dv_1*r1^4 + dv_1*r2^4 - 2*a_1*r1^4*t - 2*a_1*r2^4*t + 4*dv_1*r1*r2^3 + 4*dv_1*r1^3*r2 + l1*r2^4*t + l2*r1^4*t + 6*dv_1*r1^2*r2^2 - 4*a_1*l1*r1^3*t + 4*a_1*l2*r1^3*t - 8*a_1*r1*r2^3*t - 8*a_1*r1^3*r2*t + 3*l1*r1*r2^3*t + l1*r1^3*r2*t + l2*r1*r2^3*t + 3*l2*r1^3*r2*t - 12*a_1*r1^2*r2^2*t + 3*l1*r1^2*r2^2*t + 2*l1^2*r1*r2^2*t + 2*l1^2*r1^2*r2*t + 3*l2*r1^2*r2^2*t + 2*l2^2*r1*r2^2*t + 2*l2^2*r1^2*r2*t + 8*a_1*l1*r1^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 8*a_1*l2*r1^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l1^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l2^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*a_1*l1*r1*r2^2*t - 8*a_1*l1*r1^2*r2*t + 4*a_1*l2*r1*r2^2*t + 8*a_1*l2*r1^2*r2*t - 4*l1*l2*r1*r2^2*t - 4*l1*l2*r1^2*r2*t + 8*a_1*l1*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 8*a_1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 8*l1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/(4*(r1 + r2)*(l2^2*r1^2*t - 2*l2^2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1^2*t + l2^2*r1*r2*t + 2*l1*l2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1*r2*t));
end

for i = 1:(m-1)
    mmmpp{2+i} = diag([q(1,i) q(2,i)]) .* mmpp{2};    
end
mmmpp{2+m} = diag([1-sum(q(1,:)) 1-sum(q(2,:))]) .* mmpp{2};