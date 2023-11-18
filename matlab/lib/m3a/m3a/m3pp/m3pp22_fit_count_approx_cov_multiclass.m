function [FIT] = m3pp22_fit_count_approx_cov_multiclass(mmpp, ai, st3, t3)
% Fits a M3PP(2,2) given the underlying MMPP(2).
% INPUT
% - mmpp: underlying mmpp(2)
% - ai: rates of the two classes
% - st3: count covariance between the two classes at scale t3
% - t3: third time scale

% number of classes
m = size(ai,1);

if m > 2
   error('No more than two classes supported'); 
end

FIT = mmpp;

% degenerate case: poisson process
if size(FIT{1},1) == 1
    D0 = FIT{1};
    D1 = FIT{2};
    FIT = cell(1,2+m);
    FIT{1} = D0;
    FIT{2} = D1;
    a = sum(ai);
    pi = ai/a;
    for i = 1:m
        FIT{2+i} = pi(i)*D1;
    end
    return;
end

% just a single class!
if m == 1
    FIT = {FIT{1},FIT{2},FIT{2}};
    return;
end

l1 = FIT{2}(1,1);
l2 = FIT{2}(2,2);
r1 = FIT{1}(1,2);
r2 = FIT{1}(2,1);
t = t3;
a1 = ai(1);
w0 =  (2*r1*(1 - exp(-(r1+r2)*t)- (r1+r2)*t)*(a1^2*(r1+r2)-a1*r2*(l1-l2)))/(r2*(r1+r2)^3);
w1 = -(2*r1*(1 - exp(-(r1+r2)*t)- (r1+r2)*t)*(2*a1*l2*(r1+r2)-l2*r2*(l1-l2)))/(r2*(r1 + r2)^3);
w2 = ((2*l2^2*r2*t)*(r1+r2) + 2*l2^2*r1*(1 - exp(-(r1+r2)*t)))/(r2*(r1 + r2)^2) - (2*l2^2*t)/r2;
w3 = (r1+r2)/(l2*r1);
w4 = (l1*r2)/(l2*r1);

% bounds for the first and the second root
L1 = -inf; L2 = -inf;
U1 = +inf; U2 = +inf;

% if set to 1, the first (second) root is never feasible
infeasible1 = 0;
infeasible2 = 0;

% defined for convenience
z = w0 - w1^2/(4*w2);

% impose square root argument is >= 0
if w2 > 0
    L1 = max(L1, z);
    L2 = max(L2, z);
elseif w2 < 0
    U1 = min(U1, z);
    U2 = min(U2, z);
end

% impose q2 >= 0
% first root
if w1 >= 0
    L1 = max(L1, w0);
elseif w2 < 0
    infeasible1 = 1;
end
% second root
if w1 <= 0
    U2 = min(U2, w0);
elseif w2 > 0
    infeasible2 = 1;
end

% impose q1 >= 0
tmp = 2*a1*w3*w2 + w1;
% first root
if tmp >= 0
    U1 = min(U1, z + tmp^2/(4*w2));
elseif w2 > 0
    infeasible1 = 1;
end
% second root
if tmp <= 0
    L2 = max(L2, z + tmp^2/(4*w2));
elseif w2 < 0
    infeasible2 = 1;
end

% impose q2 <= 1
tmp = 2*w2+w1;
% first root
if tmp >= 0
    U1 = min(U1, z + tmp^2/(4*w2));
elseif w2 > 0
    infeasible1 = 1;
end
% second root
if tmp <= 0
    L2 = max(L2, z + tmp^2/(4*w2));
elseif w2 < 0
    infeasible2 = 1;
end

% impose q1 <= 1
tmp = (2*a1*w2*w3-2*w2*w4+w1);
% first root
if tmp >= 0
    L1 = max(L1, z + tmp^2/(4*w2));
elseif w2 < 0
    infeasible1 = 1;
end
% second root
if tmp <= 0
    U2 = min(U2, z + tmp^2/(4*w2));
elseif  w2 > 0
    infeasible2 = 1;
end

if infeasible1 && infeasible2
    error('Empty feasibility region. This should not happen.');
end

% compute feasible covariance
if infeasible2
    sigma = max(min(st3,U1), L1);
    root = 1;
elseif infeasible1
    sigma = max(min(st3,U2), L2);
    root = 2;
else
    sigma1 = max(min(st3,U1), L1);
    sigma2 = max(min(st3,U2), L2);
    if abs(sigma1-st3) < abs(sigma2-st3)
        sigma = sigma1;
        root = 1;
    else
        sigma = sigma2;
        root = 2;
    end
end

% compute parameters
if root == 1
    q2 = (-w1 + sqrt(w1^2 - 4*w2*(w0-sigma)) )/(2*w2);
else
    q2 = (-w1 - sqrt(w1^2 - 4*w2*(w0-sigma)) )/(2*w2);    
end
q1 = ( a1*(r1 + r2) - l2*q2*r1 )/(l1*r2);

% check feasibility just to be sure
tol = 1e-8;
if (q1 >= tol && q1 <= 1+tol) && (q2 >= tol && q2 <= 1+tol)
    q1 = min(max(q1,0),1);
    q2 = min(max(q2,0),1);
else
    error('Parameters are infeasible. This should not happen.');
end

% assemble M3PP[2]
D0 = FIT{1};
D1 = FIT{2};
FIT = {D0, D1, D1 .* [q1 0; 0 q2], D1 .* [(1-q1) 0; 0 (1-q2)]};

% print per-class rates and covariance
% fai = mmap_count_mean(FIT,1);
% for i = 1:m
%     fprintf('Rate class %d: input = %.3f, output = %.3f\n', ...
%             i, ai(i), fai(i));
% end
% fsigma = 1/2*( map_count_var(FIT,t3) - sum(mmap_count_var(FIT,t3)) );
% fprintf('Covariance(t3): input = %.4f, output = %.4f\n', st3, fsigma);

end