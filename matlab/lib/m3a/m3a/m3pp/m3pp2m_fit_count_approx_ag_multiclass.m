function [FIT] = m3pp2m_fit_count_approx_ag_multiclass(mmpp, ai, gt3, t3)
% Fits a M3PP(2,m) given the underlying MMPP(2).
% INPUT
% - mmpp: underlying mmpp(2)
% - ai: i-th element is the rate of class i
% - gt3: i-th element is the sum of the variance of class i and its
%        covariance with all the other classes combined.
% - t3: third time scale

% number of classes
m = length(ai);

FIT = mmpp;

% total rate
a = map_count_mean(mmpp,1);

% check per-class rates consistency
if abs(a - sum(ai)) > 1e-8
    error('Inconsistent per-class arrival rates.');
end

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

f1num = l1*r2*(2*l2*r1 - 2*l1*r1 + r1^3*t + r2^3*t + 2*l1*r1^2*t - 2*l2*r1^2*t + 3*r1*r2^2*t + 3*r1^2*r2*t + 2*l1*r1*exp(- r1*t - r2*t) - 2*l2*r1*exp(- r1*t - r2*t) + 2*l1*r1*r2*t - 2*l2*r1*r2*t);
f1 = f1num/(r1+r2)^4;
f2num = l2*r1*(2*l1*r2 - 2*l2*r2 + r1^3*t + r2^3*t - 2*l1*r2^2*t + 2*l2*r2^2*t + 3*r1*r2^2*t + 3*r1^2*r2*t - 2*l1*r2*exp(- r1*t - r2*t) + 2*l2*r2*exp(- r1*t - r2*t) - 2*l1*r1*r2*t + 2*l2*r1*r2*t);
f2 = f2num/(r1+r2)^4;
tmp = (f1*l2*r1 - f2*l1*r2);
q1i_ai = -(f2*(r1 + r2))/tmp;
q1i_gi = (l2*r1)/tmp;
q1i_const = 0;
q2i_ai = (f1*(r1 + r2))/tmp;
q2i_gi = -(l1*r2)/tmp;
q2i_const = 0;

%%

H = zeros(m, m);
f = zeros(m, 1);
for i = 1:m
    H(i,i) = 2/gt3(i)^2;
    f(i) = -2/gt3(i);
end

A = zeros(2*m,m);
b = zeros(2*m,1);
for i = 1:m
    A(1+2*(i-1),i) = -q1i_gi;
    b(1+2*(i-1))   =  q1i_const + q1i_ai * ai(i);
    A(2+2*(i-1),i) = -q2i_gi;
    b(2+2*(i-1))   =  q2i_const + q2i_ai * ai(i);
end

Aeq = zeros(2,m);
beq = zeros(2,1);
for i = 1:m
    Aeq(1,i) = q1i_gi;
    Aeq(2,i) = q2i_gi;
end
beq(1) = 1 - m*q1i_const - q1i_ai * a;
beq(2) = 1 - m*q2i_const - q2i_ai * a;

fprintf('Fitting per-class counting process...\n');
options = optimset('Algorithm','interior-point-convex ',...
    'Display','none');
f(~isfinite(f))=0; % remove infnans
H(~isfinite(H))=0;
A(~isfinite(A))=0;
b(~isfinite(b))=0;
Aeq(~isfinite(Aeq))=0;
beq(~isfinite(beq))=0;
[x,fx] = quadprog(H, f, A, b, Aeq, beq, [], [], [], options);
lb = zeros( size(A,2),1);
ub = 1e6*ones( size(A,2),1);
%[x,fx]=QP(H, f, A, b, Aeq, beq, lb, ub, options)
fit_error = fx + m;
fprintf('Per-class fitting error: %f\n', fit_error);

q = zeros(2,m);
for i = 1:m
    Ai = ai(i); % rate fitted exactly
    Gi = x(i); % result of optimization (optimal)
    q(1,i) = q1i_ai * Ai + q1i_gi * Gi + q1i_const;
    q(2,i) = q2i_ai * Ai + q2i_gi * Gi + q2i_const;
end

D0 = FIT{1};
D1 = FIT{2};
FIT = cell(1,2+m);
FIT{1} = D0;
FIT{2} = D1;
for i = 1:m
    FIT{2+i} = FIT{2} .* [q(1,i) 0; 0 q(2,i)];
end

if ~mmap_isfeasible(FIT)
    FIT{:}
    warning('Infeasible fitted M3PP');
end

Ai = mmap_count_mean(FIT,1);
Gt3 = zeros(m,1);
for i = 1:m
    mmap2 = {FIT{1},FIT{2},FIT{2+i},FIT{2}-FIT{2+i}};
    v2t3 = mmap_count_var(mmap2,t3);
    s2t3 = mmap_count_mcov(mmap2,t3);
    Gt3(i) = v2t3(1) + s2t3(1,2);
end

for i = 1:m
    fprintf('Rate class %d: input = %.4f, output = %.4f\n', ...
        i, ai(i), Ai(i));
end
for i = 1:m
    fprintf('g%d(t3): input = %.4f, output = %.4f\n', ...
        i, gt3(i), Gt3(i));
end

end