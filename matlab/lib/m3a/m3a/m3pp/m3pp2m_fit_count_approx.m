function [FIT] = m3pp2m_fit_count_approx(a, bt1, bt2, binf, m3t2, ...
    t1, t2, ...
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

if abs(a - sum(ai)) > 1e-8
    error('Inconsistent per-class arrival rates.');
end

% number of classes
m = size(ai,1);

% fit underlying MMPP(2)
FIT = mmpp2_fit_count_approx(a, bt1, bt2, binf, m3t2, t1, t2);

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

q1i_ai = ((r1^4*t)/2 + (r2^4*t)/2 - l1*r2^3*t + l2*r2^3*t + 2*r1*r2^3*t + 2*r1^3*r2*t + 3*r1^2*r2^2*t + 2*l1*r2^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*r2^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l1*r1*r2^2*t - l1*r1^2*r2*t + 2*l2*r1*r2^2*t + l2*r1^2*r2*t + 2*l1*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/(l1*r2*(r1 + r2)*(2*l1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*r1*t - l1*r2*t + l2*r1*t + l2*r2*t));
q1i_dvi = -(r1^4 + 4*r1^3*r2 + 6*r1^2*r2^2 + 4*r1*r2^3 + r2^4)/(4*l1*r2*(r1 + r2)*(2*l1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*r1*t - l1*r2*t + l2*r1*t + l2*r2*t));
q1i_const = -(l1*r2^4*t + l2*r1^4*t + 3*l1*r1*r2^3*t + l1*r1^3*r2*t + l2*r1*r2^3*t + 3*l2*r1^3*r2*t + 3*l1*r1^2*r2^2*t + 2*l1^2*r1*r2^2*t + 2*l1^2*r1^2*r2*t + 3*l2*r1^2*r2^2*t + 2*l2^2*r1*r2^2*t + 2*l2^2*r1^2*r2*t - 4*l1^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l2^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l1*l2*r1*r2^2*t - 4*l1*l2*r1^2*r2*t + 8*l1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/(4*l1*r2*(r1 + r2)*(2*l1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 2*l2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*r1*t - l1*r2*t + l2*r1*t + l2*r2*t));
q2i_ai = -((r1^4*t)/2 + (r2^4*t)/2 + l1*r1^3*t - l2*r1^3*t + 2*r1*r2^3*t + 2*r1^3*r2*t + 3*r1^2*r2^2*t - 2*l1*r1^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 2*l2*r1^2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + l1*r1*r2^2*t + 2*l1*r1^2*r2*t - l2*r1*r2^2*t - 2*l2*r1^2*r2*t - 2*l1*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) + 2*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/((r1 + r2)*(l2^2*r1^2*t - 2*l2^2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1^2*t + l2^2*r1*r2*t + 2*l1*l2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1*r2*t));
q2i_dvi = (r1^4 + 4*r1^3*r2 + 6*r1^2*r2^2 + 4*r1*r2^3 + r2^4)/(4*(r1 + r2)*(l2^2*r1^2*t - 2*l2^2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1^2*t + l2^2*r1*r2*t + 2*l1*l2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1*r2*t));
q2i_const = (l1*r2^4*t + l2*r1^4*t + 3*l1*r1*r2^3*t + l1*r1^3*r2*t + l2*r1*r2^3*t + 3*l2*r1^3*r2*t + 3*l1*r1^2*r2^2*t + 2*l1^2*r1*r2^2*t + 2*l1^2*r1^2*r2*t + 3*l2*r1^2*r2^2*t + 2*l2^2*r1*r2^2*t + 2*l2^2*r1^2*r2*t - 4*l1^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l2^2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - 4*l1*l2*r1*r2^2*t - 4*l1*l2*r1^2*r2*t + 8*l1*l2*r1*r2*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2))/(4*(r1 + r2)*(l2^2*r1^2*t - 2*l2^2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1^2*t + l2^2*r1*r2*t + 2*l1*l2*r1*sinh((r1*t)/2 + (r2*t)/2)*exp(- (r1*t)/2 - (r2*t)/2) - l1*l2*r1*r2*t));

H = zeros(m, m);
f = zeros(m, 1);
for i = 1:m
    H(i,i) = 2/dvt3(i)^2;
    f(i) = -2/dvt3(i);
end

A = zeros(2*m,m);
b = zeros(2*m,1);
for i = 1:m
    A(1+2*(i-1),i) = -q1i_dvi;
    b(1+2*(i-1))   =  q1i_const + q1i_ai * ai(i);
    A(2+2*(i-1),i) = -q2i_dvi;
    b(2+2*(i-1))   =  q2i_const + q2i_ai * ai(i);
end

Aeq = zeros(2,m);
beq = zeros(2,1);
for i = 1:m
    Aeq(1,i) = q1i_dvi;
    Aeq(2,i) = q2i_dvi;
end
beq(1) = 1 - m*q1i_const - q1i_ai * a;
beq(2) = 1 - m*q2i_const - q2i_ai * a;

%fprintf('Fitting per-class counting process...\n');
options = optimset('Algorithm','interior-point-convex ',...
    'Display','none');
%[x,fx] = quadprog(H, f, A, b, Aeq, beq, [], [], [], options);
lb = 1e-6*ones( size(A,2),1);
ub = 1e6*ones( size(A,2),1);
[x,fx]=QP(H, h, A, b, Aeq, beq, lb, ub, options);
%fit_error = fx + m;
%fprintf('Per-class fitting error: %f\n', fit_error);

q = zeros(2,m);
for i = 1:m
    Ai = ai(i); % rate fitted exactly
    Dvi = x(i); % result of optimization (optimal)
    q(1,i) = q1i_ai * Ai + q1i_dvi * Dvi + q1i_const;
    q(2,i) = q2i_ai * Ai + q2i_dvi * Dvi + q2i_const;
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
    error('Infeasible fitted M3PP');
end

Ai = mmap_count_mean(FIT,1);
Dvt3 = zeros(m,1);
for i = 1:m
    mmap2 = {FIT{1},FIT{2},FIT{2+i},FIT{2}-FIT{2+i}};
    v2t3 = mmap_count_var(mmap2,t3);
    Dvt3(i) = v2t3(1)-v2t3(2);
end
%
% for i = 1:m
%     fprintf('Rate class %d: input = %.3f, output = %.3f\n', ...
%             i, ai(i), Ai(i));
% end
% for i = 1:m
%     fprintf('DV%d(t3): input = %.3f, output = %.3f\n', ...
%             i, dvt3(i), Dvt3(i));
% end

end