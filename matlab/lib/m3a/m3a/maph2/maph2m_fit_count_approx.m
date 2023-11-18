function [FIT] = maph2m_fit_count_approx(a, bt1, binf, t1, ...
                                         ai, dvt3, t3)
% Fits a second-order Marked MMPP.
% a: arrival rate
% bt1: IDC at scale t1
% bt2: IDC at scale t2
% binf: IDC for t->inf
% t1: first time scale
% t2: second time scale
% ai: i-th element is the rate of class i
% dvt3: i-th element is the delta of variance of class i and the variance
%       of all other classes combined, at resolution t3
% t3: third time scale

method_d0d1 = 'exact'; 

options = struct('Algorithm', 'active-set', ...
                 'Display', 'none', ...
                 'MaxFunEvals', 100000);
% l1, l2, r1 (rate is fitted exactly)
x0 = [a/2, 2*a, a];
lb = [1e-6, 1e-6, 1e-6];
ub = [inf, inf, inf];
fprintf('Fitting unified counting process...\n');
if strcmp(method_d0d1,'exact')
    FIT = aph2_fit_count(a, bt1, binf, t1);
elseif strcmp(method_d0d1,'fmincon')
    xopt = fmincon(@compute_obj, ...
                   x0, ...
                   [], [], ... % linear inequalities
                   [], [], ... % linear equalities
                   lb, ub, ... % bounds
                   [], options);
    FIT = assemble_mmap(xopt);
elseif strcmp(method_d0d1,'gs_fmincon')
    problem = createOptimProblem('fmincon', ...
                                 'objective', @compute_obj, ...
                                 'x0', x0, ...
                                 'lb', lb, ...
                                 'ub', ub, ...
                                 'options', options);
    xopt = run(GlobalSearch,problem);
    FIT = assemble_mmap(xopt);
else
    error('Invalid method ''%s''\n', method);
end

FIT{1}
fa = map_count_mean(FIT,1);
fbt1 = map_count_var(FIT,t1)/map_count_mean(FIT,t1);
fbinf = map_count_var(FIT,1e8)/map_count_mean(FIT,1e8);
ferror = (fa/a-1)^2 + (fbt1/bt1-1)^2 + (fbinf/binf-1)^2;
fprintf('Unified fitting error: %f\n', ferror);
fprintf('Rate: input = %.3f, output = %.3f\n', a, fa);
fprintf('IDC(t1): input = %.3f, output = %.3f\n', bt1, fbt1);
fprintf('IDC(inf): input = %.3f, output = %.3f\n', binf, fbinf);

FIT = map_scale(FIT, 1/a);
fa = map_count_mean(FIT,1);
fbt1 = map_count_var(FIT,t1)/map_count_mean(FIT,t1);
fbinf = map_count_var(FIT,1e8)/map_count_mean(FIT,1e8);
fprintf('Rate: input = %.3f, output = %.3f\n', a, fa);
fprintf('IDC(t1): input = %.3f, output = %.3f\n', bt1, fbt1);
fprintf('IDC(inf): input = %.3f, output = %.3f\n', binf, fbinf);

l1 = FIT{2}(1,1);
l2 = FIT{2}(2,1);
r1 = FIT{1}(1,2);
t = t3;

m = size(ai,1);

q1i_ai = -(l1*r1 - l1*l2 - 2*l2*r1 + l1*l2*exp(l2*t + r1*t) - l1*r1*exp(l2*t + r1*t) + 2*l2*r1*exp(l2*t + r1*t) + l2^3*t*exp(l2*t + r1*t) + r1^3*t*exp(l2*t + r1*t) - l1*l2^2*t*exp(l2*t + r1*t) + l1*r1^2*t*exp(l2*t + r1*t) + l2*r1^2*t*exp(l2*t + r1*t) + l2^2*r1*t*exp(l2*t + r1*t))/(l1*l2*(l1 + r1)*(l2*t*exp(l2*t + r1*t) - exp(l2*t + r1*t) + r1*t*exp(l2*t + r1*t) + 1));
q1i_dvi = (l2^4 + 4*l2^3*r1 + 6*l2^2*r1^2 + 4*l2*r1^3 + r1^4)/(2*l1*l2*(l2 + r1)*(r1^2*t - 2*l1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) - 2*r1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) + l1*l2*t + l1*r1*t + l2*r1*t));
q1i_const = (l2*r1 - l1*r1 + (l2^3*t)/2 + (r1^3*t)/2 + l1*r1^2*t + (l2*r1^2*t)/2 + (l2^2*r1*t)/2 + l1*l2*r1*t)/(l1*(l2 + r1)*(l2*t + r1*t - 1)) - ((l2*r1 - l1*r1 + (l2^3*t)/2 + (r1^3*t)/2 + l1*r1^2*t + (l2*r1^2*t)/2 + (l2^2*r1*t)/2 + l1*l2*r1*t)/(l2*t + r1*t - 1) - l1*r1 + l2*r1)/(l1*(l2 + r1)*(l2*t*exp(l2*t + r1*t) - exp(l2*t + r1*t) + r1*t*exp(l2*t + r1*t) + 1));
q2i_ai = (l2^4*t + 2*r1^4*t + 2*l1*r1^3*t + 5*l2*r1^3*t + 3*l2^3*r1*t + 5*l2^2*r1^2*t - 2*r1^3*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) - 4*l1*r1^2*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) + 2*l2^2*r1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) + 4*l1*l2*r1^2*t + 2*l1*l2^2*r1*t - 4*l1*l2*r1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2))/((l2 + r1)*(l2*r1^3*t + l2^2*r1^2*t - 2*l2*r1^2*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) + l1*l2*r1^2*t + l1*l2^2*r1*t - 2*l1*l2*r1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2)));
q2i_dvi = -(l2^4 + 4*l2^3*r1 + 6*l2^2*r1^2 + 4*l2*r1^3 + r1^4)/(2*(l2 + r1)*(l2*r1^3*t + l2^2*r1^2*t - 2*l2*r1^2*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2) + l1*l2*r1^2*t + l1*l2^2*r1*t - 2*l1*l2*r1*sinh((l2*t)/2 + (r1*t)/2)*exp(- (l2*t)/2 - (r1*t)/2)));
q2i_const = ((l2*r1 - l1*r1 + (l2^3*t)/2 + (r1^3*t)/2 + l1*r1^2*t + (l2*r1^2*t)/2 + (l2^2*r1*t)/2 + l1*l2*r1*t)/(l2*t + r1*t - 1) - l1*r1 + l2*r1)/(r1*(l2 + r1)*(l2*t*exp(l2*t + r1*t) - exp(l2*t + r1*t) + r1*t*exp(l2*t + r1*t) + 1)) - (l2*r1 - l1*r1 + (l2^3*t)/2 + (r1^3*t)/2 + l1*r1^2*t + (l2*r1^2*t)/2 + (l2^2*r1*t)/2 + l1*l2*r1*t)/(r1*(l2 + r1)*(l2*t + r1*t - 1));

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
for i = 1:(m-1)
    Aeq(1,i) = q1i_dvi;
    Aeq(2,i) = q2i_dvi;
end
Aeq(1,m) = q1i_dvi;
Aeq(2,m) = q2i_dvi;
beq(1) = 1 - m*q1i_const - q1i_ai * a;
beq(2) = 1 - m*q2i_const - q2i_ai * a;

fprintf('Fitting per-class counting process...\n');
options = optimset('Algorithm','interior-point-convex ',...
                   'Display','none');
[x,fx] = quadprog(H, f, A, b, Aeq, beq, [], [], [], options);
fit_error = fx + m;
fprintf('Per-class fitting error: %f\n', fit_error);

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
    FIT{2+i} = FIT{2} .* [q(1,i) 0; q(2,i) 0];
end

Ai = mmap_count_mean(FIT,1);
Dvt3 = zeros(m,1);
for i = 1:m
    mmap2 = {FIT{1},FIT{2},FIT{2+i},FIT{2}-FIT{2+i}};
    v2t3 = mmap_count_var(mmap2,t3);
    Dvt3(i) = v2t3(1)-v2t3(2);
end

for i = 1:m
    fprintf('Rate class %d: input = %.3f, output = %.3f\n', ...
            i, ai(i), Ai(i));
end
for i = 1:m
    fprintf('DV%d(t3): input = %.3f, output = %.3f\n', ...
            i, dvt3(i), Dvt3(i));
end

    function obj = compute_obj(x)
        % compute characteristics
        xa = compute_rate(x);
        factor = a/xa;
        xbt1 = compute_idc(x, t1*factor);
        xbinf = compute_idc_limit(x);
        % compute objective
        obj = 0;
        obj = obj + (xa/a-1)^2;
        obj = obj + (xbt1/bt1-1)^2;
        obj = obj + (xbinf/binf-1)^2;
    end

    function xa = compute_rate(x)
        l1 = x(1);
        l2 = x(2);
        r1 = x(3);
        xa = (l2*(l1 + r1))/(l2 + r1);
    end

    function xbt = compute_idc(x,t)
        l1 = x(1);
        l2 = x(2);
        r1 = x(3);
        xbt = l2^3/(l2 + r1)^3 + r1^3/(l2 + r1)^3 + (2*l1*r1^2)/(l2 + r1)^3 + (l2*r1^2)/(l2 + r1)^3 + (l2^2*r1)/(l2 + r1)^3 + (2*l1*l2*r1)/(l2 + r1)^3 - (2*l1*r1)/(t*(l2 + r1)^3) + (2*l2*r1)/(t*(l2 + r1)^3) + (2*l1*r1*exp(- l2*t - r1*t))/(t*(l2 + r1)^3) - (2*l2*r1*exp(- l2*t - r1*t))/(t*(l2 + r1)^3);
    end

    function xbinf = compute_idc_limit(x)
        l1 = x(1);
        l2 = x(2);
        r1 = x(3);
        xbinf = (r1*(2*l1 - 2*l2))/(l2 + r1)^2 + 1;
    end

    function mmap = assemble_mmap(x)
        % extract parameters
        l1 = x(1);
        l2 = x(2);
        r1 = x(3);
        % assemble
        D0 = [-(l1+r1), r1; 0, -l2];
        D1 = [l1 0; l2 0];
        mmap = {D0, D1};
    end

end