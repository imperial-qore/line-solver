function [fit,m3pps] = m3pp22_interleave_fit_count(av, btv, binfv, stv, t)
% Fits L pairs of classes into a single MMAP, obtained by lumped
% superposition of L different M3PP[2] processes.
% INPUT
% - av:    Matrix of size Lx2 containing the per-class rates.
% - btv:   Vector of length L containing the IDC at resolution t for each
%          pairs of classes.
% - binfv: Vector of length L containing the asymptotic IDC for each pair
%          of classes.
% - stv:   Vector of length L containing the count covariance, at
%          resolution t, between each pair of classes.
% OUTPUT
% - fit:   The lumped superposition of the M3PP[2] processes.
% - m3pps: Cell-array of length L containing the M3PP[2] processes.

% number of pairs
L = size(av,1);

% vector of lower bounds for the upper off-diagonal elements
uv = zeros(L,1);
% vector of upper bounds for the upper off-diagonal elements
dv = zeros(L,1);

% compute bounds for the upper off-diagonal elment of each mmpp(2)
for i = 1:L
    a = sum(av(i,:));
    d = compute_d(btv(i), binfv(i), t);
    z = (binfv(i)-1)*d^3*a;
    u = d*z/(2*a^2*d^2+z);
    uv(i) = u;
    dv(i) = d;
end

% find feasible values for off diagonal elements
A = zeros(2*L,2*L);
b = zeros(2*L,1);
for i = 1:L
    base = 2*(i-1);
    for j = 1:L
        if j >= i
            A(base+1,j) =  1;
            A(base+2,j) = -1;
        end
    end
    b(base+1) =  dv(i)-1e-6;
    b(base+2) = -uv(i)-1e-6;
end
Aeq = zeros(L,2*L);
beq = zeros(L,1);
for i = 1:L
    for j = 1:L
        if j >= i
            Aeq(i,j) = 1;
        end
    end
    for j = 1:L
        if j <= i
            Aeq(i,j+L) = 1;
        end
    end
    beq(i) = dv(i);
end
f = zeros(2*L,1);
x = linprog(f, A, b, Aeq, beq, zeros(2*L,1), [], []);
r = zeros(2,L);
r(1,:) = x(1:L);
r(2,:) = x((L+1):(2*L));

% fit m3pp[2] processes
m3pps = cell(L,1);
for i = 1:L
    % fit underling mmpp(2)
    r1 = 0;
    r2 = 0;
    for j = 1:L
        if j >= i
            r1 = r1 + r(1,j);
        end
        if j <= i
            r2 = r2 + r(2,j);
        end
    end
    mmpp = cell(1,2);
    mmpp{1}(1,2) = r1;
    mmpp{1}(2,1) = r2;
    a = sum(av(i,:));
    d = r1 + r2;
    z = (binfv(i)-1)*d^3*a;
    delta = sqrt(z/(2*r1*r2));
    l2 = a - r2/d * delta;
    l1 = l2 + delta;
    mmpp{2}(1,1) = l1;
    mmpp{2}(2,2) = l2;
    drates = -sum(mmpp{1}+mmpp{2},2);
    for j = 1:2
        mmpp{1}(j,j) = drates(j);
    end
    % variance
    v = map_count_var(mmpp, t);
    fprintf('MMPP %d - Var(t): %d\n', i, v);
    % fit m3pp(2,2)
    m3pps{i} = m3pp22_fit_count_approx_cov_multiclass(mmpp, av(i,:)', stv(i), t);
end

% fit lumped
fit = m3pp2m_interleave(m3pps);

    % solve non-linear equation to fit underlying mmpp(2)
    function d = compute_d(bt1,binf,t1)
        if ~(binf>bt1 && bt1>1)
            error('No solution, infeasible IDC(t): IDC(%.2f) = %.3f, IDC(inf) = %.3f\n',...
                   t1, bt1, binf);
        end
        % d = r1 + r2
        c = (binf-1)/(binf-bt1);
        % according to WolframAlpha
        % the solution can be written as (ProductLog[-c e^(-c)]+c)/t1
        z = -c*exp(-c);
        w = fsolve(@(w) z-w*exp(w),  1, optimset('Display','none',...
                                                 'MaxIter',10000,...
                                                 'MaxFunEvals',10000,...
                                                 'TolFun', 1.0e-12,...
                                                 'TolX',1.0e-12));
        d = (w + c)/t1;
    end

end