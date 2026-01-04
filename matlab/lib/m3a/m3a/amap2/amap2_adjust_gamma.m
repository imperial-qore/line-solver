function [M2a, M3a, GAMMAa] = amap2_adjust_gamma(M1, M2, M3, GAMMA, weights, method, constraints)
% Computes a feasible set of characteristics for the MAP(2) that is as
% close as possible to the desired set of characterstics.
% Input:
% - M1,M2,M3: moments of the marginal distribution
% - GAMMA: auto-correlation decay rate
% - weights: an optional three-element vector with the weights associated
%            to M2, M3, GAMMA. Default is [10 1 10];
% - method: the method used to compute the feasible set of characteristics
%           1) pattern search on M2, M3, GAMMA
%           2) fit M3 as closely as possible, then perform pattern search
%              on M3, GAMMA (weight of M2 is ignored)
%           3) (default) prioritize M2, then M3, then GAMMA
%           4) fit M2 as closely as possible, then perform PSwarm
%              on M3 with bounds [M3_lb(M2), M3_ub(M2)] minimizing the 
%              objective function:
%              f(x) = weight(2)*(M3a/M3-1)^2 + weight(3)*(GAMMAa/GAMMA-1)^2
%              where GAMMAa = max(GAMMA_LB(M3a), min(GAMMA, GAMMA_UB(M3a))
% Output:
% - M2a,M3a: feasible moments of the marginal distribution (M1 is always
%            feasible)
% - GAMMAa:  feasible auto-correlation decay rate

if nargin == 4 || isempty(weights)
    % ignored with method 3 and 4
    weights = [10 1 10];
end

if nargin <= 5
    method = 3;
end
if nargin <= 6
    constraints = 2;
end

if constraints == 1
    if method == 1
        nonlcon = @nonlcon_safe;
    else
        nonlcon = @nonlcon_safe2;
    end
elseif constraints == 2
    if method == 1
        nonlcon = @nonlcon_theoretical;
    else
        nonlcon = @nonlcon_theoretical2;
    end
else
    error('Invalid option');
end

if method == 1 || method == 2
    options = psoptimset;
    options = psoptimset(options,'MaxIter', 1e5);
    options = psoptimset(options,'MaxFunEvals', 1e5);
    options = psoptimset(options,'PollingOrder', 'Success');
    options = psoptimset(options,'SearchMethod', {  @searchneldermead [] [] });
    options = psoptimset(options,'CompleteSearch', 'on');
    options = psoptimset(options,'Display', 'iter');
    options = psoptimset(options,'TolCon', 1e-100);
end

% tolerance for strict inequalities
tol = 1e-2;

if method == 1
    if GAMMA > 0
        FEASIBLE = map_scale(amap2_assemble(1,1/3,1/2,2/3,1),M1);
    else
        FEASIBLE = map_scale(amap2_assemble(1,2/3,1/2,2/3,2),M1);
    end
    FM2 = map_moment(FEASIBLE,2);
    FM3 = map_moment(FEASIBLE,3);
    FGAMMA = map_acf(FEASIBLE,4)/map_acf(FEASIBLE,3);
    x = patternsearch(@fun,[FM2 FM3 FGAMMA],[],[],[],[],[0 0 -1],[inf inf, 1-tol], nonlcon, options);
    M2a = x(1);
    M3a = x(2);
    GAMMAa = x(3);
elseif method == 2
    % force feasibility of the second moment
    M2a = max(3/2*M1^2, M2);
    % find a feasible value of the third moment
    n2 = M2/(M1^2);
    p2 = 3*(n2-2)/(3*n2) * (-2*sqrt(3)/sqrt(12-6*n2) - 1);
    a2 = (n2-2)/(p2*(1-n2) + sqrt(p2^2 + (2*p2*(n2-2))));
    l2 = 3*(a2+1)/(a2*p2+1) - (6*a2)/(2+a2*p2*(2*a2+2));
    u2 = 6*(n2-1)/n2;
    if 3/2 <= n2 && n2 < 2
        FM3 = (l2 + (u2-l2)/2) * M1* M2a;
        M3_LB = l2 * M1 * M2a;
        M3_UB = u2 * M1 * M2a;
    else
        FM3 = (3/2) * M2a^2 / M1 + tol;
        M3_LB = FM3;
        M3_UB = inf;
    end
    % null decay rate is always feasible
    FGAMMA = 0;
    % find best feasible values of M3 and GAMMA
    x = patternsearch(@fun2,[FM3 FGAMMA],[],[],[],[],[M3_LB -1],[M3_UB, 1-tol], nonlcon, options);
    %x = fmincon(@fun2,[FM3 FGAMMA],[],[],[],[],[M3_LB -1],[M3_UB, 1], nonlcon);
    M3a = x(1);
    GAMMAa = x(2);
elseif method == 3
    % priorities are M2 > M3 > GAMMA
    [M2a,M3a] = aph2_adjust(M1, M2, M3);
    [lb, ub] = compute_gamma_bounds(M2a, M3a);
    GAMMAa = max(lb, min(GAMMA, ub));
elseif method == 4
    % priorities are M2 > GAMMA > M3
    % adjust second moment
    M1sq = M1^2;
    scv = ((M2-M1sq)/M1sq);
    if scv < 1/2
        % (M2 - M1^2)/M1^2 = 1/2
        % M2 - M1^2 = 1/2 M1^2
        % M2 = 3/2 M1^2
        M2a = 3/2 * M1^2;
        scva = ((M2a-M1sq)/M1sq);
    else
        M2a = M2;
        scva = scv;
    end
    % get bounds on third moment
    if scva <= 1
        M3_lb = 3*M1^3*(3*scva-1+sqrt(2)*(1-scva)^(3/2));
        M3_ub = 6*M1^3*scva;
    elseif scva > 1
        M3_lb = 3/2*M1^3*(1+scva)^2;
        M3_ub = inf; 
    end
    % compute M3a and GAMMAa
    if abs(M3_lb - M3_ub) < tol
       % exponential
       M3a = (M3_lb + M3_ub)/2;
       GAMMAa = 0;
    else
        % optimization problem formulation
        problem.Variables = 1;
        problem.LB = M3_lb + tol;
        problem.UB = M3_ub;
        problem.ObjFunction = @gamma_objective;
        M3a = PSwarm(problem);
        % compute adjusted gamma
        [lb, ub] = compute_gamma_bounds(M2a, M3a);
        GAMMAa = max(lb, min(GAMMA, ub));
    end
else
    error('Invalid method for adjusting AMAP(2) characteristics');
end

    function obj = gamma_objective(M3a)
        [lb,ub] = compute_gamma_bounds(M2a, M3a);
        GAMMAa = max(lb, min(GAMMA, ub));
        obj = weights(2) * (M3a/M3-1)^2 + weights(3) * (GAMMAa/GAMMA-1)^2;
    end

    function [lb, ub] = compute_gamma_bounds(M2a, M3a)
        N2a = M2a/(M1)^2;
        N3a = M3a/(M2a*M1);
        if N2a < 2
            lb = -(N2a*(N3a-6)+6)/(3*N2a-6);
            ub = -(2* (0.5*(N2a-2)+0.5*sqrt(N2a^2 - 2*N2a*N3a/3))^2 )/(N2a-2);
            ub = ub * (1 - tol);
        elseif N3a < 9 - 12/N2a
            lb = -(N2a*(N3a-6)+6)/(3*N2a-6);
            ub = 1-tol;
        else
            x1 = sqrt(N2a*(N2a*(18*N2a+N3a*(N3a-18)-27)+24*N3a));
            x2 = N2a*(N3a-9);
            lb = (x2-x1+12)/(x2+x1+12);
            ub = 1-tol;
        end
    end

    function obj = fun(x)
        v = (x - [M2 M3 GAMMA]) ./ [M2 M3 GAMMA] .* weights;
        obj = norm(v);
    end

    function obj = fun2(x)
        v = (x - [M3 GAMMA]) ./ [M3 GAMMA] .* weights(2:3);
        obj = norm(v);
    end

    function [c, ceq] = nonlcon_safe(x)
        xM2 = x(1);
        xM3 = x(2);
        xGAMMA = x(3);
        AMAPS = amap2_fit_decay(M1,xM2,xM3,xGAMMA);
        if isempty(AMAPS)
            c = 1;
        else
            c = 0;
        end
        ceq = [];
    end

    function [c, ceq] = nonlcon_safe2(x)
        xM3 = x(1);
        xGAMMA = x(2);
        AMAPS = amap2_fit_decay(M1,M2,xM3,xGAMMA);
        if isempty(AMAPS)
            c = 1;
        else
            c = 0;
        end
        ceq = [];
    end

    function [c, ceq] = nonlcon_theoretical(x)
        xM2 = x(1);
        xM3 = x(2);
        xGAMMA = x(3);
        
        % compute normalized moments
        n2 = xM2/(M1^2);
        n3 = xM3/(M1*xM2);

        %fprintf('M1 = %f, M2 = %f, M3 = %f, GAMMA = %f\n', M1, xM2, xM3, xGAMMA);
        %fprintf('n2 = %f, n3 = %f\n', n2, n3);

        % simplifying notations for moments bounds
        p2 = 3*(n2-2)/(3*n2) * (-2*sqrt(3)/sqrt(12-6*n2) - 1);
        a2 = (n2-2)/(p2*(1-n2) + sqrt(p2^2 + (2*p2*(n2-2))));
        l2 = 3*(a2+1)/(a2*p2+1) - (6*a2)/(2+a2*p2*(2*a2+2));
        u2 = 6*(n2-1)/n2;

        %if (M2-M1^2)/M1^2 == 1
        %    fprintf('Exponential\n');
        %end

        eps = 1e-8;
        
        c = zeros(5,1);
        % CONSTRAINT 1: 3/2 <= n
        c(1) = 3/2 - n2;
        if 3/2 <= n2 && n2 < 2
            % CONSTRAINT 2: l2 <= n3
            % CONSTRAINT 3: n3 <= u2
            c(2) = l2 - n3;
            c(3) = n3 - u2;
        elseif n2 > 2
            % CONSTRAINT 2: 3/2 n2 < n3
            % CONSTRAINT 3: disabled
            c(2) = 3/2 * n2 - n3 + eps;
            c(3) = 0;
        end
        
        % BOUNDS FOR AUTOCORRELATION DECAY
        lb1 = -(n2*(n3-6)+6)/(3*n2-6);
        ub1 = -(2*(1/2*(n2-2)+1/2*sqrt(n2^2-(2*n2*n3)/3))^2)/(n2-2);
        tmp1 = n2*(n3-9);
        tmp2 = sqrt(n2*(n2*(18*n2+n3*(n3-18)-27)+24*n3));
        lb2 = (tmp1-tmp2+12)/(tmp1+tmp2+12);
        %fprintf('LB1 = %f, LB2 = %f\n', lb1, lb2);
        if n2 < 2
            %fprintf('LB1 <= gamma: %d\n', lb1 <= xGAMMA);
            %fprintf('gamma <= UB: %d\n', xGAMMA <= ub1);
            c(4) = lb1 - xGAMMA;
            c(5) = xGAMMA - ub1;
        elseif n2 > 2 && n3 < 9 - 12/n2
            %fprintf('LB1 <= gamma: %d\n', lb1 <= xGAMMA);
            %fprintf('gamma <= 1: %d\n', xGAMMA <= 1);
            c(4) = lb1 - xGAMMA;
            c(5) = xGAMMA - 1;
        elseif n2 > 2 && n3 >= 9 - 12/n2
            %fprintf('LB2 <= gamma: %d\n', lb2 <= xGAMMA);
            %fprintf('gamma <= 1: %d\n', xGAMMA <= 1);
            c(4) = lb2 - xGAMMA;
            c(5) = xGAMMA - 1;
        end
              
        ceq = [];

    end

    function [c, ceq] = nonlcon_theoretical2(x)
    
        xM3 = x(1);
        xGAMMA = x(2);
        
        % compute normalized moments
        n2 = M2/(M1^2);
        n3 = xM3/(M1*M2);
        
        % M3 is always feasible because of the bound constraints
        
        % BOUNDS FOR AUTOCORRELATION DECAY
        lb1 = -(n2*(n3-6)+6)/(3*n2-6);
        ub1 = -(2*(1/2*(n2-2)+1/2*sqrt(n2^2-(2*n2*n3)/3))^2)/(n2-2);
        tmp1 = n2*(n3-9);
        tmp2 = sqrt(n2*(n2*(18*n2+n3*(n3-18)-27)+24*n3));
        lb2 = (tmp1-tmp2+12)/(tmp1+tmp2+12);
        %fprintf('LB1 = %f, LB2 = %f\n', lb1, lb2);
        if n2 < 2
            %fprintf('LB1 <= gamma: %d\n', lb1 <= xGAMMA);
            %fprintf('gamma <= UB: %d\n', xGAMMA <= ub1);
            c(4) = lb1 - xGAMMA;
            c(5) = xGAMMA - ub1;
        elseif n2 > 2 && n3 < 9 - 12/n2
            %fprintf('LB1 <= gamma: %d\n', lb1 <= xGAMMA);
            %fprintf('gamma <= 1: %d\n', xGAMMA <= 1);
            c(4) = lb1 - xGAMMA;
            c(5) = xGAMMA - 1;
        elseif n2 > 2 && n3 >= 9 - 12/n2
            %fprintf('LB2 <= gamma: %d\n', lb2 <= xGAMMA);
            %fprintf('gamma <= 1: %d\n', xGAMMA <= 1);
            c(4) = lb2 - xGAMMA;
            c(5) = xGAMMA - 1;
        end
              
        ceq = [];

    end

end
