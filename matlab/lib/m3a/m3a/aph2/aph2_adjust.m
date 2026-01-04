function [M2a, M3a] = aph2_adjust(M1, M2, M3, method)
% Finds feasible values of M2 and M3 for an APH(2) distribution as close as
% possible to the values provided as an input.
% Input:
%   M1, M2, M3: the first three moments of the distribution
%   method:     string with the method name, can be one of the following
%               simple (default): [Telek and Heindl, 2002]
%               opt_param: optimize in parameter space
%               opt_param_gads: optimize in parameter space (global)
%               opt_char: optimize in characteristics space
%               opt_char_gads: optimize in characteristics space (global)
% Output:
%   M2a:        adjusted value for M2
%   M3a:        adjusted value for M3

if nargin < 4
    method = 'simple';
end

% tolerance for strict inequalities
tol = 1e-4;

if strcmp(method, 'simple')
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
    if scva <= 1
        lb = 3*M1^3*(3*scva-1+sqrt(2)*(1-scva)^(3/2));
        ub = 6*M1^3*scva;
        if M3 < lb
            M3a = lb;
        elseif M3 > ub
            M3a = ub;
        else
            M3a = M3;
        end
    else
        lb = 3/2*M1^3*(1+scva)^2;
        if M3 <= lb
            M3a = lb * (1 + tol);
        else
            M3a = M3;
        end
    end
elseif strcmp(method, 'opt_char') || strcmp(method, 'opt_char_gads')
    % optimize in moment space
    if strcmp(method, 'opt_char_gads')
        prob1 = createOptimProblem('fmincon', ...
                                  'objective', @fun, ...
                                  'x0', [M2 M3], ...
                                  'lb', [0 0], ...
                                  'ub', [inf inf], ...
                                  'nonlcon', @nonlcon1);
        prob2 = createOptimProblem('fmincon', ...
                                  'objective', @fun, ...
                                  'x0', [M2 M3], ...
                                  'lb', [0 0], ...
                                  'ub', [inf inf], ...
                                  'nonlcon', @nonlcon2);
        % adjustment to make the first fit feasible
        x = run(GlobalSearch('Display','iter'), prob1);
        M2a1 = x(1);
        M3a1 = x(2);
        % adjustment to make the second fit feasible
        x = run(GlobalSearch('Display','iter'), prob2);
        M2a2 = x(1);
        M3a2 = x(2);
    else
        options = optimset('Display', 'iter','Algorithm','active-set');
        x = fmincon(@fun, [M2 M3], [], [], [], [], [0 0], [], @nonlcon1, options);
        M2a1 = x(1);
        M3a1 = x(2);
        x = fmincon(@fun, [M2 M3], [], [], [], [], [0 0], [], @nonlcon2, options);
        M2a2 = x(1);
        M3a2 = x(2);
    end
    % pick the smallest adjustment
    if norm([M2a1 M3a1] - [M2 M3]) < norm([M2a2 M3a2] - [M2 M3])
        M2a = M2a1;
        M3a = M3a1;
    else
        M2a = M2a2;
        M3a = M3a2;
    end
elseif strcmp(method, 'opt_param') || strcmp(method, 'opt_param_gads')
    % optimize in parameter space
    feastol = 1e-6;
    degentol = 1e-8; % set to zero to allow exponential
    x0 = [M1 1/2]; % initial solution that fits M1 exactly
    lb = [feastol degentol];
    ub = [inf 1-degentol];
    if strcmp(method, 'opt_param_gads')
        prob = createOptimProblem('fmincon', ...
                                 'objective', @fun2, ...
                                 'x0', x0, ...
                                 'lb', lb, ...
                                 'ub', ub, ...
                                 'nonlcon', @nonlcon3);
        x = run(GlobalSearch('Display','none'), prob);
    else
        options = optimset('Display', 'none','Algorithm','active-set');
        x = fmincon(@fun2, x0, [], [], [], [], lb, ub, @nonlcon3, options);
    end
    l2 = x(1);
    r1 = x(2);
    l1 = M1 - l2*r1;
    M2a = 2*l1^2 + 2*r1*l1*l2 + 2*r1*l2^2;
    M3a = 6*l1^3 + 6*r1*l1^2*l2 + 6*r1*l1*l2^2 + 6*r1*l2^3;
else
    error('Invalid method: %s', method);
end

    function obj = fun2(x)
       l2 = x(1);
       r1 = x(2);
       l1 = M1 - l2*r1;
       xM2 = 2*l1^2 + 2*r1*l1*l2 + 2*r1*l2^2;
       xM3 = 6*l1^3 + 6*r1*l1^2*l2 + 6*r1*l1*l2^2 + 6*r1*l2^3;
       obj = norm([M2 M3] - [xM2 xM3]);
    end

    function [c, ceq] = nonlcon3(x)
        l2 = x(1);
        r1 = x(2);
        c = l2*r1 - M1 + feastol;
        ceq = [];
    end

    % minimize euclidean distance to sample values of M2 and M3
    function obj = fun(x)
        obj = norm(x - [M2 M3]);
    end

    function [p1,l1,l2] = aph2_fit1(xM2,xM3)
        tmp0 = (8*M1^3*xM3)/3 - 3*M1^2*xM2^2 - 2*M1*xM2*xM3 + 2*xM2^3 + xM3^2/9;
        tmp1 = 3*sqrt(tmp0);
        tmp2 = xM3 - 3*M1*xM2;
        tmp3 = (6*xM2 - 12*M1^2);
        l1 = (tmp2 + tmp1)/tmp3;
        l2 = (tmp2 - tmp1)/tmp3;
        p1 = (M1 - l1)/l2;
    end

    function [p1,l1,l2] = aph2_fit2(xM2,xM3)
        tmp0 = (8*M1^3*xM3)/3 - 3*M1^2*xM2^2 - 2*M1*xM2*xM3 + 2*xM2^3 + xM3^2/9;
        tmp1 = 3*sqrt(tmp0);
        tmp2 = xM3 - 3*M1*xM2;
        tmp3 = (6*xM2 - 12*M1^2);
        l1 = (tmp2 - tmp1)/tmp3;
        l2 = (tmp2 + tmp1)/tmp3;
        p1 = (M1 - l1)/l2;
    end

    % check feasibility for the first fit
    function [c, ceq] = nonlcon1(x)
        xM2 = x(1);
        xM3 = x(2);
        [p1,l1,l2] = aph2_fit1(xM2,xM3);
        c = zeros(4,1);
        c(1) = -tmp0; % non-negative square root argument
        c(2) = -l1 + 1e-6; % non-negative l1
        c(3) = -l2 + 1e-6; % non-negative l2
        c(4) = -p1; % non-negative p1
        ceq = [];
    end

    % check feasibility for the second fit
    function [c, ceq] = nonlcon2(x)
        xM2 = x(1);
        xM3 = x(2);
        [p1,l1,l2] = aph2_fit2(xM2,xM3);
        c = zeros(4,1);
        c(1) = -tmp0; % non-negative square root argument
        c(2) = -l1 + 1e-6; % non-negative l1
        c(3) = -l2 + 1e-6; % non-negative l2
        c(4) = -p1; % non-negative p1
        ceq = [];
    end

end