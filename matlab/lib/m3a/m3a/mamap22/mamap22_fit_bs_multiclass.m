function [mmap,fB,fS,exact] = mamap22_fit_bs_multiclass(map,p,B,S,classWeights,bsWeights,adjust)
% Performs approximate fitting of a MMAP given the underlying AMAP(2),
% the class probabilities (always fitted exactly), the backward moments,
% and the one-step class transition probabilities.
% Input
% - map:  second-order AMAP underlying the MAMAP[2]
% - p:    vector of class probabilities
% - B:    vector of backward moments
% - S:    matrix of the class transition probabilities
% - classWeights: optional vector of weights for each class
% - bsWeights:    optional 2-vector of weights of backward moments and the
%                 class transition probabilities
% Output
% - mmap: fitted MAMAP[2]
% - fB:   vector of optimal feasible backward moments
% - fS:   vector of optimal feasible class transition probabilities

if (size(map{1},1) ~= 2)
    error('Underlying MAP must be of second-order.');
end
if (map{1}(2,1) ~= 0)
    error('Underlying MAP must be acyclic');
end
if (map{2}(1,2) == 0)
    form = 1;
elseif (map{2}(1,1) == 0)
    form = 2;
else
    error('Underlying MAP must be in canonical acyclic form');
end

fprintf('Fitting MAMAP(2,2) B+S: form = %d\n', form);

% number of classes
k = length(p);
if k ~= 2
    error('Fitting MAMAP(2,2) B+S: fitting backward and transition probabilities only supports two classes');
end

% default weights to use in the objective function
if nargin < 5 || isempty(classWeights)
    classWeights = ones(k,1);
end
if nargin < 6 || isempty(bsWeights)
    bsWeights = ones(2,1);
end
if nargin < 7
    adjust = 1;
end

% result
mmap = cell(1,2+k);
mmap{1} = map{1};
mmap{2} = map{2};

% get parameters of the underlying AMAP(2)
h1 = -1/map{1}(1,1);
h2 = -1/map{1}(2,2);
r1 = map{1}(1,2) * h1;
r2 = map{2}(2,2) * h2;

% set tolerance constants
degentol = 1e-6;
feastol  = 1e-4;
denumtol = 1e-12;

% generate required coefficients
if form == 1
    [G,U,Y] = mamap2m_can1_coefficients(h1,h2,r1,r2);
    fit = @fit_can1;
elseif form == 2
    [E,V,Z] = mamap2m_can2_coefficients(h1,h2,r1,r2);
    fit = @fit_can2;
end

exact = 0;

if (form == 1 && (r1 < degentol || r2 > 1-degentol || abs(h2-h1*r2) < degentol )) || ...
        (form == 2 && (r2 > 1-degentol || abs(h1 - h2 - h1*r1 + h1*r1*r2) < degentol ))
    
    % TODO: transform into a valid second-order Poisson process
    
    % DEGENERATE PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) B+S: detected Poisson process\n');
    
    % return marked poisson process
    h = map_mean(mmap);
    mmap = cell(1,2+k);
    mmap{1} = -1/h;
    mmap{2} =  1/h;
    for c = 1:k
        mmap{2+c} = mmap{2} * p(c);
    end
    
    fB = mmap_backward_moment(mmap,1);
    fS = mmap_sigma(mmap);
    
    return;
    
elseif form == 2 && r2 < degentol && abs(1-r1) < degentol
    
    % DEGENERATE PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) B+S: detected degenerate phase-type form\n');
    
    % only one degree of freedom: match class probabilites
    q1 = p(1);
    q2 = p(1);
    q3 = p(1);
    
elseif form == 1 && r2 < degentol
    
    % CANONICAL PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) B+S: detected canonical phase-type form\n');
    
    % convert to phase-type
    aph = map;
    aph{2}(2,2) = 0;
    aph = map_normalize(aph);
    
    mmap = maph2m_fit_multiclass(aph, p, B, classWeights);
    
    fB = mmap_backward_moment(mmap, 1);
    fS = mmap_sigma(mmap);
    
    return;
    
elseif (form == 1 && abs(1-r1) < degentol) || ...
        (form == 2 && abs(1-r1) < degentol)
    
    % NON-CANONICAL PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) B+S: detected non-canonical phase-type form, converting to canonical form\n');
    
    aph = aph2_fit_map(map);
    
    mmap = maph2m_fit_multiclass(aph, p, B, classWeights);
    
    fB = mmap_backward_moment(mmap, 1);
    fS = mmap_sigma(mmap);
    
    return;
    
elseif form == 2 && r2 < degentol
    
    % DEGENERATE CASE FOR gamma < 0
    fprintf('Fitting MAMAP(2,2) B+S: detected degenerate MMAP form\n');
    
    if bsWeights(1) > bsWeights(2)
        fprintf('Fitting MAMAP(2,2) B+S: fitting B\n');
        [q1,q2,q3] = fit_can2_degen_backward(B(1));
        if ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
            % four inequalities, one variable
            A = zeros(4,1);
            b = zeros(4,1);
            % coefficients
            q1_B = - p(1) * (r1-2)/((r1-1)*(h2-h1+h1*r1));
            q1_0 = + p(1) * (r1-2)*(h2+h1*r1)/((r1-1)*(h2-h1+h1*r1));
            q2_B = - p(1) * (r1-2)/(h2-h1+h1*r1);
            q2_0 = + p(1) * (r1-2)*h1/(h2-h1+h1*r1);
            % q1 >= 0
            A(1,1) = -q1_B;
            b(1)   =  q1_0;
            % q1 <= 1
            A(2,1) = +q1_B;
            b(2)   =  1 - q1_0;
            % q2 >= 0
            A(3,1) = -q2_B;
            b(3)   =  q2_0;
            % q2 <= 1
            A(4,1) = +q2_B;
            b(4)   =  1 - q2_0;
            % objective function
            H =  2 / B(1)^2;
            h = -2 / B(1);
            % solve
            fB1 = solve_quadprog();
            % compute coefficient
            [q1,q2,q3] = fit_can2_degen_backward(fB1);
        end
    else
        fprintf('Fitting MAMAP(2,2) B+S: fitting S\n');
        [q1,q2,q3] = fit_can2_degen_transition(S(1,1));
        if ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
            % bound constraints on S11
            safetytol = 1e-10;
            q1lb = p(1)^2 * (1 - (1-r1)^2);
            q2ub = p(1)^2 - (1 - p(1))^2;
            if S(1,1) <= q1lb
                fS11 = q1lb + safetytol;
            elseif S(1,1) >= q2ub
                fS11 = q2ub - safetytol;
            else
                fS11 = S(1,1);
            end
            % compute parameters
            [q1,q2,q3] = fit_can2_degen_transition(fS11);
        end
    end
    
else
    
    % FULL FORM or "GOOD" poisson process
    
    if (form == 1 && abs(h1 - h2 + h2*r1) < degentol) || (form == 2 && abs(h1 - h2 + h2*r1) < degentol)
        fprintf('Fitting MAMAP(2,2) B+S: detected fittable Poisson process\n')
    end
    
    [q1,q2,q3] = fit(B(1), S(1,1));
    
    % options used to solve nonlinear problem
    tol = 1e-6;
    
    yalmip_nonlinear_opt = sdpsettings(...
        'verbose',2,...
        'solver','bmibnb',...
        'debug',0,...
        'showprogress',0,...
        'usex0',1,...
        'bmibnb.relgaptol',1e-3,...
        'bmibnb.absgaptol',1e-6,...
        'bmibnb.maxiter',1000,...
        'bmibnb.lowrank', 1,...
        'bmibnb.lpreduce',1,... % without this, crappy lower bounds for some problems
        'bmibnb.pdtol',-1e-8,... % x >= if x > -pdtol
        'bmibnb.eqtol',+1e-10,... % x == 0 if abs(x) < +eqtol
        'bmibnb.roottight',1,...
        'bmibnb.uppersolver','fmincon-standard',...
        'fmincon.Algorithm','sqp',...
        'fmincon.TolCon',1e-6,...
        'fmincon.TolX',1e-6,...
        'fmincon.GradObj','on',...
        'fmincon.GradConstr','on',...
        'fmincon.Hessian','off',...
        'fmincon.LargeScale','off');
    
    % constants
    M1 = map_mean(map);
    
    if isfeasible(q1) && isfeasible(q2) && isfeasible(q3)
        fprintf('Fitting MAMAP(2,2) B+S: exact fit found\n');
        exact = 1;
    elseif adjust
        [q1,q2,q3] = solve_nonlinear();
    end
    
end % end if/then/else for degenerate forms

% check parameter feasibility
if adjust && ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
    error('Fitting MAMAP(2,2) B+S: Feasibility could not be restored');
end
% parameters feasible within feastol: restrict to [0,1]
q1 = fix(q1);
q2 = fix(q2);
q3 = fix(q3);

% compute D11,D12,...,D1k
if form == 1
    mmap{2+1} = mmap{2} .* [q1 0; q2 q3];
    mmap{2+2} = mmap{2} .* [1-q1 0; 1-q2 1-q3];
else
    mmap{2+1} = mmap{2} .* [0 q1; q2 q3];
    mmap{2+2} = mmap{2} .* [0 1-q1; 1-q2 1-q3];
end

fB = mmap_backward_moment(mmap, 1);
fS = mmap_sigma(mmap);

    function feas = isfeasible(q)
        feas = q >= -feastol && q <= (1+feastol);
    end

    function qfix = fix(q)
        qfix = max(min(q,1),0);
    end

    function [q1,q2,q3] = fit_can1(vB1,vS11)
        denum = U(11) * vB1 * p(1) + U(12) * p(1);
        if abs(denum) < denumtol
            q1 = p(1);
            q2 = p(1);
            q3 = p(1);
        else
            q2 = (U(7) * vB1^2 * p(1)^2 + U(8) * vB1 * p(1)^2 + U(9) * vS11 + U(10) * p(1)^2)/denum;
            q1 = -(G(12)*p(1) - vB1*G(3)*p(1) + (G(3)*G(11) - G(2)*G(12))*q2)/Y(3);
            q3 = +(G(10)*p(1) - vB1*G(1)*p(1) + (G(1)*G(11) - G(2)*G(10))*q2)/Y(3);
        end
    end

    function [q1,q2,q3] = fit_can2(vB1,vS11)
        denum = (V(11) * vB1 * p(1) + V(12) * p(1));
        if abs(denum) < denumtol
            q1 = p(1);
            q2 = p(1);
            q3 = p(1);
        else
            q3 = (V(7) * vB1^2 * p(1)^2 + V(8) * vB1 * p(1)^2 + V(9) * p(1)^2 + V(10) * vS11)/denum;
            q1 = +(E(10)*p(1) - vB1*E(2)*p(1) + (E(2)*E(11) - E(3)*E(10))*q3)/Z(3);
            q2 = -(E(9)*p(1) - vB1*E(1)*p(1) + (E(1)*E(11) - E(3)*E(9))*q3)/Z(3);
        end
    end

    function [q1,q2,q3] = fit_can2_degen_backward(vB1)
        q1 = (p(1)*(r1 - 2)*(h2 - vB1 + h1*r1))/((r1 - 1)*(h2 - h1 + h1*r1));
        q2 = -(p(1)*(vB1 - h1)*(r1 - 2))/(h2 - h1 + h1*r1);
        q3 = p(1); % meangingless if really degenerate
    end

    function [q1,q2,q3] = fit_can2_degen_transition(vS11)
        q1 = p(1) + (p(1)^2 - vS11)^(1/2)/(r1 - 1);
        q2 = p(1) + (p(1)^2 - vS11)^(1/2);
        q3 = p(1); % meangingless if really degenerate
    end

    function [q1,q2,q3,obj] = solve_nonlinear()
        
        fprintf('Fitting MAMAP(2,2) B+S: running approximate fitting (nonlinear)\n');
        
        % for now we support only one of the formulations
        solve_nonlinear_func = @solve_nonlinear_hybrid;
        
        % trivial solution
        eobj = make_objective(M1, p(1)^2);
        fprintf('Fitting MAMAP(2,2) B+S: B1 == M1, objective = %e\n', eobj);
        % left solution
        [lfeas,lobj,lF1,lS11] = solve_nonlinear_func('left');
        if lfeas
            fprintf('Fitting MAMAP(2,2) B+S: B1 < M1, objective = %e\n', lobj);
        else
            fprintf('Fitting MAMAP(2,2) B+S: B1 < M1 infeasible\n');
        end
        % right solution
        [rfeas,robj,rF1,rS11] = solve_nonlinear_func('right');
        if rfeas
            fprintf('Fitting MAMAP(2,2) B+S: B1 > M1, objective = %e\n', robj);
        else
            fprintf('Fitting MAMAP(2,2) B+S: B1 > M1 infeasible\n');
        end
        % pick best solution
        if (~lfeas && ~rfeas) || (eobj < lobj && eobj < robj)
            obj = eobj;
            q1 = p(1);
            q2 = p(1);
            q3 = p(1);
        elseif ~rfeas || lobj < robj
            obj = lobj;
            [q1,q2,q3] = fit(lF1, lS11);
        else
            obj = robj;
            [q1,q2,q3] = fit(rF1, rS11);
        end
        
    end

% used for non-degenerate cases: hybrid space
    function [feas,bobj,bB1,bS11] = solve_nonlinear_hybrid(side)
        vB1 = sdpvar(1,1);
        vS11 = sdpvar(1,1);
        if form == 1
            vq2 = sdpvar(1,1);
            scale = U(11);
            q1exp = -(G(12)*p(1) - vB1*G(3)*p(1) + vq2*(G(11)*G(3) - G(12)*G(2)))/Y(3);
            q2num = U(7)*vB1^2*p(1)^2 + U(8)*vB1*p(1)^2 + U(9)*vS11 + U(10)*p(1)^2;
            q2den = p(1)*(U(11)*vB1 + U(12));
            q3exp = +(G(10)*p(1) - vB1*G(1)*p(1) + vq2*(G(11)*G(1) - G(10)*G(2)))/Y(3);
            hcstr = [0 <= q1exp <= 1, 0 <= vq2 <= 1, vq2*q2den/scale == q2num/scale, 0 <= q3exp <= 1, 0 <= vS11 <= 1, vB1 >= 0];
        else
            vq3 = sdpvar(1,1);
            scale = V(11);
            q1exp = +(E(10)*p(1) - vB1*E(2)*p(1) + vq3*(E(11)*E(2) - E(10)*E(3)))/Z(3);
            q2exp = -(E(9)*p(1) - vB1*E(1)*p(1) + vq3*(E(11)*E(1) - E(3)*E(9)))/Z(3);
            q3num = V(7)*vB1^2*p(1)^2 + V(8)*vB1*p(1)^2 + V(9)*p(1)^2 + V(10)*vS11;
            q3den = p(1)*(V(11)*vB1 + V(12));
            hcstr = [0 <= q1exp <= 1, 0 <= q2exp <= 1, vq3*q3den/scale == q3num/scale, 0 <= vq3 <= 1, 0 <= vS11 <= 1, vB1 >= 0];
        end
        if strcmp(side,'left')
            hcstr = [hcstr, vB1 <= M1-tol];
        else
            hcstr = [hcstr, vB1 >= M1+tol];
        end
        hobj = make_objective(vB1,vS11);
        hsol = solvesdp(hcstr, hobj, yalmip_nonlinear_opt);
        if hsol.problem == 0
            fprintf('Fitting MAMAP(2,2) B+S: solver (hybrid space) objective = %f\n', double(hobj));
            feas = 1;
            bobj = double(hobj);
            bB1 = double(vB1);
            bS11 = double(vS11);
        elseif hsol.problem == 1 || hsol.problem == 12
            fprintf('Fitting MAMAP(2,2) B+S: program is infeasible\n');
            feas = 0;
            bobj = nan;
            bB1 = nan;
            bS11 = nan;
        else
            fname = tempname;
            save(fname,'map','p','B','S');
            error('Fitting MAMAP(2,2) B+S: solver (hybrid space) error: %s, input saved to %s\n', hsol.info, fname);
        end
    end

    function obj = make_objective(vB1, vS11)
        obj = classWeights(1) * bsWeights(1) * (vB1/B(1) - 1)^2 + ...
            classWeights(1) * bsWeights(2) * (vS11/S(1,1) - 1)^2;
    end

    function x = solve_quadprog()
        fprintf('Fitting MAMAP(2,2) B+S: running quadratic programming solver...\n');
        options = optimset('Algorithm','interior-point-convex','Display','none');
        %[x,fx,xflag] = quadprog(H, h, A, b, [], [], [], [], [], options);
        lb = 1e-6*ones( size(A,2),1);
        ub = 1e6*ones( size(A,2),1);
        [x,fx,xflag]=QP(H, h, A, b, Aeq, beq, lb, ub, options);
        
        if xflag ~= 1
            error('Quadratic programming solver failed: %d\n', exit);
        end
        fit_error = fx + length(x);
        fprintf('Fitting MAMAP(2,2) B+S: error = %f\n', fit_error);
    end

end
