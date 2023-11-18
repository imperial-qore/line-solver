function [mmap,fF,fS,exact] = mamap22_fit_fs_multiclass(map,p,F,S,classWeights,fsWeights,adjust)
% Performs approximate fitting of a MMAP given the underlying AMAP(2),
% the class probabilities (always fitted exactly), the forward moments,
% and the one-step class transition probabilities.
% Input
% - map:  second-order AMAP underlying the MAMAP[2]
% - p:    vector of class probabilities
% - F:    vector of forward moments
% - S:    matrix of the class transition probabilities
% - classWeights: optional vector of weights for each class
% - fsWeights:    optional 2-vector of weights of forward moments and the
%                 class transition probabilities
% Output
% - mmap: fitted MAMAP[2]
% - fF:   vector of optimal feasible forward moments
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

fprintf('Fitting MAMAP(2,2) F+S: form = %d\n', form);

% number of classes
k = length(p);
if k ~= 2
    error('Fitting MAMAP(2,2) F+S: only two classes supported');
end

% default weights to use in the objective function
if nargin < 5 || isempty(classWeights)
    classWeights = ones(k,1);
end
if nargin < 6 || isempty(fsWeights)
    fsWeights = ones(2,1);
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
feastol  = 1e-3;
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

if (form == 1 && (r1 < degentol || r2 > 1-degentol || abs(h1 - h2 + h2*r1) < degentol )) || ...
        (form == 2 && (r2 > 1-degentol || abs(h1 - h2 + h2*r1) < degentol ))
    
    % TODO: transform into a valid second-order Poisson process
    
    % DEGENERATE PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) F+S: detected Poisson process\n');
    
    % return marked poisson process
    h = map_mean(mmap);
    mmap = cell(1,2+k);
    mmap{1} = -1/h;
    mmap{2} =  1/h;
    for c = 1:k
        mmap{2+c} = mmap{2} * p(c);
    end
    
    fF = mmap_forward_moment(mmap,1);
    fS = mmap_sigma(mmap);
    
    return;
    
elseif form == 2 && r2 < degentol && abs(1-r1) < degentol
    
    % DEGENERATE PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) F+S: detected degenerate phase-type form\n');
    
    % only one degree of freedom: match class probabilites
    q1 = p(1);
    q2 = p(1);
    q3 = p(1);
    
elseif form == 1 && r2 < degentol
    
    % CANONICAL PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) F+S: detected canonical phase-type form\n');
    
    % convert to phase-type
    aph = map;
    aph{2}(2,2) = 0;
    aph = map_normalize(aph);
    
    % the forward moments are always equal to the ordinary moments
    % since the backward moments are not provided as arguments,
    % we set the backward moments to the ordinary moments as well
    % TODO: if hypoexponential convert to non-canonical phase-type
    B = zeros(2,1);
    B(1) = map_mean(aph);
    B(2) = map_mean(aph);
    warning('Fitting MAMAP(2,2) F+S: setting backward moments to ordinary moments, you should try to fit B+S');
    
    mmap = maph2m_fit_multiclass(aph, p, B, classWeights);
    
    fF = mmap_forward_moment(mmap, 1);
    fS = mmap_sigma(mmap);
    
    return;
    
elseif (form == 1 && abs(1-r1) < degentol) || ...
        (form == 2 && abs(1-r1) < degentol)
    
    % NON-CANONICAL PHASE_TYPE
    fprintf('Fitting MAMAP(2,2) F+S: detected non-canonical phase-type form\n');
    
    % two degrees of freedom: fit p and F, not S
    [q1,q2,q3] = fit_can1_degen(F(1));
    
    % if unfeasible, perform approximate fitting
    if ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
        % four inequalities, one variable
        A = zeros(4,1);
        b = zeros(4,1);
        % coefficients
        q2_F = - p(1) / (h1 * (r2-1));
        q2_0 = + p(1) * h2 / (h1 * (r2-1));
        q3_F = - p(1) / (h1 * r2);
        q3_0 = + p(1) *  (h1 + h2) / (h1 * r2);
        % q2 >= 0
        A(1,1) = -q2_F;
        b(1)   =  q2_0;
        % q2 <= 1
        A(2,1) = +q2_F;
        b(2)   =  1 - q2_0;
        % q3 >= 0
        A(3,1) = -q3_F;
        b(3)   =  q3_0;
        % q3 <= 1
        A(4,1) = +q3_F;
        b(4)   =  1 - q3_0;
        % objective function
        H =  2 / F(1)^2;
        h = -2 / F(1);
        % solve
        fF1 = solve_quadprog();
        % compute coefficient
        [q1, q2, q3] = fit_can1_degen(fF1);
    end
    
elseif form == 2 && r2 < degentol
    
    % DEGENERATE CASE FOR gamma < 0
    fprintf('Fitting MAMAP(2,2) F+S: detected degenerate MMAP form\n');
    
    if fsWeights(1) > fsWeights(2)
        fprintf('Fitting MAMAP(2,2) F+S: fitting F\n');
        [q1,q2,q3] = fit_can2_degen_forward(F(1));
        if ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
            % four inequalities, one variable
            A = zeros(4,1);
            b = zeros(4,1);
            % coefficients
            q1_F = - p(1) * (r1 - 2) / ((r1-1)*(h1-h2+h2*r1));
            q1_0 = + p(1) * (r1 - 2)*(h1 + h2*r1) / ((r1-1)*(h1-h2+h2*r1));
            q2_F = - p(1) * (r1 - 2) / ((h1-h2+h2*r1));
            q2_0 = + p(1) * (r1 - 2)*h2 / ((h1-h2+h2*r1));
            % q1 >= 0
            A(1,1) = -q1_F;
            b(1)   =  q1_0;
            % q1 <= 1
            A(2,1) = +q1_F;
            b(2)   =  1 - q1_0;
            % q2 >= 0
            A(3,1) = -q2_F;
            b(3)   =  q2_0;
            % q2 <= 1
            A(4,1) = +q2_F;
            b(4)   =  1 - q2_0;
            % objective function
            H =  2 / F(1)^2;
            h = -2 / F(1);
            % solve
            fF1 = solve_quadprog();
            % compute coefficient
            [q1,q2,q3] = fit_can2_degen_forward(fF1);
        end
    else
        fprintf('Fitting MAMAP(2,2) F+S: fitting S\n');
        [q1,q2,q3] = fit_can2_degen_transition(S(1,1));
        if ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
            % nonlinear constraints: use yalmip
            varS11 = sdpvar(1,1);
            options = sdpsettings('verbose',0);
            boundS11 = [varS11 >= 0, varS11 <= p(1)^2];
            q1lb = varS11 >= p(1)^2 * (1 - (1-r1)^2);
            q2ub = varS11 >= p(1)^2 - (1 - p(1))^2;
            constraints = [boundS11, q1lb, q2ub];
            obj = (varS11 - S(1,1))^2;
            sol = solvesdp(constraints, obj, options);
            if sol.problem ~= 0
                error('Fitting MAMAP(2,2) F+S: YALMIP = %s\n', sol.info);
            else
                fprintf('Fitting MAMAP(2,2) F+S: YALMIP objective = %f\n', double(obj));
            end
            % compute coefficient
            [q1,q2,q3] = fit_can2_degen_transition(double(varS11));
        end
    end
    
else
    
    % FULL FORM or "GOOD" poisson process
    
    if (form == 1 && abs(h2-h1*r2) < degentol) || (form == 2 && abs(h1 - h2 - h1*r1 + h1*r1*r2) < degentol)
        fprintf('Fitting MAMAP(2,2) F+S: detected fittable Poisson process\n')
    end
    
    [q1,q2,q3] = fit(F(1), S(1,1));
    
    % options used to solve the nonlinear problem
    use_linprog_prepass = 1;
    tol = 1e-6;
    optim_space = 'hybrid';
    %optim_space = 'hybrid-z';
    %optim_space = 'param';
    %optim_space = 'char';
    
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
        fprintf('Fitting MAMAP(2,2) F+S: exact fit found\n');
        exact = 1;
    elseif adjust
        [q1,q2,q3] = solve_nonlinear(optim_space);
    end
    
end % end if/then/else for degenerate forms

% check parameter feasibility
if adjust && ~(isfeasible(q1) && isfeasible(q2) && isfeasible(q3))
    error('Fitting MAMAP(2,2) F+S: Feasibility could not be restored: q1 = %e, q2 = %e, q3 = %e', q1, q2, q3);
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

fF = mmap_forward_moment(mmap, 1);
fS = mmap_sigma(mmap);

    function feas = isfeasible(q)
        feas = q >= -feastol && q <= (1+feastol);
    end

    function qfix = fix(q)
        qfix = max(min(q,1),0);
    end

    function [q1,q2,q3] = fit_can1(vF1,vS11)
        denum = p(1) * (U(5) * vF1 + U(6));
        if abs(denum) < denumtol
            q1 = p(1);
            q2 = p(1);
            q3 = p(1);
        else
            q2 = (U(1) * vF1^2 * p(1)^2 + U(2) * vF1 * p(1)^2 + U(3) * vS11 + U(4) * p(1)^2)/denum;
            q1 = -(G(15)*p(1) - vF1*G(3)*p(1) + (G(3)*G(14) - G(2)*G(15)) * q2)/Y(2);
            q3 = +(G(13)*p(1) - vF1*G(1)*p(1) + (G(1)*G(14) - G(2)*G(13)) * q2)/Y(2);
        end
    end

    function [q1,q2,q3] = fit_can1_degen(vF1)
        q1 = p(1); % meangingless if really degenerate
        q2 = p(1) * (h2 - vF1) / (h1 * (r2-1));
        q3 = p(1) * (h1 + h2 - vF1) / (h1 * r2);
    end

    function [q1,q2,q3] = fit_can2_degen_forward(vF1)
        q1 = + p(1) * (r1 - 2)*(h1 + h2*r1 - vF1)/((r1-1)*(h1-h2+h2*r1));
        q2 = - p(1) * (vF1 - h2)*(r1 - 2)/((h1-h2+h2*r1));
        q3 = p(1); % meangingless if really degenerate
    end

    function [q1,q2,q3] = fit_can2_degen_transition(vS11)
        q1 = p(1) + (p(1)^2 - vS11)^(1/2)/(r1 - 1);
        q2 = p(1) + (p(1)^2 - vS11)^(1/2);
        q3 = p(1); % meangingless if really degenerate
    end

    function [q1,q2,q3] = fit_can2(vF1,vS11)
        denum = (V(5) * vF1 * p(1) + V(6) * p(1));
        if abs(denum) < denumtol
            q1 = p(1);
            q2 = p(1);
            q3 = p(1);
        else
            q3 = (V(1) * vF1^2 * p(1)^2 + V(2) * vF1 * p(1)^2 + V(3) * p(1)^2 + V(4) * vS11)/denum;
            q1 = -(E(13)*p(1) - vF1*E(2)*p(1) + (E(2)*E(14) - E(3)*E(13)) * q3)/Z(2);
            q2 = +(E(12)*p(1) - vF1*E(1)*p(1) + (E(1)*E(14) - E(3)*E(12)) * q3)/Z(2);
        end
    end

    function [pexp,pFexp,Sexp] = get_characteristics(q1,q2,q3)
        if form == 1
            pexp = G(1)*q1+G(2)*q2+G(3)*q3;
            pFexp = G(13)*q1+G(14)*q2+G(15)*q3;
            Sexp = G(4)*q1^2 + G(5)*q1*q2 + G(6)*q1*q3 + G(7)*q2^2 + G(8)*q2*q3 + G(9)*q3^2;
        else
            pexp = E(1)*q1+E(2)*q2+E(3)*q3;
            pFexp = (E(12)*q1+E(13)*q2+E(14)*q3);
            Sexp = E(4)*q1*q2 + E(5)*q1*q3 + E(6)*q2^2 + E(7)*q2*q3 + E(8)*q3^2;
        end
    end

    function opt = get_yalmip_linear_options()
        % we do not set solver so that YALMIP picks the best
        opt = sdpsettings('verbose',0);
    end

    function [q1,q2,q3,obj] = solve_nonlinear(optim_space)
        
        fprintf('Fitting MAMAP(2,2) F+S: running approximate fitting (nonlinear)\n');
        
        if strcmp(optim_space, 'param')
            [q1,q2,q3,obj] = solve_nonlinear_param();
        else
            if strcmp(optim_space, 'hybrid')
                solve_nonlinear_func = @solve_nonlinear_hybrid;
            elseif strcmp(optim_space, 'hybrid-z')
                solve_nonlinear_func = @solve_nonlinear_hybrid_z;
            elseif strcmp(optim_space, 'char')
                solve_nonlinear_func = @solve_nonlinear_char;
            else
                error('Fitting MAMAP(2,2) F+S: invalid value for option "optim_space"');
            end
            % trivial solution
            eobj = make_objective(M1, p(1)^2);
            fprintf('Fitting MAMAP(2,2) F+S: F1 == M1, objective = %e\n', eobj);
            % left solution
            [lfeas,lobj,lF1,lS11] = solve_nonlinear_func('left');
            if lfeas
                fprintf('Fitting MAMAP(2,2) F+S: F1 < M1, objective = %e\n', lobj);
            else
                fprintf('Fitting MAMAP(2,2) F+S: F1 < M1 infeasible\n');
            end
            % right solution
            [rfeas,robj,rF1,rS11] = solve_nonlinear_func('right');
            if rfeas
                fprintf('Fitting MAMAP(2,2) F+S: F1 > M1, objective = %e\n', robj);
            else
                fprintf('Fitting MAMAP(2,2) F+S: F1 > M1 infeasible\n');
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
        
    end

% used for non-degenerate cases: parmeter space
    function [q1,q2,q3,obj] = solve_nonlinear_param()
        % find initial feasible solution
        vq = sdpvar(3,1);
        z = sdpvar(2,1);
        % define expressions
        [pexp,pFexp,Sexp] = get_characteristics(vq(1),vq(2),vq(3));
        % set constraints
        cstr_exact_p = pexp == p(1);
        cstr_tighten_F = -z(1) <= (pFexp - p(1)*F(1))/(p(1)*F(1)) <= z(1);
        cstr_tighten_S = -z(2) <= (Sexp - S(1,1))/S(1,1) <= z(2);
        cstr = [0 <= vq(1) <= 1, ...
            0 <= vq(2) <= 1, ...
            0 <= vq(3) <= 1, ...
            cstr_exact_p, ...
            cstr_tighten_F,...
            cstr_tighten_S,...
            z(1) >= 0, z(2) >= 0];
        % set objective
        zobj = fsWeights(1) * z(1)^2 + fsWeights(2) * z(2)^2;
        % run solver
        zsol = solvesdp(cstr, zobj, yalmip_nonlinear_opt);
        % check output
        if zsol.problem == 0
            q1 = double(vq(1));
            q2 = double(vq(2));
            q3 = double(vq(3));
            obj = double(zobj);
            fprintf('Fitting MAMAP(2,2) F+S: solver (parameter space) objective = %e\n', obj);
        else
            fname = tempname;
            save(fname,'map','p','F','S');
            error('Fitting MAMAP(2,2) F+S: solver (parameter space) error: %s, input saved to %s\n', zsol.info, fname);
        end
    end

% used for non-degenerate cases: hybrid space
    function [feas,bobj,bF1,bS11] = solve_nonlinear_hybrid(side)
        vF1 = sdpvar(1,1);
        vS11 = sdpvar(1,1);
        if form == 1
            vq2 = sdpvar(1,1);
            scale = U(5);
            q1exp = -(G(15)*p(1) - vF1*G(3)*p(1) + vq2*(G(14)*G(3)-G(15)*G(2)))/Y(2);
            q2num = U(1)*vF1^2*p(1)^2 + U(2)*vF1*p(1)^2 + U(3)*vS11 + U(4)*p(1)^2;
            q2den = p(1)*(U(5)*vF1 + U(6));
            q3exp = +(G(13)*p(1) - vF1*G(1)*p(1) + vq2*(G(14)*G(1)-G(13)*G(2)))/Y(2);
            hcstr = [0 <= q1exp <= 1, 0 <= vq2 <= 1, vq2*q2den/scale == q2num/scale, 0 <= q3exp <= 1, 0 <= vS11 <= 1, vF1 >= 0];
        else
            vq3 = sdpvar(1,1);
            scale = V(5);
            q1exp = -(E(13)*p(1) - vF1*E(2)*p(1) + vq3*(E(14)*E(2)-E(13)*E(3)))/Z(2);
            q2exp = +(E(12)*p(1) - vF1*E(1)*p(1) + vq3*(E(14)*E(1)-E(12)*E(3)))/Z(2);
            q3num = V(1)*vF1^2*p(1)^2 + V(2)*vF1*p(1)^2 + V(3)*p(1)^2 + V(4)*vS11;
            q3den = p(1)*(V(5)*vF1 + V(6));
            hcstr = [0 <= q1exp <= 1, 0 <= q2exp <= 1, vq3*q3den/scale == q3num/scale, 0 <= vq3 <= 1, 0 <= vS11 <= 1, vF1 >= 0];
        end
        if strcmp(side,'left')
            hcstr = [hcstr, vF1 <= M1-tol];
        else
            hcstr = [hcstr, vF1 >= M1+tol];
        end
        hobj = make_objective(vF1,vS11);
        hsol = solvesdp(hcstr, hobj, yalmip_nonlinear_opt);
        if hsol.problem == 0
            fprintf('Fitting MAMAP(2,2) F+S: solver (hybrid space) objective = %f\n', double(hobj));
            feas = 1;
            bobj = double(hobj);
            bF1 = double(vF1);
            bS11 = double(vS11);
        elseif hsol.problem == 1 || hsol.problem == 12
            fprintf('Fitting MAMAP(2,2) F+S: program is infeasible\n');
            feas = 0;
            bobj = nan;
            bF1 = nan;
            bS11 = nan;
        else
            fname = tempname;
            save(fname,'map','p','F','S');
            error('Fitting MAMAP(2,2) F+S: solver (hybrid space) error: %s, input saved to %s\n', hsol.info, fname);
        end
    end

% used for non-degenerate cases: hybrid space with z variable and implicit vS11
    function [feas,bobj,bF1,bS11] = solve_nonlinear_hybrid_z(side)
        vF1 = sdpvar(1,1);
        z = sdpvar(1,1);
        if form == 1
            vq2 = sdpvar(1,1);
            q1exp = -(G(15)*p(1) - vF1*G(3)*p(1) + vq2*(G(14)*G(3)-G(15)*G(2)))/Y(2);
            q2den = p(1)*(U(5)*vF1 + U(6));
            q3exp = +(G(13)*p(1) - vF1*G(1)*p(1) + vq2*(G(14)*G(1)-G(13)*G(2)))/Y(2);
            S11exp = (vq2*q2den - p(1)^2*(U(1)*vF1^2 + U(2)*vF1 + U(4)))/U(3);
            tighten = -z <= S11exp - S(1,1) <= +z;
            hcstr = [0 <= q1exp <= 1, 0 <= vq2 <= 1, 0 <= q3exp <= 1, vF1 >= 0, tighten];
        else
            vq3 = sdpvar(1,1);
            q1exp = -(E(13)*p(1) - vF1*E(2)*p(1) + vq3*(E(14)*E(2)-E(13)*E(3)))/Z(2);
            q2exp = +(E(12)*p(1) - vF1*E(1)*p(1) + vq3*(E(14)*E(1)-E(12)*E(3)))/Z(2);
            q3den = p(1)*(V(5)*vF1 + V(6));
            S11exp = (vq3*q3den - p(1)^2*(V(1)*vF1^2 + V(2)*vF1 + V(3)))/V(4);
            tighten = -z <= S11exp - S(1,1) <= +z;
            hcstr = [0 <= q1exp <= 1, 0 <= q2exp <= 1, 0 <= vq3 <= 1, vF1 >= 0, tighten];
        end
        if strcmp(side,'left')
            hcstr = [hcstr, vF1 <= M1-tol];
        else
            hcstr = [hcstr, vF1 >= M1+tol];
        end
        hobj = fsWeights(1)*(vF1/F(1)-1)^2 + fsWeights(2)*(z/S(1,1))^2;
        hsol = solvesdp(hcstr, hobj, yalmip_nonlinear_opt);
        if hsol.problem == 0
            fprintf('Fitting MAMAP(2,2) F+S: solver (hybrid-z space) objective = %f\n', double(hobj));
            feas = 1;
            bobj = double(hobj);
            bF1 = double(vF1);
            bS11 = double(S11exp);
        elseif hsol.problem == 1 || hsol.problem == 12
            fprintf('Fitting MAMAP(2,2) F+S: program is infeasible\n');
            feas = 0;
            bobj = nan;
            bF1 = nan;
            bS11 = nan;
        else
            fname = tempname;
            save(fname,'map','p','F','S');
            error('Fitting MAMAP(2,2) F+S: solver (hybrid-z space) error: %s, input saved to %s\n', hsol.info, fname);
        end
    end

% used for non-degenerate cases: characteristic space
    function [feas,bobj,bF1,bS11] = solve_nonlinear_char(side)
        is_left = 0;
        if strcmp(side,'left')
            fprintf('Fitting MAMAP(2,2) F+S: optimizing for F11 < M1\n');
            is_left = 1;
        else
            fprintf('Fitting MAMAP(2,2) F+S: optimizing for F11 > M1\n');
        end
        % find initial feasible solution
        if use_linprog_prepass
            % parameter variables vq and auxiliary variables
            vq = sdpvar(3,1);
            z = sdpvar(1,1);
            % define expressions
            [pexp,pFexp,Sexp] = get_characteristics(vq(1),vq(2),vq(3));
            % set constraints
            if is_left
                cstr_side = pFexp <= pexp*(M1-tol);
            else
                cstr_side = pFexp >= pexp*(M1+tol);
            end
            cstr_exact_p = pexp == p(1);
            cstr_tighten = -z <= pFexp-p(1)*F(1) <=z;
            fcstr = [0 <= vq(1) <= 1, ...
                0 <= vq(2) <= 1, ...
                0 <= vq(3) <= 1, ...
                cstr_exact_p, ...
                cstr_side,...
                cstr_tighten,...
                z >= 0];
            % run linear programming solver
            lsol = solvesdp(fcstr,z, get_yalmip_linear_options());
            % check output
            if lsol.problem == 0
                feas = 1;
                iF1 = double(pFexp / pexp);
                iS11 = double(Sexp);
            elseif lsol.problem == 1 || lsol.problem == 12
                fprintf('Fitting MAMAP(2,2) F+S: linear program is infeasible\n');
                feas = 0;
            else
                fname = tempname;
                save(fname,'map','p','F','S');
                error('Fitting MAMAP(2,2) F+S: linear program error: %s, input saved to %s\n', lsol.info, fname);
            end
        else
            % assume feasible
            feas = 1;
        end
        % if feasible, solve nonlinear program
        if feas
            % declare variables for optimization in characteristic space
            vF1 = sdpvar(1,1);
            vS11 = sdpvar(1,1);
            % initialize variables
            if use_linprog_prepass
                assign(vF1,iF1);
                assign(vS11,iS11);
            end
            % define objective function
            cobj = make_objective(vF1, vS11);
            % define constraints
            if is_left
                if form == 1
                    [cstr,~] = make_constraints_can1(p(1), vF1, vS11);
                else
                    [cstr,~] = make_constraints_can2(p(1), vF1, vS11);
                end
            else
                if form == 1
                    [~,cstr] = make_constraints_can1(p(1), vF1, vS11);
                else
                    [~,cstr] = make_constraints_can2(p(1), vF1, vS11);
                end
            end
            % solve nonlinear problem
            csol = solvesdp(cstr, cobj, yalmip_nonlinear_opt);
            % check result
            if csol.problem ~= 0
                fname = tempname;
                save(fname,'map','p','F','S');
                error('Fitting MAMAP(2,2) F+S: solver (characteristic space) error: %s, input saved to %s\n', csol.info, fname);
            else
                bobj = double(cobj);
                bF1 = double(vF1);
                bS11 = double(vS11);
                fprintf('Fitting MAMAP(2,2) F+S: solver (characteristic space) objective = %e\n', double(cobj));
            end
        else
            % infeasible subproblem
            bobj = nan;
            bF1 = nan;
            bS11 = nan;
        end
    end

    function obj = make_objective(vF1, vS11)
        obj = classWeights(1) * fsWeights(1) * (vF1/F(1) - 1)^2 + ...
            classWeights(1) * fsWeights(2) * (vS11/S(1,1) - 1)^2;
    end

    function [q1n,q1d,q2n,q2d,q3] = make_q_can1(p, vF1, vS11)
        % define coefficients of polynomials
        q1n_F = p^2*(r1 - 1)*(r1*r2 - r2 + 1);
        q1n_s = -r1*(h1 - h2 + h2*r1);
        q1n_1 = h1*p^2*(r1*r2 - r2 + 1);
        q1d_F = p*(r1 - 1)*(r1*r2 - r2 + 1);
        q1d_1 = -p*(r1 - 1)*(h1 - h1*r2 + h2*r1);
        q2n_Fsq = -p^2*(r1*r2 - r2 + 1)^2;
        q2n_F = p^2*(r1*r2 - r2 + 1)*(2*h1 - h1*r1 - 2*h1*r2 + 3*h2*r1 - h2*r1^2 + h2*r1^2*r2 + h1*r1*r2 - h2*r1*r2);
        q2n_s = -r1*(r2 - 1)*(h1 - h2 + h2*r1)^2;
        q2n_1 = -p^2*(r1*r2 - r2 + 1)*(h2^2*r1 - h1^2*r2 + h1^2 + h1*h2*r1 - h1*h2*r1*r2);
        q2d_F = p*r1*(r2 - 1)*(r1*r2 - r2 + 1)*(h1 - h2 + h2*r1);
        q2d_1 = -p*r1*(r2 - 1)*(h1 - h1*r2 + h2*r1)*(h1 - h2 + h2*r1);
        q3_F = -(p*(r1*r2 - r2 + 1))/(r1*r2*(h1 + h2*(r1 - 1)));
        q3_1 = (p*(h1 + h2*r1)*(r1*r2 - r2 + 1))/(r1*r2*(h1 - h2 + h2*r1));
        % return expressions
        q1n = q1n_F * vF1 + q1n_s * vS11 + q1n_1;
        q1d = q1d_F * vF1 + q1d_1;
        q2n = q2n_Fsq * vF1^2 + q2n_F * vF1 + q2n_s * vS11 + q2n_1;
        q2d = q2d_F * vF1 + q2d_1;
        q3  = q3_F * vF1 + q3_1;
    end

    function [lcstr,rcstr] = make_constraints_can1(p, vF1, vS11)
        
        % get fitting expressions
        [q1n,q1d,q2n,q2d,q3] = make_q_can1(p, vF1, vS11);
        
        % these are in common for the two subproblems
        q3lb = q3 >= 0;
        q3ub = q3 <= 1;
        
        scale = 1/abs(U(5));
        %scale = 1/abs(q2d{M1-tol});
        
        % left constraints: F(1) < M1
        q1lb = q1n >= 0;
        q1ub = q1n <= q1d;
        if h2 < h1/(1-r1)
            % do not change direction of inequalities for q2
            q2lb = scale*(q2n) >= 0;
            q2ub = scale*(q2n) <= scale*(q2d);
        else
            % change direction of inequalities for q2
            q2lb = scale*(q2n) <= 0;
            q2ub = scale*(q2n) >= scale*(q2d);
        end
        F1bound = vF1 <= M1-tol;
        lcstr = [q1lb,q1ub,q2lb,q2ub,q3lb,q3ub,F1bound,vS11 >= 0, vS11 <= 1];
        
        % right constraints: F(1) > M1
        q1lb = q1n <= 0;
        q1ub = q1n >= q1d;
        if h2 > h1/(1-r1)
            % do not change direction of inequalities for q2
            q2lb = scale*(q2n) >= 0;
            q2ub = scale*(q2n) <= scale*(q2d);
        else
            % change direction of inequalities for q2
            q2lb = scale*(q2n) <= 0;
            q2ub = scale*(q2n) >= scale*(q2d);
        end
        F1bound = vF1 >= M1+tol;
        rcstr = [q1lb,q1ub,q2lb,q2ub,q3lb,q3ub,F1bound,vS11 >= 0, vS11 <= 1];
    end

    function [q1n,q1d,q2,q3n,q3d] = make_q_can2(p, vF1, vS11)
        % define coefficients of polynomials
        q1n_F = p^2*(r1 - 1)*(r1 + r2 - r1*r2 - 2);
        q1n_s = h1 - h2 + h2*r1;
        q1n_1 = h1*p^2*(r1 + r2 - r1*r2 - 2);
        q1d_F = p*(r1 - 1)*(r1 + r2 - r1*r2 - 2);
        q1d_1 = p*(r1 - 1)*(h1 + h2 - h1*r2);
        q2_F = p*(r1 + r2 - r1*r2 - 2)/(h2 - h1 + h1*r2 - h2*r1 - h2*r2 + h2*r1*r2);
        q2_1 = -h2*p*(r1 + r2 - r1*r2 - 2)/(h2 - h1 + h1*r2 - h2*r1 - h2*r2 + h2*r1*r2);
        q3n_Fsq = p^2*(r1 + r2 - r1*r2 - 2)^2;
        q3n_F = p^2*(r1 + r2 - r1*r2 - 2)*(2*h1 + 2*h2 - h1*r2 - h2*r2 + h2*r1*r2);
        q3n_s = -(r2 - 1)*(h1 - h2 + h2*r1)^2;
        q3n_1 = -h2*p^2*(2*h1 - h1*r2 + h2*r1)*(r1 + r2 - r1*r2 - 2);
        q3d_F = p*r2*(h1 - h2 + h2*r1)*(r1 + r2 - r1*r2 - 2);
        q3d_1 = p*r2*(h1 + h2 - h1*r2)*(h1 - h2 + h2*r1);
        % return expressions
        q1n = q1n_F * vF1 + q1n_s * vS11 + q1n_1;
        q1d = q1d_F * vF1 + q1d_1;
        q2  = q2_F * vF1 + q2_1;
        q3n = q3n_Fsq * vF1^2 + q3n_F * vF1 + q3n_s * vS11 + q3n_1;
        q3d = q3d_F * vF1 + q3d_1;
    end

    function [lcstr,rcstr] = make_constraints_can2(p, vF1, vS11)
        
        % get fitting expressions
        [q1n,q1d,q2,q3n,q3d] = make_q_can2(p, vF1, vS11);
        
        % these are in common for the two subproblems
        q2lb = q2 >= 0;
        q2ub = q2 <= 1;
        
        scale = 1/abs(V(5));
        %scale = 1/abs(q3d{M1-tol});
        
        % left constraints: F1 < M1
        q1lb = q1n <= 0;
        q1ub = q1n >= q1d;
        if (h1 - h2 + h2*r1) < 0
            q3lb = scale*(q3n) <= 0;
            q3ub = scale*(q3n) >= scale*(q3d);
        else
            q3lb = scale*(q3n) >= 0;
            q3ub = scale*(q3n) <= scale*(q3d);
        end
        F1bound = vF1 <= M1-tol;
        lcstr = [q1lb,q1ub,q2lb,q2ub,q3lb,q3ub,F1bound,vS11 >= 0, vS11 <= 1];
        
        % rigt costraints: F1 > M1
        q1lb = q1n >= 0;
        q1ub = q1n <= q1d;
        if (h1 - h2 + h2*r1) < 0
            q3lb = scale*(q3n) >= 0;
            q3ub = scale*(q3n) <= scale*(q3d);
        else
            q3lb = scale*(q3n) <= 0;
            q3ub = scale*(q3n) >= scale*(q3d);
        end
        F1bound = vF1 >= M1+tol;
        rcstr = [q1lb,q1ub,q2lb,q2ub,q3lb,q3ub,F1bound,vS11 >= 0, vS11 <= 1];
    end

% used for degenerate case
    function x = solve_quadprog()
        fprintf('Fitting MAMAP(2,2) F+S: running quadratic programming solver...\n');
        options = optimset('Algorithm','interior-point-convex','Display','none');
        %[x,fx,xflag] = quadprog(H, h, A, b, [], [], [], [], [], options);
        lb = 1e-6*ones( size(A,2),1);
        ub = 1e6*ones( size(A,2),1);
        [x,fx,xflag]=QP(H, h, A, b, Aeq, beq, lb, ub, options);
        
        if xflag ~= 1
            error('Quadratic programming solver failed: %d\n', exit);
        end
        fit_error = fx + length(x);
        fprintf('Fitting MAMAP(2,2) F+S: error = %f\n', fit_error);
    end

end