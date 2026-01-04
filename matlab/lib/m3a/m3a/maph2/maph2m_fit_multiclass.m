function [maph,fB] = maph2m_fit_multiclass(aph,p,B,classWeights)
% Performs approximate fitting of a MAPH given the underlying APH(2),
% the class probabilities (always fitted exactly), the backward moments,
% and the one-step class transition probabilities. 
% Input
% - aph:  second-order APH underlying the MAPH[m]
% - p:    vector of class probabilities
% - B:    vector of backward moments
% - classWeights: optional vector of weights for each class
% Output
% - maph: fitted MAPH[m]
% - fB:   vector of optimal feasible backward moments

if (size(aph{1},1) ~= 2)
	error('Underlying APH must be of second-order.');
end
if (aph{1}(2,1) ~= 0)
    error('Underlying APH must be acyclic');
end
if (aph{2}(1,2) ~= 0 || aph{2}(2,2) ~= 0)
    error('Underlying APH must be in canonical acyclic form');
end

fprintf('Fitting MAPH(2,m)\n');

% number of classes
k = length(p);

% default weights
if nargin == 3 || isempty(classWeights)
    classWeights = ones(k,1);
end

% result
maph = cell(1,2+k);
maph{1} = aph{1};
maph{2} = aph{2};

% get parameters of the underlying AMAP(2)
h1 = -1/aph{1}(1,1);
h2 = -1/aph{1}(2,2);
r1 = aph{1}(1,2) * h1;

% set tolerance constants
degentol = 1e-6;
feastol  = 1e-8;

q = zeros(2,k);

if abs(1-r1) < degentol
    
    % DEGENERATE
    fprintf('Fitting MAPH(2,m): detected degenerate form\n');
    
    % only one degree of freedom: match class probabilites
    for c = 1:k
        q(1,c) = p(c);
        q(2,c) = p(c);
    end
   
else
    
    % FULL
    
    % coefficients: q(j,c) = F(c) * q_b(j,c) + q_0(j,c)
    q_b = zeros(2,k);
    q_0 = zeros(2,k);
    for c = 1:k
        % q1
        q_b(1,c) = p(c) * ( 1/(h2*(r1 - 1)) );
        q_0(1,c) = p(c) * ( -(h1 + h2)/(h2*(r1 - 1)) );
        % q2
        q_b(2,c) = p(c) * ( 1/(h2*r1) );
        q_0(2,c) = p(c) * ( -h1/(h2*r1) );
    end

    % inequality constraints
    A = zeros(4*k,k);
    b = zeros(4*k,1);
    for c = 1:k
        for j = 1:2
            row = ((c-1)*4+(j-1)*2);
            col = (c-1);
            % q_j <= 1
            A(row+1,col+1) =     q_b(j,c);
            b(row+1)       = 1 - q_0(j,c);
            % q_j >= 0
            A(row+2,col+1) =   - q_b(j,c);
            b(row+2)       =     q_0(j,c);
        end
    end

    % equality constraints
    Aeq = zeros(2,k);
    beq = ones(2,1);
    for c = 1:k
        % equality constraints
        for j = 1:2
            Aeq(j,c) = q_b(j,c);
            beq(j)   = beq(j) - q_0(j,c);
        end
    end

    % objective function
    H = zeros(k, k);
    h = zeros(k, 1);
    for c = 1:k
        base = (c-1);
        backwardWeight  = classWeights(c);
        H(base+1,base+1) =  2/B(c)^2 * backwardWeight;
        h(base+1)        = -2/B(c)   * backwardWeight;
    end

    % solve optimization problem
    fB = solve();

    for c = 1:k
        fprintf('Fitting MAPH(2,m): B(%d) = %f -> %f\n', c, B(c), fB(c));
    end

    % compute parameters of D11,D12,...,D1k
    q = zeros(2,k);
    for c = 1:k
        for j = 1:2
            q(j,c) = fB(c) * q_b(j,c) + q_0(j,c);
        end
    end

end

% check parameter feasibility
for c = 1:k
    if ~(isfeasible(q(1,:)) && isfeasible(q(2,:)))
        error('Fitting MAPH(2,m): Feasibility could not be restored');
    end
end
% parameters feasible within feastol: restrict to [0,1]
q(1,:) = fix(q(1,:));
q(2,:) = fix(q(2,:));

% compute D11,D12,...,D1k
for c = 1:k
    maph{2+c} = maph{2} .* [q(1,c) 0; q(2,c) 0];
end

    function feas = isfeasible(qj)
        feas = min(qj) >= -feastol && sum(qj) <= (1+feastol);
    end

    function QJ = fix(qj)
        QJ = zeros(1,k);
        for cc = 1:k
            QJ(cc) = max(qj(cc),0);
        end
        QJ = QJ ./ sum(QJ);
    end

    function x = solve()
        fprintf('Fitting MAPH(2,m): running quadratic programming solver...\n');
        options = optimset('Algorithm','interior-point-convex','Display','none');
        %[x,fx,xflag] = quadprog(H, h, A, b, Aeq, beq, [], [], [], options);
        lb = 1e-6*ones( size(A,2),1);
        ub = 1e6*ones( size(A,2),1);
        [x,fx,xflag]=QP(H, h, A, b, Aeq, beq, lb, ub, options);        
        if xflag ~= 1
            error('Quadratic programming solver failed: %d\n', exit);
        end
        fit_error = fx + length(x);
        fprintf('Fitting MAPH(2,m): error = %f\n', fit_error);
    end

end % end function