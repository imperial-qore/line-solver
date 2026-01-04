function [mmap,fF,fB] = mamap2m_fit_fb_multiclass(map,p,F,B,classWeights,fbWeights)
% Performs approximate fitting of a MMAP given the underlying MAP,
% the class probabilities (always fitted exactly), the forward moments,
% and the backward moments.
% Input
% - map:  second-order AMAP underlying the MAMAP[m]
% - p:    vector of class probabilities
% - F:    vector of forward moments
% - B:    vector of backward moments
% - classWeights: optional vector of weights for each class
% - fbWeights: optional 2-vector of weights of forward and backward moments
% Output
% - mmap: fitted MAMAP[m]
% - fF:   vector of optimal feasible forward moments
% - fB:   vector of optimal feasible backward moments

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

%fprintf('Fitting MAMAP(2,m) F+B: form = %d\n', form);

% number of classes
k = length(p);

% default weights to use in the objective function
if nargin < 5 || isempty(classWeights)
	classWeights = ones(k,1);
end
if nargin < 6
    fbWeights = ones(2,1);
end

% result
mmap = cell(1,2+k);
mmap{1} = map{1};
mmap{2} = map{2};

h1 = -1/map{1}(1,1);
h2 = -1/map{1}(2,2);
r1 = map{1}(1,2) * h1;
r2 = map{2}(2,2) * h2;

degentol = 1e-8;

if (form == 1 && (r1 < degentol || r2 > 1-degentol || abs(h2-h1*r2) < degentol || abs(h1 - h2 + h2*r1) < degentol )) || ...
   (form == 2 && (r2 > 1-degentol || abs(h1 - h2 + h2*r1) < degentol || abs(h1 - h2 - h1*r1 + h1*r1*r2) < degentol ))
    
    % POISSON PROCESS
 %   fprintf('Fitting MAMAP(2,m) F+B: detected Poisson process\n');
    
    % return marked poisson process
    h = map_mean(mmap);
    mmap = cell(1,2+k);
    mmap{1} = -1/h;
    mmap{2} =  1/h;
    for c = 1:k
        mmap{2+c} = mmap{2} * p(c);
    end
    
    fF = mmap_forward_moment(mmap,1);
    fB = mmap_backward_moment(mmap,1);
    
    return;
    
elseif (form == 2 && r2 < degentol && abs(1-r1) < degentol)
    
    % DEGENERATE PHASE_TYPE
  %  fprintf('Fitting MAMAP(2,m) F+B: detected degenerate phase-type form\n');
    
    % compute parameters of D11,D12,...,D1k
    q = zeros(3,k);
    for c = 1:k
        q(1,c) = p(c);
        q(2,c) = p(c);
        q(3,c) = p(c);
    end
    
elseif form == 1 && r2 < degentol
    
    % CANONICAL PHASE_TYPE
%    fprintf('Fitting MAMAP(2,m) F+B: detected canonical phase-type form\n');
    
%    fprintf('Fitting MAMAP(2,m) F+B: fitting backward\n');
    
    % convert to phase-type
    aph = map;
    aph{2}(2,2) = 0;
    aph = map_normalize(aph);
    
    mmap = maph2m_fit_multiclass(aph, p, B, classWeights);
    
    fF = mmap_forward_moment(mmap, 1);
    fB = mmap_backward_moment(mmap, 1);

    return;

elseif (form == 1 && abs(1-r1) < degentol) || ...
       (form == 2 && abs(1-r1) < degentol)
    
    % NON-CANONICAL PHASE_TYPE
 %   fprintf('Fitting MAMAP(2,m) F+B: detected non-canonical phase-type form\n');
    
    % coefficients: q(j,c) = F(c) * q_f(j,c) + q_0(j,c)
    q_f = zeros(2,k);
    q_0 = zeros(2,k);
    for c = 1:k
        % q2
        q_f(1,c) = p(c) * ( -1/((h1 + h2*(r1 - 1))*(r2 - 1)*(r1 + r2 - r1*r2)) );
        q_0(1,c) = p(c) * ( h2/((r2 - 1)*(r1 + r2 - r1*r2)*(h1 - h2 + h2*r1)) );
        % q3
        q_f(2,c) = p(c) * ( -1/(r2*(h1 + h2*(r1 - 1))*(r1 + r2 - r1*r2)) );
        q_0(2,c) = p(c) * ( (h1 + h2*r1)/(r2*(r1 + r2 - r1*r2)*(h1 - h2 + h2*r1)) );
    end

    % inequality constraints
    A = zeros(4*k,k);
    b = zeros(4*k,1);
    for c = 1:k
        for j = 1:2
            row = ((c-1)*4+(j-1)*2);
            col = (c-1);
            % q_j <= 1
            A(row+1,col+1) =     q_f(j,c);
            b(row+1)       = 1 - q_0(j,c);
            % q_j >= 0
            A(row+2,col+1) =   - q_f(j,c);
            b(row+2)       =     q_0(j,c);
        end
    end

    % equality constraints
    Aeq = zeros(2,k);
    beq = ones(2,1);
    for c = 1:k
        % equality constraints
        for j = 1:2
            Aeq(j,c) = q_f(j,c);
            beq(j)   = beq(j) - q_0(j,c);
        end
    end

    % objective function
    H = zeros(k, k);
    h = zeros(k, 1);
    for c = 1:k
        base = (c-1);
        forwardWeight  = (classWeights(c) * fbWeights(1));
        H(base+1,base+1) =  2/F(c)^2 * forwardWeight;
        h(base+1)        = -2/F(c)   * forwardWeight;
    end

    % solve optimization problem
    fF = solve();

    for c = 1:k
  %      fprintf('Fitting MAMAP(2,m) F+B: F(%d) = %f -> %f\n', c, F(c), fF(c));
    end

    % compute parameters of D11,D12,...,D1k
    q = zeros(3,k);
    for c = 1:k
        q(1,c) = 1/k;
    end
    for c = 1:k
        q(2,c) = fF(c) * q_f(1,c) + q_0(1,c);
        q(3,c) = fF(c) * q_f(2,c) + q_0(2,c);
    end
    
elseif form == 2 && r2 < degentol
    
    % DEGENERATE CASE FOR gamma < 0
   % fprintf('Fitting MAMAP(2,m) F+B: detected degenerate MMAP form\n');
    
    if fbWeights(1) >= fbWeights(2)
        
    %    fprintf('Fitting MAMAP(2,m) F+B: fitting forward\n');

        % coefficients: q(j,c) = F(c) * q_f(j,c) + q_0(j,c)
        q_f = zeros(2,k);
        q_0 = zeros(2,k);
        for c = 1:k
            % q1
            q_f(1,c) = p(c) * ( -(r1 - 2)/((h1 + h2*(r1 - 1))*(r1 - 1)) );
            q_0(1,c) = p(c) * ( 1 - (h1 + h2)/((r1 - 1)*(h1 - h2 + h2*r1)) );
            % q2
            q_f(2,c) = p(c) * ( -(r1 - 2)/(h1 + h2*(r1 - 1)) );
            q_0(2,c) = p(c) * ( (h2*(r1 - 2))/(h1 - h2 + h2*r1) );
        end

        % inequality constraints
        A = zeros(4*k,k);
        b = zeros(4*k,1);
        for c = 1:k
            for j = 1:2
                row = ((c-1)*4+(j-1)*2);
                col = (c-1);
                % q_j <= 1
                A(row+1,col+1) =     q_f(j,c);
                b(row+1)       = 1 - q_0(j,c);
                % q_j >= 0
                A(row+2,col+1) =   - q_f(j,c);
                b(row+2)       =     q_0(j,c);
            end
        end

        % equality constraints
        Aeq = zeros(2,k);
        beq = ones(2,1);
        for c = 1:k
            % equality constraints
            for j = 1:2
                Aeq(j,c) = q_f(j,c);
                beq(j)   = beq(j) - q_0(j,c);
            end
        end

        % objective function
        H = zeros(k, k);
        h = zeros(k, 1);
        for c = 1:k
            base = (c-1);
            forwardWeight  = (classWeights(c) * fbWeights(1));
            H(base+1,base+1) =  2/F(c)^2 * forwardWeight;
            h(base+1)        = -2/F(c)   * forwardWeight;
        end

        % solve optimization problem
        fF = solve();

        for c = 1:k
     %       fprintf('Fitting MAMAP(2,m) F+B: F(%d) = %f -> %f\n', c, F(c), fF(c));
        end

        % compute parameters of D11,D12,...,D1k
        q = zeros(3,k);
        for c = 1:k
            for j = 1:2
                q(j,c) = fF(c) * q_f(j,c) + q_0(j,c);
            end
        end
        for c = 1:k
            q(3,c) = 1/k;
        end
        
    else
        
      %  fprintf('Fitting MAMAP(2,m) F+B: fitting backward\n');
    
        % coefficients: q(j,c) = F(c) * q_b(j,c) + q_0(j,c)
        q_b = zeros(2,k);
        q_0 = zeros(2,k);
        for c = 1:k
            % q1
            q_b(1,c) = p(c) * ( -(r1 - 2)/((h2 + h1*(r1 - 1))*(r1 - 1)) );
            q_0(1,c) = p(c) * ( 1 - (h1 + h2)/((r1 - 1)*(h2 - h1 + h1*r1)) );
            % q2
            q_b(2,c) = p(c) * ( -(r1 - 2)/(h2 + h1*(r1 - 1)) );
            q_0(2,c) = p(c) * ( (h1*(r1 - 2))/(h2 - h1 + h1*r1) );
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
            backwardWeight  = (classWeights(c) * fbWeights(2));
            H(base+1,base+1) =  2/B(c)^2 * backwardWeight;
            h(base+1)        = -2/B(c)   * backwardWeight;
        end

        % solve optimization problem
        fB = solve();

        for c = 1:k
       %     fprintf('Fitting MAMAP(2,m) F+B: B(%d) = %f -> %f\n', c, B(c), fB(c));
        end

        % compute parameters of D11,D12,...,D1k
        q = zeros(3,k);
        for c = 1:k
            for j = 1:2
                q(j,c) = fB(c) * q_b(j,c) + q_0(j,c);
            end
        end
        for c = 1:k
            q(3,c) = 1/k;
        end
        
    end
    
else

    % coefficients: q(j,c) = F(c) * q_f(j,c) + B(c) * q_b(j,c) + q_0(j,c)
    q_f = zeros(3,k);
    q_b = zeros(3,k);
    q_0 = zeros(3,k);
    for c = 1:k
        if form == 1
            % first canonical form (positive auto-correlation decay)
            q_f(1,c) = 0;
            q_b(1,c) = -(p(c)*(r1*r2 - r2 + 1))/((h2 - h1*r2)*(r1 - 1)*(r2 - 1));
            q_0(1,c) = (p(c)*(h1 + h2 - h1*r2)*(r1*r2 - r2 + 1))/((h2 - h1*r2)*(r1 - 1)*(r2 - 1));
            q_f(2,c) = -(p(c)*(r1*r2 - r2 + 1))/(r1*(h1 + h2*(r1 - 1))*(r2 - 1));
            q_b(2,c) = -(p(c)*(r1*r2 - r2 + 1))/(r1*(h2 - h1*r2)*(r2 - 1));
            q_0(2,c) = (p(c)*(r1*r2 - r2 + 1))/((r1 - 1)*(r2 - 1)) + (h1*p(c)*(r1*r2 - r2 + 1))/(r1*(h2 - h1*r2)*(r2 - 1)) - (h1*p(c)*(r1*r2 - r2 + 1))/(r1*(h1 + h2*(r1 - 1))*(r1 - 1)*(r2 - 1));
            q_f(3,c) = -(p(c)*(r1*r2 - r2 + 1))/(r1*r2*(h1 - h2 + h2*r1));
            q_b(3,c) = 0;
            q_0(3,c) = (p(c)*(h1 + h2*r1)*(r1*r2 - r2 + 1))/(r1*r2*(h1 - h2 + h2*r1));
        else
            % second canonical form (negative auto-correlation decay)
            q_f(1,c) = 0;
            q_b(1,c) = -(p(c)*(r1 + r2 - r1*r2 - 2))/((r1 - 1)*(r2 - 1)*(h1 - h2 - h1*r1 + h1*r1*r2));
            q_0(1,c) = (p(c)*(h2 + h1*r1 - h1*r1*r2)*(r1 + r2 - r1*r2 - 2))/((r1 - 1)*(r2 - 1)*(h1 - h2 - h1*r1 + h1*r1*r2));
            q_f(2,c) = (p(c)*(r1 + r2 - r1*r2 - 2))/((r2 - 1)*(h1 - h2 + h2*r1));
            q_b(2,c) = 0;
            q_0(2,c) = -(h2*p(c)*(r1 + r2 - r1*r2 - 2))/((r2 - 1)*(h1 - h2 + h2*r1)); 
            q_f(3,c) = (p(c)*(r1 + r2 - r1*r2 - 2))/(r2*(h1 + h2*(r1 - 1)));
            q_b(3,c) = (p(c)*(r1 + r2 - r1*r2 - 2))/(r2*(h1 - h2 - h1*r1 + h1*r1*r2));
            q_0(3,c) = (h1*p(c)*(r1 + r2 - r1*r2 - 2))/(r2*(h1 + h2*(r1 - 1))*(r1 - 1)) - (h1*p(c)*(r1 + r2 - r1*r2 - 2))/(r2*(h1 - h2 - h1*r1 + h1*r1*r2)) - (p(c)*(r1 + r2 - r1*r2 - 2))/(r2*(r1 - 1));
        end
    end

    % inequality constraints
    A = zeros(6*k,2*k);
    b = zeros(6*k,1);
    for c = 1:k
        for j = 1:3
            row = ((c-1)*6+(j-1)*2);
            col = (c-1)*2;
            % q_j <= 1
            A(row+1,col+1) =     q_f(j,c);
            A(row+1,col+2) =     q_b(j,c);
            b(row+1)       = 1 - q_0(j,c);
            % q_j >= 0
            A(row+2,col+1) =   - q_f(j,c);
            A(row+2,col+2) =   - q_b(j,c);
            b(row+2)       =     q_0(j,c);
        end
    end

    % equality constraints
    Aeq = zeros(3,2*k);
    beq = ones(3,1);
    for c = 1:k
        % equality constraints
        for j = 1:3
            Aeq(j,(c-1)*2+1) = q_f(j,c);
            Aeq(j,(c-1)*2+2) = q_b(j,c);
            beq(j)           = beq(j) - q_0(j,c);
        end
    end

    % objective function
    H = zeros(2*k, 2*k);
    h = zeros(2*k, 1);
    for c = 1:k
        base = (c-1)*2;
        forwardWeight  = (classWeights(c) * fbWeights(1));
        backwardweight = (classWeights(c) * fbWeights(2));
        H(base+1,base+1) =  2/F(c)^2 * forwardWeight;
        H(base+2,base+2) =  2/B(c)^2 * backwardweight;
        h(base+1)        = -2/F(c)   * forwardWeight;
        h(base+2)        = -2/B(c)   * backwardweight;
    end

    % solve optimization problem
    x = solve();

    % feasible set of moments
    fF = zeros(k,1);
    fB = zeros(k,1);
    for c = 1:k
        fF(c) = x((c-1)*2+1);
        fB(c) = x((c-1)*2+2);
    end

    for c = 1:k
      %  fprintf('Fitting MAMAP(2,m) F+B: F(%d) = %f -> %f\n', c, F(c), fF(c));
      %  fprintf('Fitting MAMAP(2,m) F+B: B(%d) = %f -> %f\n', c, B(c), fB(c));
    end

    % compute parameters of D11,D12,...,D1k
    q = zeros(3,k);
    for c = 1:k
        for j = 1:3
            q(j,c) = fF(c) * q_f(j,c) + fB(c) * q_b(j,c) + q_0(j,c);
        end
    end
end

% compute D11,D12,...,D1k
if form == 1
    for c = 1:k
        mmap{2+c} = mmap{2} .* [q(1,c) 0; q(2,c) q(3,c)];
    end
else
    for c = 1:k
        mmap{2+c} = mmap{2} .* [0 q(1,c); q(2,c) q(3,c)];
    end
end

fF = mmap_forward_moment(mmap, 1);
fB = mmap_backward_moment(mmap, 1);

    function x = solve()
       % fprintf('Fitting MAMAP(2,m) F+B: running quadratic programming solver...\n');
        options = optimset('Algorithm','interior-point-convex','Display','none','MaxIter',3000);
        %[x,fx,xflag] = quadprog(H, h, A, b, Aeq, beq, [], [], [], options);        
        lb = 1e-6*ones( size(A,2),1);
        ub = 1e6*ones( size(A,2),1);
        [x,fx,xflag]=QP(H, h, A, b, Aeq, beq, lb, ub, options);
        if xflag ~= 1
            error('Quadratic programming solver failed: %d\n', xflag);
        end
        fit_error = fx + length(x);
        %fprintf('Fitting MAMAP(2,m) F+B: error = %e\n', fit_error);
    end

end
