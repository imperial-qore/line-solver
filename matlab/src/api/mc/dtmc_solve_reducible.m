function [pi,pis,pi0,scc,isrec,Pl,pil] = dtmc_solve_reducible(P, pin, options)
% [PI,PIS,PI0,SCC,ISREC,PL,PIL] = DTMC_SOLVE_REDUCIBLE(P, PIN, OPTIONS)
%
% Estimate limiting distribution for a DTMC P that may have reducible
% components
%
% Input:
% P: dtmc transition matrix
% pin: initial vector, to be set to [] if not available
% options: struct where options.tol sets the tolerance
%
% Output:
% pi: - For an ergodic DTMC, this is the unique limiting distribution.
%     - For a reducible DTMC:
%           - if there is a single transient SCC then this
%             is assumed to be the starting state with uniform probability
%             across its states.
%           - if there are multiple transient SCCs then pi is the average
%             of the limiting distributions in pis
%
% pis: limiting distribution given initialization pi0 in a single SCC
% pi0: start vector of the lumped DTMC for each row of pis. For ergodic
%      DTMCs, this is the empty vector.
% scc: mapping of state of P to SCCs
% isrec: element i is true if SCC i is recurrent or false otherwise
% Pl: lumped DTMC where each SCC is replaced by a single state
% pil: limiting distribution of the lumped DTMC
%
% WARNING: the script does not consider explicitly periodic SCCs

[scc, isrec] = stronglyconncomp(P);

% Number of SCCs
numSCC = max(scc);
if numSCC==1
    pi = dtmc_solve(P);
    pis = pi;
    pi0 = [];
    Pl = P;
    pil = pi;
    return
else
    % Lump SCCs
    Pl = zeros(numSCC);

    % Lump states according to SCCs
    scc_idx = cell(1,numSCC);
    for i = 1:numSCC
        scc_idx{i} = scc == i;
    end

    for i = 1:numSCC
        scc_i = scc_idx{i};
        for j = 1:numSCC
            if i ~= j
                scc_j = scc_idx{j};
                Pl(i,j) = sum(P(scc_i, scc_j),"all");
            end
        end
    end

    Pl = dtmc_makestochastic(Pl);

    % Ensure that recurrent SCCs have self-loops
    for i = 1:numSCC
        if isrec(i)
            Pl(i,i) = 1;
        end
    end

    % acompute pinl(i), the probability of starting in scc i
    if nargin<2 || isempty(pin)
        pin = [];
        pinl = ones(1,numSCC);
        % Find columns (states) with near-zero column sums (absorbing states)
        z_cols = find(sum(P,1) < 1e-12);
        % Map state indices to SCC indices for zeroing pinl
        for j = z_cols
            pinl(scc(j)) = 0;
        end
        pinl = pinl / sum(pinl);
    else
        % in this case, an initial state vector is provided
        % lump initial vector
        pinl = zeros(1,numSCC);
        for i = 1:numSCC
            pinl(i) = sum(pin(scc_idx{i}));
        end
    end

    % using pinl, compute all limiting probabilities and average them
    pi = zeros(1,length(P));

    % Compute limiting matrix PI via spectral decomposition
    % First check if matrix is numerically stable for spectral decomposition
    % by examining the eigenvector matrix condition number
    useSpectral = true;
    try
        [V, ~] = eig(Pl);
        condV = cond(V);
        if condV > 1e10 || isnan(condV) || isinf(condV)
            % Matrix has repeated eigenvalues or is numerically degenerate
            % Use power method directly to avoid slow/failing spectd
            useSpectral = false;
        end
    catch
        useSpectral = false;
    end

    if useSpectral
        try
            % Suppress singular matrix warnings during spectd call
            warnState = warning('off', 'MATLAB:singularMatrix');
            warnState2 = warning('off', 'MATLAB:nearlySingularMatrix');
            [eigs,projectors]=spectd(Pl);
            warning(warnState);
            warning(warnState2);
            PI = zeros(size(projectors{1}));
            for e=1:length(eigs)
                if eigs(e)>1-1e-12
                    PI = PI + projectors{e};
                end
            end
            % Check if PI has NaN values (spectd failed on singular matrix)
            if any(isnan(PI(:)))
                PI = compute_limiting_matrix_power(Pl);
            end
        catch
            % Fallback to power method if spectd fails
            PI = compute_limiting_matrix_power(Pl);
        end
    else
        % Use power method for numerically degenerate matrices
        PI = compute_limiting_matrix_power(Pl);
    end
    pi0 = zeros(numSCC,numSCC);
    pil = zeros(numSCC,numSCC);
    pis = zeros(numSCC,length(P));
    for i = 1:numSCC
        if pinl(i)>0
            pi0(i,:) = zeros(1,numSCC);
            pi0(i,i) = 1;
            pil(i,:) = pi0(i,:)*PI;
            pis(i,1:length(P)) = 0;
            for j=1:numSCC
                scc_j = scc_idx{j};
                if nargin<3
                    pis(i, scc_j) = pil(i,j)*dtmc_solve(P(scc_j, scc_j));
                else
                    pis(i, scc_j) = pil(i,j)*dtmc_solve(P(scc_j, scc_j), options);
                end
            end
            pi = pi + pis(i,:) * pinl(i);
        end
    end

    transient_states = find(isrec==false);
    if isscalar(transient_states) && isempty(pin)
        pi = pis(transient_states,:);
    else
        % leave pi as it is
    end
end
end

function PI = compute_limiting_matrix_power(P)
% COMPUTE_LIMITING_MATRIX_POWER Compute limiting matrix using power iteration
%
% For ergodic DTMCs, P^n converges to a matrix where each row is the
% stationary distribution. This function serves as a fallback when
% spectral decomposition fails due to singular or defective matrices.

n = size(P, 1);
Pk = P;
maxIter = 1000;
tol = 1e-10;

for iter = 1:maxIter
    Pk1 = Pk * P;

    % Check convergence
    maxDiff = max(abs(Pk1(:) - Pk(:)));
    Pk = Pk1;

    if maxDiff < tol
        break
    end
end

PI = Pk;
end
