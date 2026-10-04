% Computes the normalizing constant of a product-form queueing network.
% [LG,X,Q,METHOD] = PFQN_NC(LAMBDA,L,N,Z,VARARGIN)
function [lG,X,Q,method] = pfqn_nc(lambda,L,N,Z,varargin)

options = Solver.parseOptions(varargin, SolverNC.defaultOptions);
method = 'exact'; % if early return is triggered
% backup initial parameters
Rin = length(N);

if any(N<0) || isempty(N)
    lG = -Inf;
    X = [];
    Q = [];
    return
end

if sum(N)==0
    lG = 0;
    X = [];
    Q = [];
    return
end

if isempty(lambda)
    lambda=0*N;
end

X=[]; Q=[];

% compute open class contributions
Qopen = [];
lGopen = 0;
Ut = zeros(1, size(L,1));
for i=1:size(L,1)
    Ut(i) = (1-lambda*L(i,:)');
    if isnan(Ut(i))
        Ut(i) = 0;
    end
    L(i,:) = L(i,:)/Ut(i);
    Qopen(i,:) = lambda.*L(i,:)/Ut(i);
    %lGopen = lGopen + log(Ut(i));
end
Qopen(isnan(Qopen))=0;
ocl = find(isinf(N));
% then erase open classes
N(isinf(N)) = 0;

% first remove empty classes
nnzClasses = find(N);
lambda = lambda(:,nnzClasses);
L = L(:,nnzClasses);
N = N(:,nnzClasses);
Z = Z(:,nnzClasses);

% see _kb/03-api-layer.md (pfqn_nc dispatch notes) for rationale
R = length(N);
scalevec = ones(1,R);
%switch options.method
%    case {'adaptive','comom','default'}
%        % no-op
%    otherwise
%end
for r=1:R
    scalevec(r) = max([L(:,r);Z(:,r)]);
end
%end
L = L ./ repmat(scalevec,size(L,1),1);
Z = Z ./ scalevec;

% remove stations with no demand
Lsum = sum(L,2);
Lmax = max(L,[],2);
demStations = find((Lmax./Lsum)>GlobalConstants.FineTol);
noDemStations = setdiff(1:size(L,1), demStations);
L = L(demStations,:);
if any(N((sum(L,1) + sum(Z,1)) == 0)>0) % if there is a class with jobs but L and Z all zero
    if options.verbose
        line_warning(mfilename,'The model has no positive demands in any class.\n');
    end
    if isempty(Z) || sum(Z(:))<options.tol
        lG = 0;
    else
        lG = - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    return
end

% update M and R
[M,R]=size(L);

% return immediately if degenerate case
if isempty(L) || sum(L(:))<options.tol % all demands are zero
    if isempty(Z) || sum(Z(:))<options.tol
        lG = lGopen;
    else
        lG = lGopen - sum(factln(N)) + sum(N.*log(sum(Z,1))) + N*log(scalevec)';
    end
    return
elseif M==1 && (isempty(Z) || sum(Z(:))<options.tol) % single node and no think time
    lG = factln(sum(N)) - sum(factln(N)) + sum(N.*log(L(1,:))) + N*log(scalevec)';
    return
elseif size(unique(L,'rows'),1)==1 && (isempty(Z) || sum(Z(:))<options.tol)  % M identical replicas
    lG = factln(sum(N)+M-1) - sum(factln(N)) + sum(N.*log(L(1,:))) + N*log(scalevec)' - factln(M-1);
    return
end

% determine contribution from jobs that permanently loop at delay
zeroDemandClasses = find(sum(L,1)<options.tol); % all jobs in delay
nonzeroDemandClasses = setdiff(1:R, zeroDemandClasses);

if isempty(sum(Z,1)) || all(sum(Z(:,zeroDemandClasses),1)<options.tol)
    lGzdem = 0;
    Nz = 0;
else
    if isempty(zeroDemandClasses) % for old MATLAB release compatibility
        lGzdem = 0;
        Nz = 0;
    else
        Nz = N(zeroDemandClasses);
        lGzdem = - sum(factln(Nz)) + sum(Nz.*log(sum(Z(:,zeroDemandClasses),1))) + Nz*log(scalevec(zeroDemandClasses))';
    end
end
L = L(:,nonzeroDemandClasses);
N = N(nonzeroDemandClasses);
Zz = Z(:,zeroDemandClasses);
Z = Z(:,nonzeroDemandClasses);
scalevecz = scalevec(nonzeroDemandClasses);
% compute G for classes No with non-zero demand. A method that supplies mean
% values by simulation ('mcmc') is told whether the caller asked for them: a
% single-output call wants lG only and must not pay for a simulation.
[lGnzdem,Xnnzdem,Qnnzdem,method] = compute_norm_const(L, N, Z, options, nargout>=2);

if isempty(Xnnzdem) % in this case the NC method does not return metrics as a by-product
    X = [];
    Q = [];
else
    zClasses = setdiff(1:Rin, nnzClasses);
    Xz = zeros(1,length(zClasses));
    Xnnz = zeros(1,length(nnzClasses));
    Xnnz(zeroDemandClasses) = Nz./ sum(Zz,1)./ scalevec(zeroDemandClasses);
    Xnnz(nonzeroDemandClasses) = Xnnzdem./ scalevec(nonzeroDemandClasses);
    X(1,[zClasses, nnzClasses]) = [Xz, Xnnz];
    X(ocl) = lambda(ocl);
    Qz = zeros(size(Qnnzdem,1),length(zClasses));
    Qnnz = zeros(size(Qnnzdem,1),length(nnzClasses));
    Qnnz(:,zeroDemandClasses) = 0; % they are all in the delay
    Qnnz(:,nonzeroDemandClasses) = Qnnzdem; % Q does not require scaling
    Q(noDemStations,:) = 0;
    Q(demStations,[zClasses, nnzClasses]) = [Qz, Qnnz];
    Q(:,ocl) = Qopen(:,ocl);
end
% scale back to original demands
lG = lGopen + lGnzdem + lGzdem + N*log(scalevecz)';
end

function [lG,X,Q,method] = compute_norm_const(L,N,Z,options,wantXQ)
% LG = COMPUTE_NORM_CONST(L,N,Z,OPTIONS,WANTXQ)
% Auxiliary script that computes LG after the initial filtering of L,N,Z.
% WANTXQ is false when the caller requested lG alone.

% Note: methods that can handle more efficiently replicas need to do so
% within the method function

% L,N,Z
[M,R] = size(L);
X=[];Q=[];
if nargin < 5, wantXQ = true; end
method = options.method;
switch options.method
    case {'ca'}
        [~,lG] = pfqn_ca(L,N,sum(Z,1));
    case {'divdiff'}
        % Divided-difference closed form, Casale (SIGMETRICS 2017), Eqs. (15)
        % and (16). Load-independent single-server queues only: a think time
        % needs the integral form of Corollary 3.4, which is not implemented.
        % Unlike the default route below this one keeps pfqn_explicit's
        % warnings, since a caller that named the method has no fallback.
        if sum(Z(:))>0
            line_error(mfilename,'The ''divdiff'' method requires a model without think time, which needs the integral form of Corollary 3.4. Use ''ca'' or ''default''.');
        end
        [lG,~,expr] = pfqn_explicit(L,N);
        method = ['divdiff/',expr];
    case {'clw'}
        % Choudhury-Leung-Whitt generating function inversion: each
        % single-server station is a multiplicity-1 queue, delay is the IS term
        [~,lG] = pfqn_clw(L,N,sum(Z,1));
    case {'ger'}
        % Residue closed form of the same generating function 'clw' inverts
        % numerically. A class eliminated by residues enters only as a pole
        % ORDER, so its population is free: this is the cheap route when one
        % population dwarfs the others, and the expensive one when the classes
        % are many, since the term count grows as C(S+M-1,M-1) per further
        % elimination. options.maxterms REFUSES rather than truncating, so an
        % oversized model errors here instead of returning a wrong lG. The
        % solver options are deliberately not forwarded: pfqn_gerasimov's tol is
        % a pole-merging threshold, not the iterative tolerance options.tol
        % carries.
        [~,lG] = pfqn_gerasimov(L,N,sum(Z,1));
    case {'adaptive','default'}
        % ONE ESTIMATOR ANSWERS THE WHOLE FAMILY. The divided-difference closed
        % form of Casale (SIGMETRICS 2017) is exact here and was briefly tried
        % first on M>1 && R==1 && sum(Z)==0, but the default route does not
        % serve a single constant: the analyzer differences it at N-e_r for X
        % and at the AUGMENTED shape for Q, one extra class holding one job at
        % station i. That shape has R+1 classes, which the closed form refuses
        % at any sizeable population (the outer sum's cancellation), so it kept
        % the cubature while G(N) turned exact. Mixing the two costs more than
        % either: on mqn_singleserver_ps the closed-form G(N) under cubature
        % numerators left sum_i Q_i at 99.500 of N=100, and the conservation
        % rescale then moved the entire cubature error into X, 0.5% against the
        % 0.06% the cubature ratio carries on its own. 'divdiff' stays a NAMED
        % method, where the caller owns the whole family.
        if M>1
            order = -1;
            if sum(N)<1e3
                % see _kb/03-api-layer.md (pfqn_nc dispatch notes) for rationale
                Cmax      = M*R*(50)^3;                  % upper cost budget
                maxorder  = min(ceil((sum(N)-1)/2),16);

                totCost   = 0;
                order     = 0;                         % will be raised as far as possible
                while order < maxorder
                    nextCost = R * nchoosek(M + 2*(order+1), M-1);  % cost of order+1
                    if totCost + nextCost <= Cmax
                        order    = order + 1;
                        totCost  = totCost + nextCost;
                    else
                        break
                    end
                end
            end
            % Cmax prices neither the Grundmann-Moeller node count nor the
            % think-time v-integration, so the order is re-priced here and
            % lowered until it fits. Lowering the order keeps the cubature;
            % switching to le instead would hand these models to a Laplace
            % expansion whose mode sits on the simplex boundary when L has
            % near-zero rows, which is the flat-layer case that trips the budget
            while order > 0 && pfqn_cub_evals(M,order,sum(Z,1)) > GlobalConstants.CubMaxEvals
                order = order - 1;
            end
            if order >= 0
                [~,lG] = pfqn_cub(L,N,sum(Z,1),order,GlobalConstants.FineTol);
                method = 'cub';
            else
                % BLE on the default path: the correction is strictly better on
                % lG and cancels in G(N-e_r)/G(N). 'le' stays the published form.
                [~,lG] = pfqn_ble(L,N,sum(Z,1));
                method = 'ble';
                % Birman-Kogan Algorithm 2 supplies the MEAN VALUES here. The
                % caller's fallback differences lG at R+M*R reduced populations,
                % which on many stations is both dearer and ~300x less accurate
                % than the load concealment fixed point. Gated on the station count,
                % since load concealment is mean field in M: see
                % _kb/06-solver-catalog.md for the measured 10-station cliff.
                if M >= 10 && R > 1 && all(N >= 0) && any(N > 0)
                    [X,Q] = pfqn_bklc(L,N,sum(Z,1),'mva',[],options.iter_max);
                    method = 'ble/lc';
                end
            end
        elseif sum(Z(:))==0 % single queue, no delay
            lG = -N*log(L)';
            method = 'exact';
        else % repairman model
            if N<10000
                % gleint is a better method but there are JAR loading issues
                % at times upon loading the txt files
                %[~,lG] = pfqn_mmint2_gausslegendre(L,N,sum(Z,1));
                %method = 'gleint';
                [lG] = pfqn_comomrm(L,N,Z,1,options.tol);
                method = 'comom';
            else
                [~,lG] = pfqn_ble(L,N,sum(Z,1));
                method = 'ble';
            end
        end
    case {'sampling'}
        if M==1
            [~,lG] = pfqn_mmsample2(L,N,sum(Z,1),options.samples);
            method = 'sampling';
        elseif M>R
            [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,'imci');
            method = 'imci';
        else
            [~,lG] = pfqn_ls(L,N,sum(Z,1),options.samples);
            method = 'ls';
        end
    case {'mmint2','gleint'}
        if size(L,1)>1
            if options.verbose
                line_warning(mfilename,sprintf('The %s method requires a model with a delay and a single queueing station.',options.method));
            end
            lG = [];
            return
        end
        [~,lG] = pfqn_mmint2_gausslegendre(L,N,sum(Z,1));
    case {'cub','gm'} % Grundmann-Mueller cubatures
        order = ceil((sum(N)-1)/2); % exact
        [~,lG] = pfqn_cub(L,N,sum(Z,1),order,GlobalConstants.FineTol);
    case 'kt'
        [~,lG] = pfqn_kt(L,N,sum(Z,1));
    case 'bkt' % KT minus the exact Stirling remainder of each Laplaced class; see _kb/03-api-layer.md
        [~,lG] = pfqn_bkt(L,N,sum(Z,1));
    case 'lekt' % the estimator ble and bkt both compute, on the cheaper side; see _kb/03-api-layer.md
        [~,lG] = pfqn_lekt(L,N,sum(Z,1));
    case 'bk'
        [~,lG] = pfqn_bk(L,N,sum(Z,1));
    case 'bkue'
        % The uniform expansion is single chain by construction: it keeps one
        % dominant pole and the saddle point in a single erfc formula. The
        % multichain fallback is the saddle point of the same paper, which is
        % also how the analyzer reaches this branch, since it conditions on a
        % station population by augmenting the model with an auxiliary class.
        if size(L,2) > 1
            [~,lG] = pfqn_bk(L,N,sum(Z,1));
            method = 'bkue/bk';
        else
            [~,lG] = pfqn_bkue(L,N,sum(Z,1));
        end
    case {'lc','lc.ue'}
        if strcmpi(options.method,'lc.ue')
            inner = 'ue';
        else
            inner = 'mva';
        end
        % the fixed point converges linearly and slowly, so a solver-level
        % reporting tolerance would stop it far from its own limit and at a
        % different sweep in each codebase: iterate to the method's accuracy
        [X,Q] = pfqn_bklc(L,N,sum(Z,1),inner,[],options.iter_max);
        % Algorithm 2 returns mean values, not a multichain constant; the
        % saddle point that seeds it supplies lG on the same asymptotics
        [~,lG] = pfqn_bk(L,N,sum(Z,1));
    case {'mcmc'}
        % Chen-O'Cinneide REGULARIZATION (TOMACS 8(3), 1998). The chain is
        % simulated on the regularized network, which shares the steady-state
        % distribution of the original one, so it returns X and Q directly
        % through the pfqn_nc X/Q channel, like 'lc'. What it does NOT
        % return is the constant itself: the algorithm estimates the RATIOS
        % G(N-e_r)/G(N), never G, so lG here is the BLE expansion and is not
        % part of the paper. It cancels out of every mean value reported by the
        % analyzer; only getProbNormConstAggr reads it.
        if wantXQ
            [X,Q] = pfqn_mcmc(L,N,sum(Z,1),[],options);
        end
        if M > 1
            [~,lG] = pfqn_ble(L,N,sum(Z,1));
        else
            lG = pfqn_comomrm(L,N,Z,1,options.tol);
        end
    case 'le'
        [~,lG] = pfqn_le(L,N,sum(Z,1));
    case 'ble' % LE plus the empirical eps->0 correction; see _kb/03-api-layer.md
        [~,lG] = pfqn_ble(L,N,sum(Z,1));
    case 'aghq' % adaptive Gauss-Hermite over the simplex; q=1 is 'le'
        % options.config.aghq_nodes: nodes per simplex direction. The rule costs
        % q^(M-1) evaluations, so the default stays small.
        aghqNodes = 3;
        if isfield(options,'config') && isfield(options.config,'aghq_nodes') ...
                && ~isempty(options.config.aghq_nodes)
            aghqNodes = max(1,round(options.config.aghq_nodes));
        end
        [~,lG] = pfqn_aghq(L,N,sum(Z,1),aghqNodes);
    case 'ls'
        [~,lG] = pfqn_ls(L,N,sum(Z,1),options.samples);
    case {'is'}
        % see _kb/03-api-layer.md (pfqn_nc dispatch notes) for rationale
        [~,lG] = pfqn_is(L,N,sum(Z,1),options);
        method = 'is';
    case {'mci','imci'}
        [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,options.method);
    case {'mva'}
        [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
    case {'exact'}
        if M>=R || sum(N)>10 || sum(Z(:))>0
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
            method = 'exact/ca';
        else
            [~,lG] = pfqn_recal(L,N,sum(Z,1));% implemented with Z=0
            method = 'exact/recal';
        end
    case {'comom'}
        if R>1
            % see _kb/03-api-layer.md (pfqn_nc dispatch notes) for rationale
            if M>1
                line_error(mfilename,'The ''comom'' method supports a single queueing station, but this model has %d. Use ''default'' or ''ca'' for an exact normalizing constant, or SolverJMT with method ''jmva.recal''.', M);
            end
            try
                % comom has a bug in computing X, sometimes the
                % order is switched
                [lG] = pfqn_comomrm(L,N,Z,1,options.tol);
            catch ME
                getReport(ME,'basic')
                % java exception, probably singular linear system
                %if options.verbose
                %line_warning(mfilename,'Numerical problems.');
                %end
                lG = [];
            end
        else
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
            method = 'ca';
        end
    case {'pana'}
        [~,lG] = pfqn_panacea(L,N,sum(Z,1));
        if isnan(lG)
            % Outside normal usage the {phi(n)} series diverges and
            % pfqn_panacea answers NaN. Bail out at ANY verbosity: solver_nc
            % reads an empty constant as out-of-domain and returns an empty
            % table, whereas a NaN reaches the reported table unnoticed. This
            % is the run-path half of the rule nc_method_refusal states for the
            % report through nc_is_normal_usage.
            if options.verbose
                line_warning(mfilename,'Model is not in normal usage, pana cannot continue.\n');
            end
            lG = [];
            return
        end
    case 'propfair'
        [~,lG] = pfqn_propfair(L,N,sum(Z,1));
    case {'recal'}
        % This arm is stated for models WITHOUT think time; the 'exact' route
        % above reaches pfqn_recal only when sum(Z)==0. Bail out at ANY
        % verbosity, as on the 'mmint2'/'gleint' arm: an empty constant is what
        % solver_nc reads as out-of-domain, while falling through would
        % evaluate the recursion off the shape this arm states.
        if sum(Z(:))>0
            if options.verbose
                line_warning(mfilename,'RECAL is currently available only for models without think times.\n');
            end
            lG = [];
            return
        end
        [~,lG] = pfqn_recal(L,N,sum(Z,1));
    case {'rgf'}
        % Recursion by generating functions. Single class: one sequence per
        % group of identically loaded stations (Coury-Harrison 1997, Property 1).
        % Multiclass: the residue recursion of Harrison-Coury 2002, Thm 1.
        if R > 1
            % Multiclass goes through the residue elimination of Harrison-Coury
            % Thm 1, with think times carried by the Bertozzi-McKenna truncation
            % that neither RGF paper has. That sum is ALTERNATING, so pfqn_rgfmc
            % refuses when the cancellation leaves no significant digits rather
            % than returning a wrong lG; the exact convolution answers those and
            % the reported method says so.
            try
                [~,lG] = pfqn_rgfmc(L,N,sum(Z,1));
            catch
                [~,lG] = pfqn_ca(L,N,sum(Z,1));
                method = 'rgf/ca';
            end
        else
            [~,lG] = pfqn_rgf(L,N,sum(Z,1));
        end
    otherwise
        % An unrecognized token must not degrade into an empty normalizing
        % constant: the caller reports the analysis as completed and hands
        % back an empty table, so the refusal is invisible.
        line_error(mfilename,sprintf('Unrecognized method: %s',options.method));
end
end