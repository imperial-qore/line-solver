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

% then scale demands in [0,1], importat that stays before the other
% simplications in case both D and Z are all very small or very large in a
% given class, in which case the may look to filter but not if all of them
% are at the same scale
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
% compute G for classes No with non-zero demand
[lGnzdem,Xnnzdem,Qnnzdem,method] = compute_norm_const(L, N, Z, options);

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

function [lG,X,Q,method] = compute_norm_const(L,N,Z,options)
% LG = COMPUTE_NORM_CONST(L,N,Z,OPTIONS)
% Auxiliary script that computes LG after the initial filtering of L,N,Z

% Note: methods that can handle more efficiently replicas need to do so
% within the method function

% L,N,Z
[M,R] = size(L);
X=[];Q=[];
method = options.method;
switch options.method
    case {'ca'}
        [~,lG] = pfqn_ca(L,N,sum(Z,1));
    case {'adaptive','default'}
        if M>1
            if sum(N)<1e3
                % Sequence of Grundmann-Mueller cubature (CUB) executions
                %
                % This is slowish since CUB is O(N^M) but competitors fails
                % 'imci' has occasionally large errors in the benchmarks and
                % 'kt' fails with self-looping customers leading to
                % degradations see eg in example_mixedModel_4
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
                
                [~,lG] = pfqn_cub(L,N,sum(Z,1),order,GlobalConstants.FineTol);
                method = 'cub';
            else
                [~,lG] = pfqn_le(L,N,sum(Z,1));
                method = 'le';
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
                [~,lG] = pfqn_le(L,N,sum(Z,1));
                method = 'le';
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
    case 'le'
        [~,lG] = pfqn_le(L,N,sum(Z,1));
    case 'ls'
        [~,lG] = pfqn_ls(L,N,sum(Z,1),options.samples);
    case {'mci','imci'}
        [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,options.method);
    case {'mva'}
        [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
        %case 'mom'
        %    line_warning(mfilename,'mom supported only in SolverJMT (method ).');
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
            try
                % comom has a bug in computing X, sometimes the
                % order is switched
                if M>1
                    if options.verbose
                        line_warning(mfilename,'comom supported only for M=1 queue. Use SolverJMT instead (method jmva.comom).');
                    end
                    lG = [];
                    return
                else % use double precision script for M=1
                    [lG] = pfqn_comomrm(L,N,Z,1,options.tol);
                end
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
    case {'panacea'}
        [~,lG] = pfqn_panacea(L,N,sum(Z,1));
        if isnan(lG)
            if options.verbose
                line_warning(mfilename,'Model is not in normal usage, panacea cannot continue.\n');
                lG = [];
                return
            end
        end
    case 'propfair'
        [~,lG] = pfqn_propfair(L,N,sum(Z,1));
    case {'recal'}
        if sum(Z(:))>0
            if options.verbose
                line_warning(mfilename,'RECAL is currently available only for models with non-zero think times.\n');
                lG = [];
                return
            end
        end
        [~,lG] = pfqn_recal(L,N,sum(Z,1));
        % case 'rgf'
        %     if sum(Z(:))>0
        %         if option.verbose
        %             line_warning(mfilename,'RGF is defined only for models with non-zero think times.\n');
        %             lG = [];
        %             return
        %         end
        %     end
        %     [~,lG] = pfqn_rgf(L,N);
    otherwise
        lG=[];
        X=[];
        Q=[];
        if options.verbose
            line_warning(mfilename,sprintf('Unrecognized method: %s',options.method));
            lG = [];
            return
        end
        return
end
end