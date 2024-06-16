function [lG,X,Q,method] = pfqn_nc(lambda,L,N,Z,varargin)
% [LG,X,Q,METHOD] = PFQN_NC(L,N,Z,VARARGIN)
%
% L: Service demand matrix
% N: Population vector
% Z: Think time vector
% varargin: solver options (e.g., SolverNC.defaultOptions)
%
% LG: Logarithm of normalizing constant
% X: System throughputs
% Q: Mean queue-lengths
% METHOD: Solution method

options = Solver.parseOptions(varargin, SolverNC.defaultOptions);
method = options.method;
% backup initial parameters
Rin = length(N);

if sum(N)==0 || isempty(N)
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
    line_warning(mfilename,'The model has no positive demands in any class.\n');
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

[M,R] = size(L);
X=[];Q=[];
method = options.method;
switch options.method
    case {'ca'}
        [~,lG] = pfqn_ca(L,N,sum(Z,1));
    case {'adaptive','default'}
        if M>1
            if R==1 || (R <= 3 && sum(N)<50)
                [~,lG] = pfqn_ca(L,N,sum(Z,1));
                method = 'ca';
            else
                if M>R
                    [~,lG] = pfqn_kt(L,N,sum(Z,1));
                    method = 'kt';
                else
                    [~,lG] = pfqn_le(L,N,sum(Z,1));
                    method = 'le';
                end
            end
        elseif sum(Z,1)==0 % single queue, no delay
            lG = -N*log(L)';
            method = 'exact';
        else % repairman model
            if N<10000
                %if max(abs(N-round(N))) < GlobalConstants.FineTol
                %[lG] = pfqn_comomrm(L,N,Z,1,options.tol);
                %else
                [~,lG] = pfqn_mmint2_gausslegendre(L,N,sum(Z,1));
                method = 'gleint';
                %end
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
            line_error(mfilename,sprintf('The %s method requires a model with a delay and a single queueing station.',options.method));
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
        [~,lG] = pfqn_mci(L,N,sum(Z,1),options.samples,'imci');
    case {'mva'}
        [~,~,~,~,lG] = pfqn_mva(L,N,sum(Z,1));
    case 'mom'
        if length(N)>1
            try
                [~,lG,X,Q] = pfqn_mom(L,N,Z);
            catch ME
                ME.getReport
                % java exception, probably singular linear system
                %line_warning(mfilename,'Numerical problems.');
                lG = NaN;
            end
        else
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
        end
    case {'exact'}
        if M>=R || sum(N)>10 || sum(Z)>0
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
            method = 'ca';
        else
            [~,lG] = pfqn_recal(L,N,sum(Z,1));% implemented with Z=0
            method = 'recal';
        end
    case {'comom'}
        if R>1
            try
                % comom has a bug in computing X, sometimes the
                % order is switched
                if M>1
                    %line_error(mfilename,'comom is temporarily disabled');
                    [~,lG,X,Q] = pfqn_comombtf_java(L,N,Z);
                else % use double precision script for M=1
                    [lG] = pfqn_comomrm(L,N,Z,1,options.tol);
                end
            catch ME
                ME.getReport
                % java exception, probably singular linear system
                %line_warning(mfilename,'Numerical problems.');
                lG = NaN;
            end
        else
            [~,lG] = pfqn_ca(L,N,sum(Z,1));
            method = 'ca';
        end
    case {'pana','panacea','pnc'}
        [~,lG] = pfqn_panacea(L,N,sum(Z,1));
        if isnan(lG)
            line_warning(mfilename,'Model is not in normal usage, panacea cannot continue.\n');
        end
    case 'propfair'
        [~,lG] = pfqn_propfair(L,N,sum(Z,1));
    case {'recal'}
        if sum(Z)>0
            line_error(mfilename,'RECAL is currently available only for models with non-zero think times.\n');
        end
        [~,lG] = pfqn_recal(L,N,sum(Z,1));
    case 'rgf'
        if sum(Z)>0
            line_error(mfilename,'RGF is defined only for models with non-zero think times.\n');
        end
        [~,lG] = pfqn_rgf(L,N);
    otherwise
        line_error(mfilename,sprintf('Unrecognized method: %s',options.method));
end
end