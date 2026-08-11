%{ @file ctmc_relsolve.m
 %  @brief Computes the equilibrium distribution relative to a reference state
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes the equilibrium distribution relative to a reference state
 %
 % @details
 % Solves the global balance equations with the normalization condition p(refstate) = 1.
 %
 % @par Syntax:
 % @code
 % p = ctmc_relsolve(Q)
 % [p, Q, nConnComp, connComp] = ctmc_relsolve(Q, refstate)
 % [p, Q, nConnComp, connComp] = ctmc_relsolve(Q, refstate, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Infinitesimal generator matrix
 % <tr><td>refstate<td>(Optional) Index of the reference state (default: 1)
 % <tr><td>options<td>(Optional) Solver options
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>p<td>Relative equilibrium distribution vector
 % <tr><td>Q<td>Processed generator matrix
 % <tr><td>nConnComp<td>Number of connected components
 % <tr><td>connComp<td>Component assignment vector
 % </table>
%}
function [p, Q, nConnComp, connComp]=ctmc_relsolve(Q,refstate,options)

if nargin<2
    refstate = 1;
end

if length(Q) > 6000 && (nargin==1 || ~options.force)
    fprintf(1,'ctmc_relsolve: the order of Q is large (%d). Press key to continue.\n',length(Q));
    pause;
end

if size(Q)==1
    p = 1;
    return
end

Q = ctmc_makeinfgen(Q); % so that spurious diagonal elements are set to 0
n = length(Q);

if issym(Q)
    symvariables = symvar(Q); % find all symbolic variables
    B = double(subs(Q+Q',symvariables,ones(size(symvariables)))); % replace all symbolic variables with 1.0
else
    B = abs(Q+Q')>0;
end
[nConnComp, connComp] = weaklyconncomp(B);

if nConnComp > 1
    % reducible generator - solve each component recursively
    if issym(Q)
        p = sym(zeros(1,n));
    else
        p = zeros(1,n);
    end
    for c=1:nConnComp
        Qc = Q(connComp==c,connComp==c);
        Qc = ctmc_makeinfgen(Qc);
        p(connComp==c) = ctmc_solve(Qc);
    end
    p = p /sum(p);
    return
end

if all(Q==0)
    p = ones(1,n)/n;
    return
end
p = zeros(1,n);
b = zeros(n,1);

nnzel = 1:n;
Qnnz = Q; bnnz = b;
Qnnz_1 = Qnnz; bnnz_1 = bnnz;

isReducible = false;
goon = true;
while goon

    %
    %     zerorow=find(sum(abs(Qnnz),2)==0);
    %         if length(zerorow)>=1
    %             if nargin==1 || options.verbose
    %                 %warning('ctmc_solve: the infinitesimal generator is reducible (zero row)');
    %                 fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
    %                 isReducible = true;
    %             end
    %         end
    %     nnzrow = setdiff(nnzel, zerorow);
    %
    %     zerocol=find(sum(abs(Qnnz),1)==0);
    %     nnzcol = setdiff(nnzel, zerocol);
    %         if length(zerocol)>=1
    %             if ~isReducible && (nargin==1 || options.verbose)
    %                 %warning('ctmc_solve: the infinitesimal generator is reducible (zero column)');
    %                 fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
    %             end
    %         end
    %     nnzel = intersect(nnzrow, nnzcol);

    nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0);
    if length(nnzel) < n && ~isReducible
        isReducible = true;
        if (nargin > 1 && options.verbose == 2) % debug
            fprintf(1,'ctmc_solve: the infinitesimal generator is reducible.\n');
        end
    end
    Qnnz = Qnnz(nnzel, nnzel);
    bnnz = bnnz(nnzel);
    Qnnz = ctmc_makeinfgen(Qnnz);
    if all(size(Qnnz_1(:)) == size(Qnnz(:))) && all(size(bnnz_1(:)) == size(bnnz(:)))
        goon = false;
    else
        Qnnz_1 = Qnnz; bnnz_1 = bnnz; nnzel = 1:length(Qnnz);
    end
end

if isempty(Qnnz)
    p = ones(1,n)/n;
    return
end
Qnnz(:,end) = 0;
Qnnz(refstate,end) = 1;
bnnz(end) = 1;

if ~isdeployed
    if issym(Q)
        p = sym(p);
    end
end

if nargin == 1
    p(nnzel)=Qnnz'\ bnnz;
else
    switch options.method
        case 'gpu'
            try
                gQnnz = gpuArray(Qnnz');
                gbnnz = gpuArray(bnnz);
                pGPU = gQnnz \ gbnnz;
                gathered_pGPU = gather(pGPU);
                p(nnzel) = gathered_pGPU; % transfer from GPU to local env
            catch
                warning('ctmc_solve: GPU either not available or execution failed. Switching to default method.');
                p(nnzel) = Qnnz'\ bnnz;
            end
        otherwise
            p(nnzel)=Qnnz'\ bnnz;
    end
end

if issym(Q)
    Q=simplify(Q);
end

end
