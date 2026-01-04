%{ @file sn_get_demands_chain.m
 %  @brief Computes aggregated demand parameters per chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes aggregated demand parameters per chain
 %
 % @details
 % This function aggregates class-level service demands and related parameters
 % to chain-level parameters, computing total demands, service times, visits,
 % and aggregation factors.
 %
 % @par Syntax:
 % @code
 % [Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Lchain<td>Total service demand per station per chain
 % <tr><td>STchain<td>Mean service time per station per chain
 % <tr><td>Vchain<td>Mean number of visits per station per chain
 % <tr><td>alpha<td>Aggregation factors mapping classes to chains
 % <tr><td>Nchain<td>Population of each chain
 % <tr><td>SCVchain<td>Squared coefficient of variation of service times per station per chain
 % <tr><td>refstatchain<td>Reference station index for each chain
 % </table>
%}
function [Lchain,STchain,Vchain,alpha,Nchain,SCVchain,refstatchain] = sn_get_demands_chain(sn)
M = sn.nstations;    %number of stations
%K = sn.nclasses;    %number of classes
C = sn.nchains;
N = sn.njobs';  % initial population per class

%PH = sn.proc;
SCV = sn.scv;
SCV(isnan(SCV))=1;
% determine service times
ST = 1./sn.rates;
ST(isnan(ST))=0;

alpha = zeros(sn.nstations,sn.nclasses);
Vchain = zeros(sn.nstations,sn.nchains);
for c=1:sn.nchains
    inchain = sn.inchain{c};    
    if sn.refclass(c)>0 % if the model has a ref class
        for i=1:sn.nstations
            Vchain(i,c) = sum(sn.visits{c}(sn.stationToStateful(i),inchain)) / sum(sn.visits{c}(sn.stationToStateful(sn.refstat(inchain(1))),sn.refclass(c)));
            for k=inchain
                alpha(i,k) = alpha(i,k) + sn.visits{c}(sn.stationToStateful(i),k) / sum(sn.visits{c}(sn.stationToStateful(i),inchain));
            end
        end
    else
        for i=1:sn.nstations
            Vchain(i,c) = sum(sn.visits{c}(sn.stationToStateful(i),inchain)) / sum(sn.visits{c}(sn.stationToStateful(sn.refstat(inchain(1))),inchain));
            for k=inchain
                alpha(i,k) = alpha(i,k) + sn.visits{c}(sn.stationToStateful(i),k) / sum(sn.visits{c}(sn.stationToStateful(i),inchain));
            end
        end
    end
end

Vchain(~isfinite(Vchain))=0;
for c=1:sn.nchains
    inchain = sn.inchain{c};
    Vchain(:,c) = Vchain(:,c) / Vchain(sn.refstat(inchain(1)),c);
end
alpha(~isfinite(alpha))=0;
alpha(alpha<GlobalConstants.Zero)=0;

Lchain = zeros(M,C);
STchain = zeros(M,C);
SCVchain = zeros(1,C);
Nchain = zeros(1,C);
refstatchain = zeros(C,1);
for c=1:sn.nchains
    inchain = sn.inchain{c};
    Nchain(c) = sum(N(inchain));
    isOpenChain = any(isinf(N(inchain)));
    for i=1:sn.nstations
        % we assume that the visits in L(i,inchain) are equal to 1
        STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        if isOpenChain && i == sn.refstat(inchain(1)) % if this is a source ST = 1 / arrival rates
            STchain(i,c) = 1 / sumfinite(sn.rates(i,inchain)); % ignore degenerate classes with zero arrival rates
        else
            STchain(i,c) = ST(i,inchain) * alpha(i,inchain)';
        end
        Lchain(i,c) = Vchain(i,c) * STchain(i,c);
        alphachain = sum(alpha(i,inchain(isfinite(SCV(i,inchain))))');
        if alphachain>0
            SCVchain(i,c) = SCV(i,inchain) * alpha(i,inchain)' / alphachain;
        end
    end
    refstatchain(c) = sn.refstat(inchain(1));
    if any((sn.refstat(inchain(1))-refstatchain(c))~=0)
        line_error(mfilename,sprintf('Classes in chain %d have different reference station.',c));
    end
end
Lchain(~isfinite(Lchain))=0;
STchain(~isfinite(STchain))=0;
end