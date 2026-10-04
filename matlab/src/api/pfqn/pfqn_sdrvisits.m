%{
%{
 % @file pfqn_sdrvisits.m
 % @brief Section 3.2 coefficients xi of a network with state-dependent routing.
%}
%}

%{
%{
 % @brief Section 3.2 coefficients xi of a network with state-dependent routing.
 % @fn pfqn_sdrvisits(sdr, P)
 % @param sdr State-dependent routing structure.
 % @param P State-independent routing probabilities per chain.
 % @return xi Coefficients xi_ij of the product form.
%}
%}
function xi = pfqn_sdrvisits(sdr, P)
% XI = PFQN_SDRVISITS(SDR, P)
%
% Coefficients xi_ij of Krzesinski (1987), Section 3.2, for a network with
% state-dependent routing. P is an MxMxJ array of the state-independent (SIR)
% routing probabilities: P(x,y,j) is the probability that a chain j customer
% leaving center x proceeds to center y. The state-dependent arcs out of the
% entry center e of Q(V,V) are not part of P and are ignored if present.
%
% The coefficients are fixed by three rules:
%   - the entry center e and the departure center d of Q(V,V) satisfy
%     xi_ej = xi_dj, and the centers of the complement M-V obey the ordinary
%     traffic equations in which the whole SDR subnetwork acts as a single arc
%     from e to d carrying probability one;
%   - the entry center e(b) and the departure center d(b) of every branch, and
%     every center inside a branch, obey the branch's own traffic equations
%     driven by an injection of xi_ej at e(b). Because a customer leaves a
%     branch only through d(b), this returns xi_{d(b)j} = xi_ej, and it
%     returns xi_{e(b)j} = xi_ej whenever e(b) receives no internal feedback,
%     which covers every single-center branch. The paper states the identity
%     xi_ij = xi_ej for the branch entry and departure centers and works out
%     only single-center branches; the traffic equations above are the reading
%     that extends it to a branch holding several centers;
%   - the normalization xi_ej = 1.
%
% These xi are not relative visit counts. Under SDR the rate at which
% customers enter a branch depends on the network state, so the ratio of two
% xi carries no flow interpretation; they are the solution of the transformed
% balance equations of Appendix A.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

c = pfqn_sdrcoeff(sdr);
M = size(P,1);
if size(P,2) ~= M
    line_error(mfilename,'The SIR routing array must be square in its first two dimensions.');
end
J = size(P,3);
xi = zeros(M,J);

inV = false(1,M);
for b = 2:c.B
    inV(c.branch{b}) = true;
end
mv = find(~inV);
ie = find(mv == c.entry, 1);
id = find(mv == c.departure, 1);
if isempty(ie) || isempty(id)
    line_error(mfilename,'The entry and departure centers of Q(V,V) must lie outside every branch.');
end

for j = 1:J
    % Complement M-V with the SDR subnetwork collapsed into the single arc e->d
    Pmv = P(mv, mv, j);
    Pmv(ie,:) = 0;
    Pmv(ie,id) = 1;
    rs = sum(Pmv,2);
    if any(abs(rs - 1) > GlobalConstants.CoarseTol)
        line_error(mfilename,sprintf('The SIR routing of chain %d does not keep customers inside the complement M-V.',j));
    end
    xmv = dtmc_solve(Pmv);
    if xmv(ie) <= 0
        line_error(mfilename,sprintf('The entry center of Q(V,V) is unreachable in chain %d.',j));
    end
    xmv = xmv / xmv(ie);
    xi(mv,j) = xmv(:);

    % Each branch, driven by an injection of xi_e at its entry center
    for b = 2:c.B
        sb = c.branch{b}(:)';
        Pbb = P(sb, sb, j);
        inj = zeros(1,numel(sb));
        inj(sb == c.entryOf(b)) = xi(c.entry,j);
        xi(sb,j) = (inj / (eye(numel(sb)) - Pbb))';
    end
end
end
