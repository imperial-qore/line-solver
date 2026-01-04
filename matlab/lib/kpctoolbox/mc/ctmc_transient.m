%{ @file ctmc_transient.m
 %  @brief Transient analysis of a continuous-time Markov chain using ODE solvers
 %
 %  @author LINE Development Team
%}

%{
 % @brief Transient analysis of a continuous-time Markov chain using ODE solvers
 %
 % @details
 % Computes the transient state probabilities using ODE solvers.
 %
 % @par Syntax:
 % @code
 % [pi, t] = ctmc_transient(Q, pi0)
 % [pi, t] = ctmc_transient(Q, pi0, t0, t1)
 % [pi, t] = ctmc_transient(Q, pi0, t0, t1, useStiff, reltol, timestep)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Infinitesimal generator matrix
 % <tr><td>pi0<td>Initial probability distribution vector
 % <tr><td>t0<td>(Optional) Start time
 % <tr><td>t1<td>(Optional) End time
 % <tr><td>useStiff<td>(Optional) Boolean to use stiff ODE solver (ode15s) instead of ode23. Default: false
 % <tr><td>reltol<td>(Optional) Relative tolerance for the ODE solver. Default: 1e-3
 % <tr><td>timestep<td>(Optional) Fixed time step for the output
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>pi<td>State probability vectors at each time point
 % <tr><td>t<td>Time points corresponding to the probabilities
 % </table>
%}
function [pi,t]=ctmc_transient(Q,pi0,t0,t1,useStiff,reltol,timestep)
if ~exist('useStiff','var')
    useStiff = false;
end
if ~exist('reltol','var')
    reltol = 1e-3;
end
if ~exist('timestep','var')
    timestep = [];
end
if nargin==2
    t1=pi0;
    t0=0;
    pi0=ones(1,length(Q));pi0=pi0/sum(pi0);
end
if nargin==3
    t1=t0;
    t0=0;
end

odeoptions = odeset('RelTol', reltol);

% Use fixed intervals if timestep parameter is specified
if ~isempty(timestep) && timestep > 0
    t = t0:timestep:t1;
    if t(end) ~= t1
        t = [t, t1]; % Ensure t1 is included
    end
    if useStiff
        [t,pi]=ode15s(@ctmc_transientode, t, pi0, odeoptions);
    else
        [t,pi]=ode23(@ctmc_transientode, t, pi0, odeoptions);
    end
else
    % Original adaptive stepping behavior
    if isinf(t1)
        nonZeroRates = abs(Q(Q~=0));
        nonZeroRates = nonZeroRates( nonZeroRates >0 );
        T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
        if useStiff
            [t,pi]=ode15s(@ctmc_transientode, [t0,T], pi0, odeoptions);
        else
            [t,pi]=ode23(@ctmc_transientode, [t0,T], pi0, odeoptions);
        end
    else
        if useStiff
            [t,pi]=ode15s(@ctmc_transientode, [t0,t1], pi0, odeoptions);
        else
            [t,pi]=ode23(@ctmc_transientode, [t0,t1], pi0, odeoptions);
        end
    end
end
%[t,pi]=ode113(@ctmc_transientode,[t0,t1],pi0);

function dpidt=ctmc_transientode(t,pi)
pi=pi(:)';
dpidt=pi*Q;
dpidt=dpidt(:);
end

end