%{ @file qbd_setupdelayoff.m
 %  @brief Queue-length analysis for system with setup delay and turn-off
 %
 %  @author LINE Development Team
%}

%{
 % @brief Analyzes queue length for system with setup and turn-off phases
 %
 % @details
 % This function performs queue-length analysis for a queueing system with
 % setup delay and turn-off periods using QBD methods.
 %
 % @par Syntax:
 % @code
 % QN = qbd_setupdelayoff(lambda, mu, alpharate, alphascv, betarate, betascv)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rate
 % <tr><td>mu<td>Service rate
 % <tr><td>alpharate<td>Rate of setup delay phase
 % <tr><td>alphascv<td>Squared coefficient of variation for setup delay
 % <tr><td>betarate<td>Rate of turn-off phase
 % <tr><td>betascv<td>Squared coefficient of variation for turn-off period
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>QN<td>Average queue length
 % </table>
%}
function QN = qbd_setupdelayoff(lambda, mu, alpharate, alphascv, betarate, betascv)

% The phases are given as rates, so an exponential one is built directly from
% its rate rather than round-tripped through its mean. APH.fitMeanAndSCV maps
% any mean at or below GlobalConstants.FineTol to Exp(Inf), and an Immediate
% setup or delay-off has rate GlobalConstants.Immediate whose mean is exactly
% FineTol; the round trip therefore turned a finite 1e8 rate into an infinite
% one, put -Inf in the generator, and collapsed QBD_pi onto the boundary level,
% so the level loop below ran past the end of pn. For SCV 1 the fit returns
% Exp(1/MEAN) anyway, so this is the same process everywhere else.
% Both phases are taken in canonical Coxian form, whose entry vector is [1 0 ..]
% for every SCV. The chain below overloads its phase indices (at the boundary
% level phase 1 is the off server and phases na+1.. are the delay-off; above it
% phases 1..na are the setup and na+1 is the busy server), so a level-up
% transition cannot also redistribute the phase: an arrival to an off server has
% to enter the setup at phase 1. APH.fitMeanAndSCV violates that for SCV > 1,
% where it returns a hyperexponential entered at phase 2 with probability 3/4,
% which the chain then silently entered at phase 1 regardless.
alpha = coxian_phase(alpharate, alphascv);
na = length(alpha{1});
% Completion rate of each phase. A Coxian may complete from any phase, not only
% from the last one, so these are read per phase rather than off the diagonal.
ta = -sum(alpha{1},2);

beta = coxian_phase(betarate, betascv);
nb = length(beta{1});
tb = -sum(beta{1},2);
n = na+nb;

F = zeros(n); % forward transitions
B = zeros(n);  % backward transitions
for i=1:na
    F(i,i) = lambda;
end
for i=1:nb
    F(na+i,na+1) = lambda;
end
F(na+1,na+1) = lambda;
B(na+1,na+1) = mu;

L = zeros(n);  % local transitions
for i=1:na
    % Whole generator row, so a phase that both advances and completes (any
    % Coxian with SCV > 1) contributes both; reading only the diagonal and the
    % strict upper triangle assumed a pure series, which holds only for SCV <= 1.
    L(i,1:na) = alpha{1}(i,1:na);
    L(i,i) = L(i,i) - lambda;
    L(i,na+1) = ta(i); % setup completes from phase i -> busy server
end
L(na+1,na+1) = -mu -lambda;
for i=2:nb
    L(na+i,na+i) = -lambda;
end
%[B,L,F]

L0 = zeros(na+nb);  % local transitions at the boundary level
for i=1:na
    L0(i,i) = -lambda;
end
for i=1:nb
    % As above, the whole generator row: the delay-off may expire from any
    % phase, and its phase-to-phase rate is the generator entry, not the
    % diagonal (they coincide only for a series phase).
    L0(na+i,na+(1:nb)) = beta{1}(i,1:nb);
    L0(na+i,na+i) = L0(na+i,na+i) - lambda;
    L0(na+i,1) = tb(i); % delay-off expires from phase i -> server off
end
%[L0,F]

[~,R,~] = QBD_CR(B,L,F);
pn = QBD_pi(B,L0,R);

% see _kb/03-api-layer.md (QBD boundary conventions) for rationale
QN = 0; % queue-lengths
j = n+1;
ni = 0;
while j+n-1 <= length(pn)
    ni = ni + 1;
    QN = QN + ni*sum(pn(j:(j+n-1)));
    j = j +n;
end
end

function p = coxian_phase(rate, scv)
% P = COXIAN_PHASE(RATE, SCV) - canonical PH representation of a phase given
% its rate, entered at phase 1 for every SCV. Exponential phases are built from
% the rate directly, which is exact and avoids the mean-based FineTol cutoff in
% the fitters.
if scv == 1.0
    p = Exp(rate).getProcess;
else
    p = Coxian.fitMeanAndSCV(1/rate, scv).getProcess;
end
end
