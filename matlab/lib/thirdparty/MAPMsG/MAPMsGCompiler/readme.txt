% O. Gursoy, K. A. Mehr, N. Akar:
% The MAP/M/s + G Call Center Model with Generally Distributed Patience Times: Steady-state Solution and First Passage Time Distribution
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Steady-state and first passage time solver for MAP/M/s + G queues
%
% INPUTS
%-----------------------------
% -SOLUTION: Solution type(steady state or first passage actual/virtual)
% -SERVERSIZE: server size of the call center problem
% -mu
% -C: C of MAP(C,D) 
% -D: D of MAP(C,D) 
% -MAPSIZE: Order of MAP(C,D) distribution, Service Rate(mu)
% -QUANTIZATION: number of regimes that will quantize abandonment function ga
% -BoundaryLevels: Boundary levels of Regimes
% -ga: Abandonment Probability in Each Regime (starts with a 0 value)
% -gabanfunction: for continious abandonment functions enter ccdf of abandonement as a function of x
%
%
% INPUTS for first passage times(these are not necessary for steady state
% calculations)
%-----------------------------
% -b: threshold level 
% -tau: time horizon
% -CMEorErlang: indicator of which distribution is chosen for approximation
% -OrderofPHCME: order of the chosen approximator (25, 51 or 101)
% -pi0: initial server occupancy distribution
% -theta0: initial MAP state distribution
%
% OUTPUTS
%---------------------------
% -steadyStateResult: 1x7 steady state result vector including PR{W=0}
% PR{W=0|S}   PR{A}    E{W|S}   Var{W|S} Fw|s,w>0(.1) Fw|s,w>0(.2)
% respectively(when SOLUTION=1).
%
% -FirstPassageTimeofVirtual: First passage time probability for the virtual waiting time(when SOLUTION=2).
%
% -FirstPassageTimeofActual: First passage time probability for the actual
% waiting time(when SOLUTION=3).
%
%
% IMPORTANT:
%-----------
% For first passage time cases the initial wait time, a=0.
%
% 
%_________________________________________________________________
%/////////////////////////////////////////////////////////////////
