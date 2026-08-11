% @brief MATLAB Coder script to generate MEX functions for rl_ module.
%
% This script generates MEX (MATLAB Executable) versions of reinforcement-
% learning helper functions for improved performance.
%
% Skipped functions (Coder-incompatible):
%   rl_env, rl_env_general          - classdef handle classes referencing
%                                     a Network model and SolverSSA; Coder
%                                     does not support handle classes that
%                                     close over OOP model objects.
%   rl_td_agent, rl_td_agent_general - classdef handle classes that
%                                     instantiate solvers and call State,
%                                     EventType, NodeType, etc.
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% No codegen-compatible functions in this domain.
