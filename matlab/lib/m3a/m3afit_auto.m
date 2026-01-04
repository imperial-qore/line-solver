function [MMAP]=m3afit_auto(mtrace, varargin)
% [MMAP] = m3afit_auto(T,'option1',val1,'option2',val2,...)
%
% DESCRIPTION
% Automatic fitting of trace T into a Marked Markovian Arrival Process
% based on the M3A toolbox.
%
% INPUT
% T        - data structure returned by m3afit_init
%
% OUTPUT
% MMAP          - fitted marked MAP
%
% EXAMPLE
%  T = m3afit_init(S,C)
%  MMAP = m3afit_auto(T,'NumStates',2)
%
% OPTION LIST
% 'Method'  - 0 for inter-arrival fitting, 1 for counting process fitting.
% 'NumStates'  - integer, mandatory parameter
%
% REFERENCES
% [1] A. Sansottera, G. Casale, P. Cremonesi. Fitting Second-Order Acyclic
%     Marked Markovian Arrival Processes. IEEE/IFIP DSN 2013.
% [2] G. Casale, A. Sansottera, P. Cremonesi. Compact Markov-Modulated
%     Models for Multiclass Trace Fitting. European Journal of Operations
%     Research, 2016.
%

%% options
OptionNames = [
    'Method    ';
    'NumStates ';
    ];

OptionTypes = [
    'numeric';
    'numeric'];

OptionValues = [];
for i = 1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end

% Default settings
options.NumStates = [];
options.Method = 1;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
%% fitting algorithm run parameterization

if options.Method==1 % counts
    options.Timescale = 10*mean(mtrace.S);
    options.TimescaleAsy = max(10*options.Timescale, (sum(mtrace.S)-mtrace.S(1)) / 100);
end

%% fitting
if mtrace.NumClasses == 2 && options.NumStates == 2 && options.Method==0 % inter-arrivals
    MMAPType=sprintf('%d-state MAMAP[%d]',options.NumStates,mtrace.NumClasses);
    fprintf(1,'Init: M3A algorithm will search for a %s\n',MMAPType);
    MMAP = mamap22_fit_gamma_fs_trace(mtrace.S,mtrace.C);
elseif  mtrace.NumClasses > 2 && options.NumStates >= 2 && options.Method==0 % inter-arrivals
    MMAPType=sprintf('%d-state MAMAP[%d]',options.NumStates,mtrace.NumClasses);
    fprintf(1,'Init: M3A algorithm will search for a %s\n',MMAPType);
    MMAP = mamap2m_fit_trace(mtrace.S,mtrace.C);
elseif  mtrace.NumClasses >= 2 && options.NumStates == 2 && options.Method==1 % counts
    MMAPType=sprintf('%d-state M3PP[%d]',options.NumStates,mtrace.NumClasses);
    fprintf(1,'Init: M3A algorithm will search for a %s\n',MMAPType);
    MMAP = m3pp2m_fitc_trace(mtrace.S,mtrace.C,'approx_ag',options.Timescale,options.TimescaleAsy);
elseif  mtrace.NumClasses >= 2 && options.NumStates > 2 && options.Method==1 % counts
    MMAP = m3pp_superpos_fitc_trace(mtrace.S,mtrace.C,options.Timescale,options.TimescaleAsy); % superposition
%    MMAP = m3pp_interleave_fitc_trace(T,A,options.Timescale,options.TimescaleAsy,[1 0 0; 0 1 0; 1 0 0; 0 0 1]);
else
    fprintf(1,'Output: M3A algorithm could *not* obtain a valid MMAP.\n');
    MMAP={};
    return
end

MMAPResType=sprintf('%d-state M3PP[%d]',length(MMAP{1}),size(MMAP,2)-2);
if mmap_isfeasible(MMAP)
    fprintf(1,'Output: M3A algorithm found a valid %s.\n',MMAPResType);
else
    fprintf(1,'Output: M3A algorithm could *not* obtain a valid MMAP.\n');
end

%     elseif strcmp(model, 'm3pp22-interleave-cov')
%         fprintf('Fitting model\n');
%         fit = m3pp_interleave_fitc_trace(T,A,100,1000,[1 0 0; 0 1 0; 1 0 0; 0 0 1]);

end