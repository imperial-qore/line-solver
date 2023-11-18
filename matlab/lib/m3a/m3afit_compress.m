function [MMAP]=m3afit_compress(MMAP, varargin)
% [MMAP] = m3afit_compress(MMAP,'option1',val1,'option2',val2,...)
%
% DESCRIPTION
% Compresses the representation of a Marked Markovian Arrival Process
% based on the M3A toolbox.
%
% INPUT
% MMAP         - a feasible MMAP={D0,D1,D11,...,D1c}
%
% OUTPUT
% MMAP          - compressed MMAP
%
% EXAMPLE
%  D0=rand(3); D11 = rand(3); D12=rand(3); D1=D11+D12;
%  MMAP = {D0,D1,D11,D12}; 
%  MMAP = mmap_normalize(MMAP);
%  MMAP_compressed = m3afit_compress(MMAP,'Method',0) % two-state AMAP compression
%
% OPTION LIST
% 'Method'  - 0, 2-state acyclic MAP (AMAP) compression
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
options.NumStates = 2;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
options.NumClasses = length(MMAP)-2; 
%% fitting algorithm run parameterization
%N/A

%% fitting
if options.Method==0 % two-state AMAP compression
    MMAPType=sprintf('%d-state AMAP[%d]',options.NumStates,options.NumClasses);
    fprintf(1,'Init: M3A will search for a %s\n',MMAPType);
    MMAP = mamap2m_fit_gamma_fb_mmap(MMAP);
else
    fprintf(1,'Unsupported option.\n');
    MMAP={};
    return
end

MMAPResType=sprintf('%d-state M3PP[%d]',length(MMAP{1}),size(MMAP,2)-2);
if mmap_isfeasible(MMAP)
    fprintf(1,'Output: M3A found a valid %s.\n',MMAPResType);
else
    fprintf(1,'Output: M3A could *not* obtain a valid MMAP.\n');
end

end