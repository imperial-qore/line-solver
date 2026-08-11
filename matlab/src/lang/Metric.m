function met = Metric(type, class, station)
% An output metric of a Solver, such as a performance index
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

persistent emptyMetric;

if nargin > 2
    if isempty(emptyMetric)
        met = struct('type',type,'class',class,'station',station,'disabled',false,'transient',false);
        emptyMetric = met;
    else
        met = emptyMetric; % copying faster than creating a new struct
        met.class = class;
        met.type = type;
        met.station = station;
    end
else
    if isempty(emptyMetric)
        met = struct('type',type,'class',class,'station',NaN,'disabled',false,'transient',false);
        emptyMetric = met;
    else
        met = emptyMetric;
        met.class = class;
        met.type = type;
        met.station = NaN;
    end
end % nargin
end