function [Q,U,R,T,A,W] = getAvgHandles(self)
% [Q,U,R,T,A,W] = GETAVGHANDLES()
% Get handles for mean performance metrics.
%
% Q(i,r): mean queue-length of class r at node i
% U(i,r): mean utilization of class r at node i
% R(i,r): mean response time of class r at node i (summed across visits)
% T(i,r): mean throughput of class r at node i
% A(i,r): mean arrival rate of class r at node i
% W(i,r): mean residence time of class r at node i
%
% For tardiness metrics, use getAvgTardHandles() and getAvgSysTardHandles()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% Q = self.getAvgQLenHandles;
% U = self.getAvgUtilHandles;
% R = self.getAvgRespTHandles;
% T = self.getAvgTputHandles;
% A = self.getAvgArvRHandles;
% W = self.getAvgResidTHandles;

M = getNumberOfStations(self);
K = getNumberOfClasses(self);
classes = self.classes;
stations = self.stations;

isSource = false(M,1);
isSink = false(M,1);
hasServiceTunnel = false(M,1);
isServiceDefined = true(M,K);
for ist=1:M
    isSource(ist) = isa(stations{ist},'Source');
    isSink(ist) = isa(stations{ist},'Sink');
    hasServiceTunnel(ist) = strcmpi(class(stations{ist}.server),'ServiceTunnel');
    if ~hasServiceTunnel(ist)
        for r=1:K
            if isempty(stations{ist}.server.serviceProcess{r}) || stations{ist}.server.serviceProcess{r}{end}.isDisabled()
                isServiceDefined(ist,r) = false;
            end
        end
    end
end

% Identify cache hit/miss classes (these should not have their throughput disabled)
isCacheClass = false(1,K);
nodes = self.nodes;
for i=1:length(nodes)
    if isa(nodes{i},'Cache')
        hitClasses = nodes{i}.server.hitClass;
        missClasses = nodes{i}.server.missClass;
        for r=1:length(hitClasses)
            if hitClasses(r) > 0 && hitClasses(r) <= K
                isCacheClass(hitClasses(r)) = true;
            end
        end
        for r=1:length(missClasses)
            if missClasses(r) > 0 && missClasses(r) <= K
                isCacheClass(missClasses(r)) = true;
            end
        end
    end
end

if isempty(self.handles) || ~isfield(self.handles,'Q')
    Q = cell(M,K); % queue-length
    for ist=1:M
        for r=1:K
            Qir = Metric(MetricType.QLen, classes{r}, stations{ist});
            if isSource(ist)
                Qir.disabled = true;
            end
            if isSink(ist)
                Qir.disabled = true;
            end
            if ~hasServiceTunnel(ist)
                if ~isServiceDefined(ist,r)
                    Qir.disabled = true;
                end
            end
            Q{ist,r} = Qir;
        end
    end
    self.handles.Q = Q;
else
    Q = self.handles.Q;
end

if isempty(self.handles) || ~isfield(self.handles,'U')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    U = cell(M,1); % utilizations
    for ist=1:M
        for r=1:K
            Uir = Metric(MetricType.Util, classes{r}, stations{ist});
            if isSource(ist)
                Uir.disabled = true;
            end
            if isSink(ist)
                Uir.disabled = true;
            end
            if isa(stations{ist},'Join') || isa(stations{ist},'Fork')
                Uir.disabled = true;
            end
            if ~hasServiceTunnel(ist)
                if ~isServiceDefined(ist,r)
                    Uir.disabled = true;
                end
            end
            U{ist,r} = Uir;
        end
    end
    self.handles.U = U;
else
    U = self.handles.U;
end

if isempty(self.handles) || ~isfield(self.handles,'R')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    R = cell(M,K); % response times
    for ist=1:M
        for r=1:K
            Rir = Metric(MetricType.RespT, classes{r}, stations{ist});
            if isSource(ist)
                Rir.disabled = true;
            end
            if isSink(ist)
                Rir.disabled = true;
            end
            if ~hasServiceTunnel(ist)
                if ~isServiceDefined(ist,r)
                    Rir.disabled = true;
                end
            end
            R{ist,r} = Rir;
        end
    end
    self.handles.R = R;
else
    R = self.handles.R;
end

if isempty(self.handles) || ~isfield(self.handles,'W')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);
    
    W = cell(M,K); % response times
    for ist=1:M
        for r=1:K
            Wir = Metric(MetricType.ResidT, classes{r}, stations{ist});
            if isSource(ist)
                Wir.disabled = true;
            end
            if isSink(ist)
                Wir.disabled = true;
            end
            if ~hasServiceTunnel(ist)
                if ~isServiceDefined(ist,r)
                    Wir.disabled = true;
                end
            end
            W{ist,r} = Wir;
        end
    end
    self.handles.W = W;
else
    W = self.handles.W;
end


if isempty(self.handles) || ~isfield(self.handles,'T')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);

    T = cell(1,K); % throughputs
    for ist=1:M
        for r=1:K
            Tir = Metric(MetricType.Tput, classes{r}, stations{ist});
            if ~hasServiceTunnel(ist)
                % Don't disable throughput for cache hit/miss classes
                if ~isServiceDefined(ist,r) && ~isCacheClass(r)
                    Tir.disabled = true;
                end
            end
            T{ist,r} = Tir;
        end
    end
    self.handles.T = T;
else
    T = self.handles.T;
end

if isempty(self.handles) || ~isfield(self.handles,'A')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);

    A = cell(1,K); % arrival rate
    for ist=1:M
        for r=1:K
            Air = Metric(MetricType.ArvR, classes{r}, stations{ist});
            if ~hasServiceTunnel(ist)
                % Don't disable arrival rate for cache hit/miss classes
                if ~isServiceDefined(ist,r) && ~isCacheClass(r)
                    Air.disabled = true;
                end
            end
            A{ist,r} = Air;
        end
    end
    self.handles.A = A;
else
    A = self.handles.A;
end

if isempty(self.handles) || ~isfield(self.handles,'Tard')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);

    Tard = cell(M,K); % tardiness
    for ist=1:M
        for r=1:K
            Tardir = Metric(MetricType.Tard, classes{r}, stations{ist});
            if isSource(ist)
                Tardir.disabled = true;
            end
            if isSink(ist)
                Tardir.disabled = true;
            end
            if ~hasServiceTunnel(ist)
                if ~isServiceDefined(ist,r)
                    Tardir.disabled = true;
                end
            end
            Tard{ist,r} = Tardir;
        end
    end
    self.handles.Tard = Tard;
else
    Tard = self.handles.Tard;
end

if isempty(self.handles) || ~isfield(self.handles,'SysTard')
    M = getNumberOfStations(self);
    K = getNumberOfClasses(self);

    SysTard = cell(1,K); % system tardiness
    for r=1:K
        SysTardr = Metric(MetricType.SysTard, classes{r}, []);
        SysTard{1,r} = SysTardr;
    end
    self.handles.SysTard = SysTard;
else
    SysTard = self.handles.SysTard;
end

end