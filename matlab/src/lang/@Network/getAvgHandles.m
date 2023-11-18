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

% Copyright (c) 2012-2023, Imperial College London
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
for i=1:M
    isSource(i) = isa(stations{i},'Source');
    isSink(i) = isa(stations{i},'Sink');
    hasServiceTunnel(i) = strcmpi(class(stations{i}.server),'ServiceTunnel');
    if ~hasServiceTunnel(i)
        for r=1:K
            if isempty(stations{i}.server.serviceProcess{r}) || stations{i}.server.serviceProcess{r}{end}.isDisabled()
                isServiceDefined(i,r) = false;
            end
        end
    end
end

if isempty(self.handles) || ~isfield(self.handles,'Q')
    Q = cell(M,K); % queue-length
    for i=1:M
        for r=1:K
            Qir = Metric(MetricType.QLen, classes{r}, stations{i});
            if isSource(i)
                Qir.disabled = true;
            end
            if isSink(i)
                Qir.disabled = true;
            end
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Qir.disabled = true;
                end
            end
            Q{i,r} = Qir;
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
    for i=1:M
        for r=1:K
            Uir = Metric(MetricType.Util, classes{r}, stations{i});
            if isSource(i)
                Uir.disabled = true;
            end
            if isSink(i)
                Uir.disabled = true;
            end
            if isa(stations{i},'Join') || isa(stations{i},'Fork')
                Uir.disabled = true;
            end
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Uir.disabled = true;
                end
            end
            U{i,r} = Uir;
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
    for i=1:M
        for r=1:K
            Rir = Metric(MetricType.RespT, classes{r}, stations{i});
            if isSource(i)
                Rir.disabled = true;
            end
            if isSink(i)
                Rir.disabled = true;
            end
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Rir.disabled = true;
                end
            end
            R{i,r} = Rir;
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
    for i=1:M
        for r=1:K
            Wir = Metric(MetricType.ResidT, classes{r}, stations{i});
            if isSource(i)
                Wir.disabled = true;
            end
            if isSink(i)
                Wir.disabled = true;
            end
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Wir.disabled = true;
                end
            end
            W{i,r} = Wir;
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
    for i=1:M
        for r=1:K
            Tir = Metric(MetricType.Tput, classes{r}, stations{i});
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Tir.disabled = true;
                end
            end
            T{i,r} = Tir;
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
    for i=1:M
        for r=1:K
            Air = Metric(MetricType.ArvR, classes{r}, stations{i});
            if ~hasServiceTunnel(i)
                if ~isServiceDefined(i,r)
                    Air.disabled = true;
                end
            end
            A{i,r} = Air;
        end
    end
    self.handles.A = A;
else
    A = self.handles.A;
end

end