function QN2MATLAB(model, modelName, fid)
% QN2MATLAB(MODEL, MODELNAME, FID)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%global GlobalConstants.CoarseTol; % Tolerance for distribution fitting

if nargin<2%~exist('modelName','var')
    modelName='myModel';
end
if nargin<3%~exist('fid','var')
    fid=1;
end
if ischar(fid)
    fid = fopen(fid,'w');
end
sn = model.getStruct;
%% initialization
fprintf(fid,'model = Network(''%s'');\n',modelName);
fprintf(fid,'\n%%%% Block 1: nodes');
rt = sn.rt;
rtnodes = sn.rtnodes;
hasSink = 0;
sourceID = 0;
PH = sn.proc;

fprintf(fid,'\n');
%% write nodes
for i= 1:sn.nnodes
    switch sn.nodetype(i)
        case NodeType.Source
            sourceID = i;
            fprintf(fid,'node{%d} = Source(model, ''%s'');\n',i,sn.nodenames{i});
            hasSink = 1;
        case NodeType.Delay
            fprintf(fid,'node{%d} = DelayStation(model, ''%s'');\n',i,sn.nodenames{i});
        case NodeType.Queue
            fprintf(fid,'node{%d} = Queue(model, ''%s'', SchedStrategy.%s);\n', i, sn.nodenames{i}, SchedStrategy.toProperty(sn.sched(sn.nodeToStation(i))));
            if sn.nservers(sn.nodeToStation(i))>1
                fprintf(fid,'node{%d}.setNumServers(%d);\n', i, sn.nservers(sn.nodeToStation(i)));
            end
        case NodeType.Router
            fprintf(fid,'node{%d} = Router(model, ''%s'');\n',i,sn.nodenames{i});
        case NodeType.Fork
            fprintf(fid,'node{%d} = Fork(model, ''%s'');\n',i,sn.nodenames{i});
        case NodeType.Join
            fprintf(fid,'node{%d} = Join(model, ''%s'', node{%d});\n',i,sn.nodenames{i},find(sn.fj(:,i)));
        case NodeType.Sink
            fprintf(fid,'node{%d} = Sink(model, ''%s'');\n',i,sn.nodenames{i});
        case NodeType.ClassSwitch
            %            csMatrix = eye(sn.nclasses);
            %            fprintf(fid,'\ncsMatrix%d = zeros(%d);\n',i,sn.nclasses);
            %             for k = 1:sn.nclasses
            %                 for c = 1:sn.nclasses
            %                     for m=1:sn.nnodes
            %                         % routing matrix for each class
            %                         csMatrix(k,c) = csMatrix(k,c) + rtnodes((i-1)*sn.nclasses+k,(m-1)*sn.nclasses+c);
            %                     end
            %                 end
            %             end
            %            for k = 1:sn.nclasses
            %                for c = 1:sn.nclasses
            %                    if csMatrix(k,c)>0
            %                        fprintf(fid,'csMatrix%d(%d,%d) = %f; %% %s -> %s\n',i,k,c,csMatrix(k,c),sn.classnames{k},sn.classnames{c});
            %                    end
            %                end
            %            end
            %fprintf(fid,'node{%d} = ClassSwitch(model, ''%s'', csMatrix%d);\n',i,sn.nodenames{i},i);
            fprintf(fid,'node{%d} = Router(model, ''%s''); %% Class switching is embedded in the routing matrix \n',i,sn.nodenames{i});
    end
end
%% write classes
fprintf(fid,'\n%%%% Block 2: classes\n');
for k = 1:sn.nclasses
    if sn.njobs(k)>0
        if isinf(sn.njobs(k))
            fprintf(fid,'jobclass{%d} = OpenClass(model, ''%s'', %d);\n',k,sn.classnames{k},sn.classprio(k));
        else
            fprintf(fid,'jobclass{%d} = ClosedClass(model, ''%s'', %d, node{%d}, %d);\n',k,sn.classnames{k},sn.njobs(k),sn.stationToNode(sn.refstat(k)),sn.classprio(k));
        end
    else
        % if the reference node is unspecified, as in artificial classes,
        % set it to the first node where the rate for this class is
        % non-null
        iref = 0;
        for i=1:sn.nstations
            if sum(nnz(sn.proc{i}{k}{1}))>0
                iref = i;
                break
            end
        end
        if isinf(sn.njobs(k))
            fprintf(fid,'jobclass{%d} = OpenClass(model, ''%s'', %d);\n',k,sn.classnames{k},sn.classprio(k));
        else
            fprintf(fid,'jobclass{%d} = ClosedClass(model, ''%s'', %d, node{%d}, %d);\n',k,sn.classnames{k},sn.njobs(k),iref,sn.classprio(k));
        end
    end
end
fprintf(fid,'\n');
%% arrival and service processes
for i=1:sn.nstations
    for k=1:sn.nclasses
        if sn.nodetype(sn.stationToNode(i)) ~= NodeType.Join
            if isprop(model.stations{i},'serviceProcess') && strcmp(class(model.stations{i}.serviceProcess{k}),'Replayer')
                switch sn.sched(i)
                    case SchedStrategy.EXT
                        fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Replayer(''%s'')); %% (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                    otherwise
                        fprintf(fid,'node{%d}.setService(jobclass{%d}, Replayer(''%s'')); %% (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                end
            elseif isprop(model.stations{i},'serviceProcess') && strcmp(class(model.stations{i}.serviceProcess{k}),'Trace')
                switch sn.sched(i)
                    case SchedStrategy.EXT
                        fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Trace(''%s'')); %% (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                    otherwise
                        fprintf(fid,'node{%d}.setService(jobclass{%d}, Trace(''%s'')); %% (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                end
            else
                SCVik = map_scv(PH{i}{k});
                if SCVik >= 0.5
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            if SCVik == 1
                                meanik = map_mean(PH{i}{k});
                                if meanik < GlobalConstants.CoarseTol
                                    fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Immediate()); %% (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k}); 
                                else
                                    fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Exp.fitMean(%f)); %% (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                end
                            else
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, APH.fitMeanAndSCV(%f,%f)); %% (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),SCVik,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                        otherwise
                            if SCVik == 1
                                meanik = map_mean(PH{i}{k});
                                if meanik < GlobalConstants.CoarseTol
                                    fprintf(fid,'node{%d}.setService(jobclass{%d}, Immediate()); %% (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k}); 
                                else
                                    fprintf(fid,'node{%d}.setService(jobclass{%d}, Exp.fitMean(%f)); %% (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                end
                            else
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, APH.fitMeanAndSCV(%f,%f)); %% (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),SCVik,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                    end
                else
                    % this could be made more precised by fitting into a 2-state
                    % APH, especially if SCV in [0.5,0.1]
                    nPhases = max(1,round(1/SCVik));
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            if isnan(PH{i}{k}{1})
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Disabled.getInstance()); %% (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setArrival(jobclass{%d}, Erlang(%f,%f)); %% (%s,%s)\n',sn.stationToNode(i),k,nPhases/map_mean(PH{i}{k}),nPhases,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                        otherwise
                            if isnan(PH{i}{k}{1})
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Disabled.getInstance()); %% (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            else
                                fprintf(fid,'node{%d}.setService(jobclass{%d}, Erlang(%f,%f)); %% (%s,%s)\n',sn.stationToNode(i),k,nPhases/map_mean(PH{i}{k}),nPhases,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                    end
                end
            end
        end
    end
end

fprintf(fid,'\n%%%% Block 3: topology');
if hasSink
    rt(sn.nstations*sn.nclasses+(1:sn.nclasses),sn.nstations*sn.nclasses+(1:sn.nclasses)) = zeros(sn.nclasses);
    for k=find(isinf(sn.njobs))' % for all open classes
        for i=1:sn.nstations
            % all open class transitions to ext station are re-routed to sink
            rt((i-1)*sn.nclasses+k, sn.nstations*sn.nclasses+k) = rt((i-1)*sn.nclasses+k, (sourceID-1)*sn.nclasses+k);
            rt((i-1)*sn.nclasses+k, (sourceID-1)*sn.nclasses+k) = 0;
        end
    end
end

fprintf(fid,'\n');
fprintf(fid,'P = model.initRoutingMatrix(); %% initialize routing matrix \n',sn.nclasses,sn.nclasses,sn.nnodes,sn.nnodes);
for k = 1:sn.nclasses
    for c = 1:sn.nclasses
        for i=1:sn.nnodes
            for m=1:sn.nnodes
                % routing matrix for each class
                myP{k,c}(i,m) = rtnodes((i-1)*sn.nclasses+k,(m-1)*sn.nclasses+c);
                if myP{k,c}(i,m) > 0 && sn.nodetype(i) ~= NodeType.Sink
                    % do not change %d into %f to avoid round-off errors in
                    % the total probability
                    if sn.nodetype(i)==NodeType.Fork
                        fprintf(fid,'P{%d,%d}(%d,%d) = 1.0; %% (%s,%s) -> (%s,%s)\n',k,c,i,m,sn.nodenames{i},sn.classnames{k},sn.nodenames{m},sn.classnames{c});
                    else
                        fprintf(fid,'P{%d,%d}(%d,%d) = %d; %% (%s,%s) -> (%s,%s)\n',k,c,i,m,myP{k,c}(i,m),sn.nodenames{i},sn.classnames{k},sn.nodenames{m},sn.classnames{c});
                    end
                end
            end
        end
        %fprintf(fid,'P{%d,%d} = %s;\n',k,c,mat2str(myP{k,c}));
    end
end

fprintf(fid,'model.link(P);\n');
%if fid~=1
%    fclose(fid);
%end
end
