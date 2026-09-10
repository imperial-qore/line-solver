function QN2JAVA(model, modelName, fid, headers)
% QN2JAVA(MODEL, MODELNAME, FID, HEADERS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%global GlobalConstants.CoarseTol; % Tolerance for distribution fitting

if nargin<2 %~exist('modelName','var')
    modelName='myModel';
end
if nargin<3 %~exist('fid','var')
    fid=1;
end
if nargin<4 %~exist('fid','var')
    headers=true;
end
if ischar(fid)
    fid = fopen(fid,'w');
end
sn = model.getStruct;
%% initialization
if headers
    fprintf(fid,'\tpublic static Network ex() {\n',modelName);
end
fprintf(fid,'\t\tNetwork model = new Network("%s");\n',modelName);
fprintf(fid,'\n\t\t// Block 1: nodes');
rt = sn.rt;
rtnodes = sn.rtnodes;
hasSink = 0;
sourceID = 0;
PH = sn.proc;

fprintf(fid,'\t\t\t\n');
%% write nodes
for i= 1:sn.nnodes
    switch sn.nodetype(i)
        case NodeType.Source
            sourceID = i;
            fprintf(fid,'\t\tSource node%d = new Source(model, "%s");\n',i,sn.nodenames{i});
            hasSink = 1;
        case NodeType.Delay
            fprintf(fid,'\t\tDelay node%d = new Delay(model, "%s");\n',i,sn.nodenames{i});
        case NodeType.Queue
            fprintf(fid,'\t\tQueue node%d = new Queue(model, "%s", SchedStrategy.%s);\n', i, sn.nodenames{i}, SchedStrategy.toProperty(SchedStrategy.toText(sn.sched(sn.nodeToStation(i)))));
            if sn.nservers(sn.nodeToStation(i))>1
                if isinf(sn.nservers(sn.nodeToStation(i)))
                    fprintf(fid,'\t\tnode%d.setNumberOfServers(Integer.MAX_VALUE);\n', i);
                else
                    fprintf(fid,'\t\tnode%d.setNumberOfServers(%d);\n', i, sn.nservers(sn.nodeToStation(i)));
                end
            end
        case NodeType.Router
            fprintf(fid,'\t\tRouter node%d = new Router(model, "%s");\n',i,sn.nodenames{i});
        case NodeType.Fork
            fprintf(fid,'\t\tFork node%d = new Fork(model, "%s");\n',i,sn.nodenames{i});
        case NodeType.Join
            fprintf(fid,'\t\tJoin node%d = new Join(model, "%s", node%d);\n',i,sn.nodenames{i},find(sn.fj(:,i)));
        case NodeType.Sink
            fprintf(fid,'\t\tSink node%d = new Sink(model, "%s");\n',i,sn.nodenames{i});
        case NodeType.ClassSwitch
            %            csMatrix = eye(sn.nclasses);
            %            fprintf(fid,'\t\t\ncsMatrix%d = zeros(%d);\n',i,sn.nclasses);
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
            %                        fprintf(fid,'\t\tcsMatrix%d(%d,%d) = %f; %% %s -> %s\n',i,k,c,csMatrix(k,c),sn.classnames{k},sn.classnames{c});
            %                    end
            %                end
            %            end
            %fprintf(fid,'\t\tnode%d = ClassSwitch(model, "%s", csMatrix%d);\n',i,sn.nodenames{i},i);
            fprintf(fid,'\t\tRouter node%d = new Router(model, "%s"); // Dummy node, class switching is embedded in the routing matrix P \n',i,sn.nodenames{i});
    end
end
%% write classes
fprintf(fid,'\n\t\t// Block 2: classes\n');
for k = 1:sn.nclasses
    if sn.njobs(k)>0
        if isinf(sn.njobs(k))
            fprintf(fid,'\t\tOpenClass jobclass%d = new OpenClass(model, "%s", %d);\n',k,sn.classnames{k},sn.classprio(k));
        else
            fprintf(fid,'\t\tClosedClass jobclass%d = new ClosedClass(model, "%s", %d, node%d, %d);\n',k,sn.classnames{k},sn.njobs(k),sn.stationToNode(sn.refstat(k)),sn.classprio(k));
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
            fprintf(fid,'\t\tOpenClass jobclass%d = new OpenClass(model, "%s", %d);\n',k,sn.classnames{k},sn.classprio(k));
        else
            fprintf(fid,'\t\tClosedClass jobclass%d = new ClosedClass(model, "%s", %d, node%d, %d);\n',k,sn.classnames{k},sn.njobs(k),iref,sn.classprio(k));
        end
    end
end
fprintf(fid,'\t\t\n');
%% arrival and service processes
for i=1:sn.nstations
    for k=1:sn.nclasses
        if sn.nodetype(sn.stationToNode(i)) ~= NodeType.Join
            if isprop(model.stations{i},'serviceProcess') && strcmp(class(model.stations{i}.serviceProcess{k}),'Replayer')
                switch sn.sched(i)
                    case SchedStrategy.EXT
                        fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Replayer("%s")); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                    otherwise
                        if sn.schedparam(i,k) ~= 1
                            fprintf(fid,'\t\tnode%d.setService(jobclass%d, Replayer("%s"), %f); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.schedparam(i,k),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                        else
                            fprintf(fid,'\t\tnode%d.setService(jobclass%d, Replayer("%s")); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                        end
                end
            else if isprop(model.stations{i},'serviceProcess') && strcmp(class(model.stations{i}.serviceProcess{k}),'Trace')
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Trace("%s")); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                        otherwise
                            if sn.schedparam(i,k) ~= 1
                                fprintf(fid,'\t\tnode%d.setService(jobclass%d, Trace("%s"), %f); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.schedparam(i,k),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            else
                                fprintf(fid,'\t\tnode%d.setService(jobclass%d, Trace("%s")); // (%s,%s)\n',sn.stationToNode(i),k,model.stations{i}.serviceProcess{k}.params{1}.paramValue,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                    end
            else
                SCVik = map_scv(PH{i}{k});
                if SCVik >= 0.5
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            if SCVik == 1
                                meanik = map_mean(PH{i}{k});
                                if meanik < GlobalConstants.CoarseTol
                                    fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Immediate.getInstance()); // (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                else
                                    fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Exp.fitMean(%f)); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                end
                            else
                                fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, APH.fitMeanAndSCV(%f,%f)); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),SCVik,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                        otherwise
                            if SCVik == 1
                                meanik = map_mean(PH{i}{k});
                                if meanik < GlobalConstants.CoarseTol
                                    fprintf(fid,'\t\tnode%d.setService(jobclass%d, Immediate.getInstance()); // (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                else
                                    if sn.schedparam(i,k) ~= 1
                                        fprintf(fid,'\t\tnode%d.setService(jobclass%d, Exp.fitMean(%f), %f); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),sn.schedparam(i,k),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                    else
                                        fprintf(fid,'\t\tnode%d.setService(jobclass%d, Exp.fitMean(%f)); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                    end
                                end
                            else
                                if sn.schedparam(i,k) ~= 1
                                    fprintf(fid,'\t\tnode%d.setService(jobclass%d, APH.fitMeanAndSCV(%f,%f), %f); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),SCVik,sn.schedparam(i,k),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                else
                                    fprintf(fid,'\t\tnode%d.setService(jobclass%d, APH.fitMeanAndSCV(%f,%f)); // (%s,%s)\n',sn.stationToNode(i),k,map_mean(PH{i}{k}),SCVik,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                end
                            end
                    end
                else
                    nPhases = max(1,round(1/SCVik));
                    switch sn.sched(i)
                        case SchedStrategy.EXT
                            if isnan(PH{i}{k}{1})
                                fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Disabled.getInstance()); // (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            else
                                fprintf(fid,'\t\tnode%d.setArrival(jobclass%d, Erlang(%f,%f)); // (%s,%s)\n',sn.stationToNode(i),k,nPhases/map_mean(PH{i}{k}),nPhases,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            end
                        otherwise
                            if isnan(PH{i}{k}{1})
                                fprintf(fid,'\t\tnode%d.setService(jobclass%d, Disabled.getInstance()); // (%s,%s)\n',sn.stationToNode(i),k,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                            else
                                if sn.schedparam(i,k) ~= 1
                                    fprintf(fid,'\t\tnode%d.setService(jobclass%d, Erlang(%f,%f), %f); // (%s,%s)\n',sn.stationToNode(i),k,nPhases/map_mean(PH{i}{k}),nPhases,sn.schedparam(i,k),sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                else
                                    fprintf(fid,'\t\tnode%d.setService(jobclass%d, Erlang(%f,%f)); // (%s,%s)\n',sn.stationToNode(i),k,nPhases/map_mean(PH{i}{k}),nPhases,sn.nodenames{sn.stationToNode(i)},sn.classnames{k});
                                end
                            end
                    end
                end
            end
            end
        end
    end
end
fprintf(fid,'\n\t\t// Block 3: topology');
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

fprintf(fid,'\t\n');
fprintf(fid,'\t\tRoutingMatrix routingMatrix = model.initRoutingMatrix(); \n');
% fprintf(fid,'\t\tRoutingMatrix routingMatrix = new RoutingMatrix(model,\n');
% fprintf(fid,'\t\t\t Arrays.asList(');
% for c = 1:sn.nclasses
%     fprintf(fid,'jobclass%d',c);
%     if c<sn.nclasses
%         fprintf(fid,', ');
%     else
%         fprintf(fid,'),\n');
%     end
% end
% fprintf(fid,'\t\t\t Arrays.asList(');
% for i = 1:sn.nnodes
%     fprintf(fid,'node%d',i);
%     if i<sn.nnodes
%         fprintf(fid,', ');
%     else
%         fprintf(fid,'));\n');
%     end
% end
fprintf(fid,'\t\n');


for k = 1:sn.nclasses
    for c = 1:sn.nclasses
        for i=1:sn.nnodes
            for m=1:sn.nnodes
                % routing matrix for each class
                myP{k,c}(i,m) = rtnodes((i-1)*sn.nclasses+k,(m-1)*sn.nclasses+c);
                if myP{k,c}(i,m) > 0
                    % do not change %d into %f to avoid round-off errors in
                    % the total probability
                    switch sn.nodetype(i)
                        case NodeType.Sink
                            % no-op
                        case NodeType.Fork
                            fprintf(fid,'\t\troutingMatrix.set(jobclass%d, jobclass%d, node%d, node%d, %f); // (%s,%s) -> (%s,%s)\n',k,c,i,m,sn.nodeparam{i}.fanOut,sn.nodenames{i},sn.classnames{k},sn.nodenames{m},sn.classnames{c});
                        otherwise
                            fprintf(fid,'\t\troutingMatrix.set(jobclass%d, jobclass%d, node%d, node%d, %f); // (%s,%s) -> (%s,%s)\n',k,c,i,m,myP{k,c}(i,m),sn.nodenames{i},sn.classnames{k},sn.nodenames{m},sn.classnames{c});
                    end
                end
            end
        end
        %fprintf(fid,'\t\tP{%d,%d} = %s;\n',k,c,mat2str(myP{k,c}));
    end
end


fprintf(fid,'\n\t\tmodel.link(routingMatrix);\n\n');
if headers
    fprintf(fid,'\t\treturn model;\n');
    fprintf(fid,'\t}\n');
end
%if fid~=1
%    fclose(fid);
%end
end
