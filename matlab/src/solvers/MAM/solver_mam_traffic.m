function ARV = solver_mam_traffic(sn, DEP, config)
% ARV = SOLVER_MAM_ESTFLOWS(QN, DEP, CONFIG)
% DEP{i,r} is the departure process of class r from i in (D0,D1) format

I = sn.nnodes;
C = sn.nchains;
R = sn.nclasses;

% In this function we use indexing over all non-ClassSwitch nodes
non_cs_classes = [];
isNCS = zeros(1,I);
nodeToNCS = zeros(1,I);
for ind=1:I
    if sn.nodetype(ind) ~= NodeType.ClassSwitch
        non_cs_classes(end+1:end+R)= ((ind-1)*R+1):(ind*R);
        isNCS(ind) = true;
        nodeToNCS(ind) = sum(isNCS);
    else
        isNCS(ind) = false;
    end
end

% Hide the nodes that are not class switches
rtncs = dtmc_stochcomp(sn.rtnodes,non_cs_classes);
Inc = I - sum(sn.nodetype == NodeType.ClassSwitch);

MMAP = DEP; % PH renewal process initially in (D0,D1) format

% We now bring into MMAP format with a single class
for ist=1:size(MMAP,1)
    for r=1:size(MMAP,2)
        if isempty(MMAP{ist,r}) || any(any(isnan(MMAP{ist,r}{1})))
            MMAP{ist,r} = {[0],[0],[0]}; % no arrivals from this class
        else
            MMAP{ist,r}{3} = MMAP{ist,r}{2};
        end
    end
end

ARV = cell(Inc,1);
DEP = cell(Inc,R);
LINKS = cell(Inc,Inc);
% first we determine all outgoing flows from all stations
for ind=1:I
    if isNCS(ind)
        inc = nodeToNCS(ind);
        switch sn.nodetype(ind)
            case {NodeType.Source, NodeType.Delay, NodeType.Queue}
                ist = sn.nodeToStation(ind);
                
                % obtain departure maps
                if R>1
                    DEP{inc} = mmap_super({MMAP{ist,1:R}});
                else
                    DEP{inc} = MMAP{ist,1};
                end
                Psplit = zeros(R,Inc*R);
                for r=1:R % superpose all classes
                    for jnd = 1:I %to
                        if isNCS(jnd)
                            jnc = nodeToNCS(jnd);
                            for s=1:R %to
                                Psplit(r,(jnc-1)*R+s) = rtncs((inc-1)*R+r, (jnc-1)*R+s);
                            end
                        end
                    end
                end
                
                [Fsplit{1:Inc}] = npfqn_traffic_split_cs(DEP{inc}, Psplit, config);
                for jnc=1:Inc
                    LINKS{inc,jnc} = Fsplit{jnc};
                    LINKS{inc,jnc} = mmap_normalize(LINKS{inc,jnc});
                end
        end
    end
end
% then we determine all incoming flows from all stations
for ind=1:I
    FLOWS={};
    if isNCS(ind) && sn.nodetype(ind) ~= NodeType.Source
        inc = nodeToNCS(ind);
        for jnd=1:Inc
            if ~isempty(LINKS{jnd,inc}) && sum(mmap_lambda(LINKS{jnd,inc}))> GlobalConstants.FineTol
                FLOWS{end+1} = LINKS{jnd,inc};
            end
        end
        if length(FLOWS)>1
            ARV{ind} = npfqn_traffic_merge(FLOWS, config);
        elseif length(FLOWS)==1
            ARV{ind} = FLOWS{1};
        else %length(FLOWS)==0
            ARV{ind} = LINKS{jnd,1}; % all links are zeros, take one
        end
    else
        ARV{ind} = [];
    end
end
end
