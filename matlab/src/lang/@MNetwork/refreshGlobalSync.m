function gsync = refreshGlobalSync(self)
% SYNC = REFRESHGLOBALSYNC()

sn = self.sn;
local = self.getNumberOfNodes+1;
nclasses = sn.nclasses;
gsync = {};
emptystate = cellzeros(sn.nnodes,1,0,0);
if any(sn.isstatedep(:))
    rtmask = self.sn.rtfun(emptystate, emptystate);
else
    rtmask = ceil(self.sn.rt);
end

for ind=1:sn.nnodes
    for r=1:sn.nclasses

        if sn.isstateful(ind)
            if sn.nodetype(ind) == NodeType.Transition
                for m=1:sn.nodeparam{ind}.nmodes
                    % mode enabling
                    enablingPlaces = find(sn.nodeparam{ind}.enabling{m});
                    gsync{end+1,1}.active{1} = ModeEvent(EventType.ENABLE, ind, m, 1.0);
                    gsync{end,1}.passive = cell(1,length(enablingPlaces));
                    for ep=1:length(enablingPlaces)
                        gsync{end,1}.passive{ep} = ModeEvent(EventType.LOCAL, enablingPlaces(ep), m, 1.0); % ID_LOCAL has no state effects
                    end
                end
                for m=1:sn.nodeparam{ind}.nmodes
                    % mode firing
                    firingPlaces = find(sn.nodeparam{ind}.firing{m});
                    enablingPlaces = find(sn.nodeparam{ind}.enabling{m});
                    gsync{end+1,1}.active{1} = ModeEvent(EventType.FIRE, ind, m);

                    gsync{end,1}.passive = {};
                    for ep=1:length(enablingPlaces)
                        % TODO: this creates a departure event for each
                        % job pulled from the enabling places, which is
                        % inefficient
                        gsync{end,1}.passive{end+1} = ModeEvent(EventType.PRE, enablingPlaces(ep), m, sn.nodeparam{ind}.enabling{m}(enablingPlaces(ep)));
                    end
                    for fp=1:length(firingPlaces)
                        % TODO: this creates an arrival event for each
                        % fired job to the destination place, which is
                        % inefficient
                        gsync{end,1}.passive{end+1} = ModeEvent(EventType.POST, firingPlaces(fp), m, sn.nodeparam{ind}.firing{m}(firingPlaces(fp)));
                    end
                end
            end

        end
    end
end
if ~isempty(self.sn) %&& isprop(self.sn,'nvars')
    self.sn.gsync = gsync;
end
end
