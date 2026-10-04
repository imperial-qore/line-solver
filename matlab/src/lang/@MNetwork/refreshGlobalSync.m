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

% see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync) for rationale
for ind=1:sn.nnodes

        if sn.isstateful(ind)
            if sn.nodetype(ind) == NodeType.Transition
                for m=1:sn.nodeparam{ind}.nmodes
                    % mode enabling
                    enablingPlaces = find(sn.nodeparam{ind}.enabling{m});
                    % inhibitor-arc input places (finite threshold). Added as
                    % LOCAL events (no state effect) so their markings are read
                    % into ep_space and the enabling test can apply the upper
                    % bound. Exclude places already listed as enabling inputs.
                    inhibitingPlaces = setdiff(find(~isinf(sn.nodeparam{ind}.inhibiting{m})), enablingPlaces);
                    gsync{end+1,1}.active{1} = ModeEvent(EventType.ENABLE, ind, m, 1.0);
                    gsync{end,1}.passive = cell(1,length(enablingPlaces)+length(inhibitingPlaces));
                    for ep=1:length(enablingPlaces)
                        gsync{end,1}.passive{ep} = ModeEvent(EventType.LOCAL, enablingPlaces(ep), m, 1.0); % ID_LOCAL has no state effects
                    end
                    for ip=1:length(inhibitingPlaces)
                        gsync{end,1}.passive{length(enablingPlaces)+ip} = ModeEvent(EventType.LOCAL, inhibitingPlaces(ip), m, 1.0);
                    end
                end
                for m=1:sn.nodeparam{ind}.nmodes
                    % mode firing
                    % see _kb/04-networkstruct.md (refreshGlobalSync/refreshSync) for rationale
                    firingPlaces = find(sn.nodeparam{ind}.firing{m} > 0);
                    enablingPlaces = find(sn.nodeparam{ind}.enabling{m});
                    % inhibitor-arc input places (finite threshold), read into
                    % ep_space for the firing enabling-degree test. Added as
                    % LOCAL events (no token effect); exclude places already
                    % present as enabling (PRE) or firing (POST) passives.
                    inhibitingPlaces = setdiff(find(~isinf(sn.nodeparam{ind}.inhibiting{m})), union(enablingPlaces, firingPlaces));
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
                    for ip=1:length(inhibitingPlaces)
                        gsync{end,1}.passive{end+1} = ModeEvent(EventType.LOCAL, inhibitingPlaces(ip), m, 1.0);
                    end
                end
            end

        end
end
if ~isempty(self.sn) %&& isprop(self.sn,'nvars')
    self.sn.gsync = gsync;
end
end
