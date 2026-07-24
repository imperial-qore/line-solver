function [outspace, outrate, outprob, eventCache] = afterEventTransition(sn, ind, inspace, K, Ks, event, class, isSimulation, eventCache, R, space_buf, space_srv, space_fired, space_var, key)
% job arrives in class, then reads and moves into hit or miss
% class, then departs

outspace = [];
outrate = [];
outprob = 1;
switch event
    case EventType.ENABLE
        % no-op this is a global event
    case EventType.FIRE
        % no-op this is a global event
    case EventType.PHASE
        mode = class; % in a Transition, Event.class is interpreted as the mode
        outspace = [];
        outrate = [];
        outprob = [];
        [ni,nir,~,kir] = State.toMarginal(sn,ind,inspace,K,Ks,space_buf,space_srv,space_var);
        % ni, nir, kir is the count of *enabled* servers while
        % they progress their execution through the phases
        if nir(mode)>0
            for k=1:K(mode)
                en = space_srv(:,Ks(mode)+k) > 0;
                if any(en)
                    for kdest=setdiff(1:K(mode),k) % new phase
                        rate = 0;
                        space_srv_k = space_srv(en,:);
                        space_buf_k = space_buf(en,:);
                        space_fired_k = space_fired(en,:);
                        space_var_k = space_var(en,:);
                        % MAP transitions currently unsupported
                        %if ismkvmodclass(class)
                        %    space_var_k(sum(sn.nvars(ind,1:class))) = kdest;
                        %end
                        space_srv_k(:,Ks(mode)+k) = space_srv_k(:,Ks(mode)+k) - 1;
                        space_srv_k(:,Ks(mode)+kdest) = space_srv_k(:,Ks(mode)+kdest) + 1;

                        rate = sn.nodeparam{ind}.firingproc{mode}{1}(k,kdest)*kir(:,mode,k); % assume active
                        % We do not support class-dependence
                        outrate = [outrate; nir(mode).*rate];
                        outspace = [outspace; space_buf_k, space_srv_k, space_fired_k, space_var_k];
                        outprob = [outprob; ones(size(rate,1),1)];
                    end
                end
            end
            if isSimulation
                if nargin>=7 && isobject(eventCache)
                    eventCache(key) = {outprob, outspace,outrate};
                end

                if size(outspace,1) > 1
                    tot_rate = sum(outrate);
                    cum_rate = cumsum(outrate) / tot_rate;
                    firing_ctr = 1 + max([0,find( rand > cum_rate' )]); % select action
                    outspace = outspace(firing_ctr,:);
                    outrate = sum(outrate);
                    outprob = outprob(firing_ctr,:);
                end
            end
        end
end
end

