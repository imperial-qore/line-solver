function [stateSpace,stateSpaceAggr,stateSpaceHashed,nodeStateSpace,sn] = ctmc_ssg(sn,options)

[stateSpace,stateSpaceHashed,qnc] = State.spaceGenerator(sn, options.cutoff, options);
nodeStateSpace = qnc.space;
sn.space = nodeStateSpace;

% if options.verbose
%     line_printf('\nCTMC state space size: %d states. ',size(stateSpace,1));
% end
if ~isfield(options, 'hide_immediate')
    options.hide_immediate = true;
end

nstateful = sn.nstateful;
nclasses = sn.nclasses;
sync = sn.sync;
A = length(sync);
stateSpaceAggr = zeros(size(stateSpaceHashed));

% for all synchronizations
for a=1:A
    stateCell = cell(nstateful,1);
    for s=1:size(stateSpaceHashed,1)
        state = stateSpaceHashed(s,:);
        % update state cell array and SSq
        for ind = 1:sn.nnodes
            if sn.isstateful(ind)
                isf = sn.nodeToStateful(ind);
                stateCell{isf} = sn.space{isf}(state(isf),:);
                if sn.isstation(ind)
                    ist = sn.nodeToStation(ind);
                    [~,nir] = State.toMarginal(sn,ind,stateCell{isf});

                    stateSpaceAggr(s,((ist-1)*nclasses+1):ist*nclasses) = nir;
                end
            end
        end
    end
end
end