function U = getEntryServiceMatrix(self)
% Matrix that returns the entry servt after multiplication with residt
% of entries and activities
eshift = self.lqn.eshift;
U = sparse(self.lqn.nidx + self.lqn.ncalls, self.lqn.nidx + self.lqn.ncalls);
for e = 1:self.lqn.nentries
    eidx = eshift + e;
    U = getEntryServiceMatrixRecursion(self.lqn, eidx, eidx, U);
end
U = double(U > 0);
end

function U = getEntryServiceMatrixRecursion(lqn, aidx, eidx, U)
% auxiliary function to getServiceMatrix
nextaidxs = find(lqn.graph(aidx,:));
G = double(full(lqn.graph)>0); % connection matrix between LQN elements
for nextaidx = nextaidxs
    if ~isempty(nextaidx)
        isLoop = false;
        % in the activity graph, the following if is entered only
        % by an edge that is the return from a LOOP activity
        if (lqn.graph(aidx,nextaidx) ~= lqn.dag(aidx,nextaidx))
            isLoop = true;
        end
        if ~(lqn.parent(aidx) == lqn.parent(nextaidx))
            %% if the successor activity is a call
            for cidx = lqn.callsof{aidx}
                switch lqn.calltype(cidx)
                    case CallType.SYNC
                        U(eidx,lqn.nidx+cidx) = 1; % mean number of calls already factored in
                    case CallType.ASYNC
                        % nop - doesn't contribute to respt
                end
            end
        end
        % here we have processed all calls, let us do the activities now
        %% if the successor activity is not a call
        if (lqn.parent(aidx) == lqn.parent(nextaidx))
            if nextaidx ~= aidx
                if ~isLoop
                    eidx=full(eidx);
                    U(eidx,nextaidx) = U(eidx,nextaidx) + G(aidx,nextaidx); % self-loops not included as already accounted in visits
                    %% now recursively build the rest of the routing matrix graph
                    U = getEntryServiceMatrixRecursion(lqn, nextaidx, eidx, U);
                end
            end
        end
    end
end % nextaidx
end
