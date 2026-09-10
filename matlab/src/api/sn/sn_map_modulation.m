%{ @file sn_map_modulation.m
 %  @brief Extracts the Markov-modulation records of the MAP/MMPP/MMAP processes
 %
 %  @author LINE Development Team
%}

%{
 % @brief Collects the (D0,D1) modulation records of every non-renewal arrival
 % or service process declared in the network
 %
 % @details
 % A MAP with matrices (D0,D1) is a Poisson-like point process modulated by the
 % CTMC with generator Q = D0 + D1 (the phase process), whose conditional
 % intensity in phase k is lambda(k) = sum_j D1(k,j). This function returns one
 % record per modulating process, so that a solver-agnostic transformation can
 % replace each of them by a random-environment stage set (see map2renv).
 %
 % Only processes declared as MAP, MMPP2 or MMAP are reported: every other
 % distribution is stored in sn.proc in (D0,D1) form as well (Erlang, Coxian,
 % APH, ...), but those are renewal processes that carry no modulation and are
 % supported natively by the phase-type solvers.
 %
 % Marked processes (MMAP) at a Source are reported as a single record whose
 % classes vector lists every marked class, since all marks share one phase
 % process; the per-class intensity comes from the mark-specific D1 matrices.
 %
 % @par Syntax:
 % @code
 % mods = sn_map_modulation(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>mods<td>Struct array with fields ist, node, kind ('arrival'|'service'),
 %                 classes, D0, D1 (cell, one per entry of classes), order, isMMPP
 % </table>
%}
function mods = sn_map_modulation(sn)

mods = struct('ist',{},'node',{},'kind',{},'classes',{},'D0',{},'D1',{},'order',{},'isMMPP',{});

if ~isfield(sn,'procid') || isempty(sn.procid)
    return
end

modTypes = [ProcessType.MAP, ProcessType.MMPP2, ProcessType.MMAP];
hasMark = isfield(sn,'markidx') && ~isempty(sn.markidx);

for ist = 1:sn.nstations
    nd = sn.stationToNode(ist);
    if sn.nodetype(nd) == NodeType.Source
        kind = 'arrival';
    else
        kind = 'service';
    end
    done = false(1,sn.nclasses);
    for r = 1:sn.nclasses
        if done(r) || ~any(sn.procid(ist,r) == modTypes)
            continue
        end
        map_ir = sn_map_of(sn, ist, r);
        if isempty(map_ir)
            continue
        end
        if hasMark && sn.markidx(ist,r) > 0
            % MMAP: one phase process shared by every marked class of the
            % station, one D1 block per mark
            marked = find(sn.markidx(ist,:) > 0);
            carrier = marked(sn.markidx(ist,marked) == min(sn.markidx(ist,marked)));
            map_c = sn_map_of(sn, ist, carrier(1));
            if numel(map_c) < 2 + numel(marked)
                line_error(mfilename, sprintf(['The marked arrival process at station %d carries %d mark ' ...
                    'matrices for %d marked classes; the (D0,D1,D1^(1),...,D1^(C)) form is required.'], ...
                    ist, max(0,numel(map_c)-2), numel(marked)));
            end
            D1c = cell(1,numel(marked));
            for k = 1:numel(marked)
                D1c{k} = map_c{2 + sn.markidx(ist,marked(k))};
            end
            mods(end+1) = struct('ist',ist,'node',nd,'kind',kind,'classes',marked, ...
                'D0',map_c{1},'D1',{D1c},'order',size(map_c{1},1), ...
                'isMMPP',all(cellfun(@sn_is_diagonal, D1c))); %#ok<AGROW>
            done(marked) = true;
        else
            D1r = map_ir(2); % single-class list: one D1 block
            mods(end+1) = struct('ist',ist,'node',nd,'kind',kind,'classes',r, ...
                'D0',map_ir{1},'D1',{D1r},'order',size(map_ir{1},1), ...
                'isMMPP',sn_is_diagonal(map_ir{2})); %#ok<AGROW>
            done(r) = true;
        end
    end
end
end

function map_ir = sn_map_of(sn, ist, r)
% Per station-class (D0,D1,...) representation, [] when absent or disabled.
map_ir = [];
if isempty(sn.proc) || numel(sn.proc) < ist || isempty(sn.proc{ist})
    return
end
if numel(sn.proc{ist}) < r || isempty(sn.proc{ist}{r})
    return
end
map_ir = sn.proc{ist}{r};
if numel(map_ir) < 2 || isempty(map_ir{1}) || any(isnan(map_ir{1}(:)))
    map_ir = [];
end
end

function bool = sn_is_diagonal(D1)
% True when D1 has no off-diagonal mass, i.e. the process is an MMPP.
bool = norm(D1 - diag(diag(D1)), 'fro') <= GlobalConstants.Zero * max(1, norm(D1,'fro'));
end
