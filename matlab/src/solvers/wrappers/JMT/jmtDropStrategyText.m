function text = jmtDropStrategyText(sn, ist, r)
% TEXT = JMTDROPSTRATEGYTEXT(SN, IST, R)
%
% The JMT dropStrategy/dropRule string for station IST, class R.
%
% Beyond DropStrategy.toText this does two things.
%
% It resolves the two ways LINE can declare BAS blocking onto the one way JMT
% can read it. JMT's queue section says what happens to an arrival that finds
% THIS buffer full, so it only understands the rule on the destination; a WAITQ
% slot that jmtIsBasDestination marks is therefore written out as 'BAS blocking'.
%
% And it keeps the written file VALID. JMT recognizes exactly four strings --
% 'drop', 'BAS blocking', 'waiting queue', 'retrial' -- so BBS, RSRD and
% retrial-with-limit are spelled 'waiting queue', JMT's own no-limit default.
% That substitution is only ever reached where the rule cannot be consulted
% (infinite size, or a closed capacity equal to the population): a buffer that
% can actually fill under one of those three is refused outright by
% jmtStationCapAssert in saveBufferCapacity.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isnan(ist) || ist < 1
    text = 'drop'; % JMT sets the field to 'drop' for nodes without a buffer
    return
end
dr = sn.droprule(ist, r);
if dr == 0
    text = 'drop'; % the "no rule recorded" slot, which JMT also spells 'drop'
    return
end
if dr == DropStrategy.WAITQ && jmtIsBasDestination(sn, ist, r)
    text = DropStrategy.toText(DropStrategy.BAS);
    return
end
switch dr
    case {DropStrategy.BBS, DropStrategy.RSRD, DropStrategy.RETRIAL_WITH_LIMIT}
        text = DropStrategy.toText(DropStrategy.WAITQ);
        return
end
text = DropStrategy.toText(dr);
end
