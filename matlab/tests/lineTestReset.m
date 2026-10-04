function lineTestReset()
% LINETESTRESET Clear the suite-wide runtests accumulator.
%
%   Called at the start of allTests so results do not carry over between runs
%   in a reused MATLAB session.
%
%   Root appdata rather than a global; see LINETESTACCUM for why.
%
%   See also LINETESTACCUM, LINETESTASSERT.

setappdata(0, 'LINETestResults', []);
end
