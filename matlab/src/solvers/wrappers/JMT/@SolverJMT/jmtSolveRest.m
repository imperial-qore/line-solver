function [status, cmdout] = jmtSolveRest(self, restUrl, mode, fname, seed, options)
% [STATUS, CMDOUT] = JMTSOLVEREST(RESTURL, MODE, FNAME, SEED, OPTIONS)
%
% Solve through a JMT REST server (the imperialqore/jmt-rest container). The
% wire format carries the same JSIM or JMVA document the CLI reads, and the
% response carries the same result document the CLI writes, which is written
% to [FNAME,'-result.jsim'] or [FNAME,'-result.jmva'] so that the local result
% parsers are unaffected by which backend ran.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch mode
    case 'sim'
        route = 'sim';
        resultExt = 'jsim';
    case 'mva'
        route = 'mva';
        resultExt = 'jmva';
    otherwise
        line_error(mfilename, sprintf('Unknown JMT analysis mode: %s', mode));
end

url = regexprep(char(restUrl), '/+$', '');
if isempty(regexp(url, '/api/v\d+/solve/(sim|mva)$', 'once'))
    url = [url, '/api/v1/solve/', route];
end

req = struct();
req.model = struct('content', fileread(fname), 'base64', false);
% JMVA takes its algorithm and tolerance from the model document, so the seed
% is only meaningful for the simulation route; sending it to /solve/mva would
% be rejected by the server's option allow-list.
if strcmp(mode, 'sim')
    req.options = struct('seed', seed);
end

timeout = 3600;
if isfield(options, 'timeout') && ~isempty(options.timeout) && isfinite(options.timeout)
    timeout = options.timeout;
end
wopts = weboptions('MediaType', 'application/json', 'ContentType', 'json', ...
    'RequestMethod', 'post', 'Timeout', timeout);

if options.verbose == VerboseLevel.DEBUG
    line_printf('JMT REST: POST %s\n', url);
end

try
    resp = webwrite(url, req, wopts);
catch ME
    line_error(mfilename, sprintf('JMT REST request to %s failed: %s', url, ME.message));
end

if ~isstruct(resp) || ~isfield(resp, 'status')
    line_error(mfilename, sprintf('JMT REST server at %s returned an unexpected payload.', url));
end
if ~strcmpi(resp.status, 'completed')
    msg = 'unspecified error';
    if isfield(resp, 'error') && ~isempty(resp.error)
        msg = resp.error;
    end
    line_error(mfilename, sprintf('JMT REST solve failed: %s', strtrim(char(msg))));
end
if ~isfield(resp, 'raw_output') || ~isstruct(resp.raw_output) ...
        || ~isfield(resp.raw_output, 'result_xml') || isempty(resp.raw_output.result_xml)
    line_error(mfilename, ['JMT REST response carries no result document. The server was ', ...
        'asked to include the raw output; check that include_raw_output is not disabled.']);
end

resultPath = [fname, '-result.', resultExt];
fid = fopen(resultPath, 'w');
if fid < 0
    line_error(mfilename, sprintf('Cannot write the JMT result file %s', resultPath));
end
fprintf(fid, '%s', resp.raw_output.result_xml);
fclose(fid);

status = 0;
cmdout = '';
if isfield(resp, 'raw_output') && isfield(resp.raw_output, 'stdout') && ~isempty(resp.raw_output.stdout)
    cmdout = char(resp.raw_output.stdout);
end
end
