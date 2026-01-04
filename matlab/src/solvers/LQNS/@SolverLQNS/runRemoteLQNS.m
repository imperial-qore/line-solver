function result = runRemoteLQNS(self, filename, options)
% RUNREMOTELQNS Execute LQNS via remote REST API
%
% @brief Executes LQNS solver via a remote Docker container running lqns-rest API
% @param filename Path to the LQNX model file (without extension)
% @param options Solver options including remote configuration
% @return result Response from the remote server
%
% This method sends the LQNX model to a remote lqns-rest server for execution.
% The server URL is configured via options.config.remote_url.
% Upon successful execution, the LQXO response is written to disk for parsing.

% Read LQNX model file
modelFile = [filename, '.lqnx'];
if ~exist(modelFile, 'file')
    % Try without extension if already includes it
    modelFile = filename;
end
modelContent = fileread(modelFile);

% Build API URL
baseUrl = options.config.remote_url;
if ~endsWith(baseUrl, '/')
    baseUrl = [baseUrl, '/'];
end

% Determine endpoint based on method
switch options.method
    case {'sim', 'lqsim'}
        url = [baseUrl, 'api/v1/solve/lqsim'];
    otherwise
        url = [baseUrl, 'api/v1/solve/lqns'];
end

% Build request structure
request = struct();
request.model = struct();
request.model.content = modelContent;
request.model.base64 = false;

% Add options
request.options = struct();
request.options.include_raw_output = true;

% Add pragmas based on method
switch options.method
    case {'sim', 'lqsim'}
        request.options.blocks = 30;
        if options.samples > 0
            request.options.run_time = options.samples;
        end
    otherwise
        request.options.pragmas = struct();

        % Multiserver pragma
        switch options.config.multiserver
            case 'default'
                request.options.pragmas.multiserver = 'rolia';
            otherwise
                request.options.pragmas.multiserver = options.config.multiserver;
        end

        % Method-specific pragmas
        switch options.method
            case 'srvn'
                request.options.pragmas.layering = 'srvn';
            case 'exactmva'
                request.options.pragmas.mva = 'exact';
            case 'srvn.exactmva'
                request.options.pragmas.layering = 'srvn';
                request.options.pragmas.mva = 'exact';
        end

        request.options.pragmas.stop_on_message_loss = false;
end

% Set HTTP options
webOpts = weboptions(...
    'MediaType', 'application/json', ...
    'Timeout', 300, ...
    'ContentType', 'json');

% Make HTTP request
try
    response = webwrite(url, request, webOpts);
catch ME
    line_error(mfilename, 'Remote LQNS execution failed: %s', ME.message);
    result = struct('status', 'error', 'error', ME.message);
    return;
end

% Check response status
if isfield(response, 'status')
    if strcmp(response.status, 'error') || strcmp(response.status, 'failed')
        errorMsg = '';
        if isfield(response, 'error')
            errorMsg = response.error;
        end
        line_error(mfilename, 'Remote solver returned error: %s', errorMsg);
        result = response;
        return;
    end
end

% Write LQXO response to file for parsing by existing code
if isfield(response, 'raw_output') && isstruct(response.raw_output)
    if isfield(response.raw_output, 'lqxo') && ~isempty(response.raw_output.lqxo)
        lqxoFile = [filename, '.lqxo'];
        fid = fopen(lqxoFile, 'w');
        if fid ~= -1
            fprintf(fid, '%s', response.raw_output.lqxo);
            fclose(fid);
        else
            line_warning(mfilename, 'Could not write LQXO file: %s', lqxoFile);
        end
    end
end

result = response;
end
