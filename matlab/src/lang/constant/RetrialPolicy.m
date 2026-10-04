classdef (Sealed) RetrialPolicy
    % Enumeration of retrial policies for an orbiting population.
    %
    % The policy fixes how the aggregate rate at which the orbit attempts to
    % re-enter the station depends on the orbit size n:
    %
    %   LINEAR   - each orbiting customer retries independently at rate nu, so
    %              the aggregate retrial rate is n*nu. This is the classical
    %              retrial queue of Falin and Templeton.
    %   CONSTANT - the orbit as a whole retries at rate nu whenever it is
    %              non-empty, independently of n. This models a single retrial
    %              controller shared by the orbit rather than per-customer
    %              timers.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        LINEAR = 1;
        CONSTANT = 2;
    end

    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case RetrialPolicy.LINEAR
                    text = 'linear';
                case RetrialPolicy.CONSTANT
                    text = 'constant';
                otherwise
                    line_error(mfilename, 'Unrecognized retrial policy type.');
            end
        end

        function type = fromText(text)
            % TYPE = FROMTEXT(TEXT)
            switch lower(text)
                case {'linear','per-customer','falin'}
                    type = RetrialPolicy.LINEAR;
                case {'constant','fixed','controller'}
                    type = RetrialPolicy.CONSTANT;
                otherwise
                    line_error(mfilename, sprintf('Unrecognized retrial policy ''%s''.', text));
            end
        end
    end
end
