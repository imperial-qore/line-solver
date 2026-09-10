classdef (Sealed) CallType
    properties (Constant)
        SYNC = 1;
        ASYNC = 2;
        FWD = 3;
    end

    methods (Static)

        function txt = toText(callType)
            switch callType
                case CallType.SYNC
                    txt = 'Synchronous';
                case CallType.ASYNC
                    txt = 'Asynchronous';
                case CallType.FWD
                    txt = 'Forwarding';
                otherwise
                    line_error(mfilename, 'Unknown call type.');
            end
        end
    end
end
