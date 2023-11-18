classdef (Sealed) CallType
    properties (Constant)
        SYNC = 'Synchronous';
        ASYNC = 'Asynchronous';
        FWD = 'Forwarding';
        
        ID_SYNC = 1;
        ID_ASYNC = 2;
        ID_FWD = 3;
    end
end