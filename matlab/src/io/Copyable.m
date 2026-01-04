classdef Copyable < handle
    % Copyable allows to perform deep-copy of objects via the copy() method.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function newObj = copy(obj)
            % NEWOBJ = COPY(OBJ)

                objByteArray = getByteStreamFromArray(obj);
                newObj = getArrayFromByteStream(objByteArray);
        end
    end          
end