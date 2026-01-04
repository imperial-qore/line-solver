classdef SelfLoopingClass < ClosedClass
    % A class of jobs that perpetually cycle at its reference station.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods

        %Constructor
        function self = SelfLoopingClass(model, name, njobs, refstat, prio)
            % SELF = SELFLOOPINGCLASS(MODEL, NAME, NJOBS, REFSTAT, PRIO)
            self@ClosedClass(model, name, njobs, refstat, prio);
            if isa(model,'JNetwork')
                self.obj = jline.lang.SelfLoopingClass(model.obj, name, njobs, refstat.obj, prio);
            end
        end
    end

end

