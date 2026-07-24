classdef Trace < Replayer
    % Empirical time series from a trace, alias for Replayer
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
        
    methods
        %Constructor
        function self = Trace(data)
            % SELF = TRACE(data)
            self@Replayer(data);
        end
        
        function [m1,m2,m3,scv,skew] = getMoments(self)
            data = self.data;
            [row,col] = size(data);
            m1 = 0; % the first moment
            m2 = 0; % the second moment
            m3 = 0; % the third moment
            for i = 1:1:row-1
                bin1 = ((data(i+1,2)-data(i,2))/2+data(i,2))*((data(i+1,1)-data(i,1)));
                bin2 = ((data(i+1,2)-data(i,2))/2+data(i,2))^2*((data(i+1,1)-data(i,1)));
                bin3 = ((data(i+1,2)-data(i,2))/2+data(i,2))^3*((data(i+1,1)-data(i,1)));
                m1 = m1+bin1;
                m2 = m2+bin2;
                m3 = m3+bin3;
            end
            scv = (m2/m1^2)-1;
            skew = (m3-3*m1*(m2-m1^2)-m1^3)/((m2-m1^2)^(3/2));
        end
    end
end

