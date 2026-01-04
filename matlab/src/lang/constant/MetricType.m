classdef (Sealed) MetricType
    % An output metric of a Solver, such as a performance index
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        ResidT = 0;
        RespT = 1;
        DropRate = 2;
        QLen = 3;
        QueueT = 4;
        FCRWeight = 5;
        FCRMemOcc = 6;
        FJQLen = 7;
        FJRespT = 8;
        RespTSink = 9;
        SysQLen = 10;
        SysRespT = 11;
        SysTput = 12;
        Tput = 13;
        ArvR = 14;
        TputSink = 15;
        Util = 16;
        TranQLen = 17;
        TranUtil = 18;
        TranTput = 19;
        TranRespT = 20;
        Tard = 21;
        SysTard = 22;
    end

    methods (Static)
        function txt = toText(metricId)
            switch metricId
                case MetricType.ResidT
                    txt = 'Residence Time';
                case MetricType.RespT
                    txt = 'Response Time';
                case MetricType.DropRate
                    txt = 'Drop Rate';
                case MetricType.QLen
                    txt = 'Number of Customers';
                case MetricType.QueueT
                    txt = 'Queue Time';
                case MetricType.FCRWeight
                    txt = 'FCR Total Weight';
                case MetricType.FCRMemOcc
                    txt = 'FCR Memory Occupation';
                case MetricType.FJQLen
                    txt = 'Fork Join Number of Customers';
                case MetricType.FJRespT
                    txt = 'Fork Join Response Time';
                case MetricType.RespTSink
                    txt = 'Response Time per Sink';
                case MetricType.SysQLen
                    txt = 'System Number of Customers';
                case MetricType.SysRespT
                    txt = 'System Response Time';
                case MetricType.SysTput
                    txt = 'System Throughput';
                case MetricType.Tput
                    txt = 'Throughput';
                case MetricType.ArvR
                    txt = 'Arrival Rate';
                case MetricType.TputSink
                    txt = 'Throughput per Sink';
                case MetricType.Util
                    txt = 'Utilization';
                case MetricType.TranQLen
                    txt = 'Tran Number of Customers';
                case MetricType.TranUtil
                    txt = 'Tran Utilization';
                case MetricType.TranTput
                    txt = 'Tran Throughput';
                case MetricType.TranRespT
                    txt = 'Tran Response Time';
                case MetricType.Tard
                    txt = 'Tardiness';
                case MetricType.SysTard
                    txt = 'System Tardiness';
                otherwise
                    txt = 'Unknown Metric';
            end
        end
    end
end
