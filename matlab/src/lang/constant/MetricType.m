classdef (Sealed) MetricType
    % An output metric of a Solver, such as a performance index
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ResidT = 'Residence Time'; % Response Time * Visits
        RespT = 'Response Time'; % Response Time for a single Visit
        DropRate = 'Drop Rate';
        QLen = 'Number of Customers';
        QueueT = 'Queue Time';
        FCRWeight = 'FCR Total Weight';
        FCRMemOcc = 'FCR Memory Occupation';
        FJQLen = 'Fork Join Response Time';
        FJRespT = 'Fork Join Response Time';
        RespTSink = 'Response Time per Sink';
        SysDropR = 'System Drop Rate';
        SysQLen = 'System Number of Customers';
        SysPower = 'System Power';
        SysRespT = 'System Response Time';
        SysTput = 'System Throughput';
        Tput = 'Throughput';
        ArvR = 'Arrival Rate';
        TputSink = 'Throughput per Sink';
        Util = 'Utilization';
        TranQLen = 'Tran Number of Customers';
        TranUtil = 'Tran Utilization';
        TranTput = 'Tran Throughput';
        TranRespT = 'Tran Response Time';
        
        ID_ResidT = 0; % Response Time * Visits
        ID_RespT = 1; % Response Time for one Visit
        ID_DropRate = 2;
        ID_QLen = 3;
        ID_QueueT = 4;
        ID_FCRWeight = 5;
        ID_FCRMemOcc = 6;
        ID_FJQLen = 7;
        ID_FJRespT = 8;
        ID_RespTSink = 9;
        ID_SysDropR = 10;
        ID_SysQLen = 11;
        ID_SysPower = 12;
        ID_SysRespT = 13;
        ID_SysTput = 14;
        ID_Tput = 15;
        ID_ArvR = 16;
        ID_TputSink = 17;
        ID_Util = 18;
        ID_TranQLen = 19;
        ID_TranUtil = 20;
        ID_TranTput = 21;
        ID_TranRespT = 22;
    end
        
end

