/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang.constant;

import java.io.Serializable;

/**
 * Constants for specifying a type of metric
 */
public enum MetricType implements Serializable {
    ResidT,
    RespT,
    DropRate,
    QLen,
    QueueT,
    FCRWeight,
    FCRMemOcc,
    FJQLen,
    FJRespT,
    RespTSink,
    SysQLen,
    SysRespT,
    SysTput,
    Tput,
    ArvR,
    TputSink,
    Util,
    TranQLen,
    TranUtil,
    TranTput,
    TranRespT,
    Tard,
    SysTard;

    public static MetricType toMetricType(String desc) {
        switch (desc) {
            case "Residence Time":  // Response Time * Visits
            case "Residence time":  // JMVA uses lowercase
                return ResidT;
            case "Response Time":   // Response Time for a single Visit
                return RespT;
            case "Drop Rate":
                return DropRate;
            case "Number of Customers":
                return QLen;
            case "Queue Time":
                return QueueT;
            case "FCR Total Weight":
                return FCRWeight;
            case "FCR Memory Occupation":
                return FCRMemOcc;
            case "Fork Join Response Time":
                return FJRespT;
            case "Fork Join Number of Customers":
                return FJQLen;
            case "Response Time per Sink":
                return RespTSink;
            case "System Number of Customers":
                return SysQLen;
            case "System Response Time":
                return SysRespT;
            case "System Throughput":
                return SysTput;
            case "Throughput":
                return Tput;
            case "Arrival Rate":
                return ArvR;
            case "Throughput per Sink":
                return TputSink;
            case "Utilization":
                return Util;
            case "Tran Number of Customers":
                return TranQLen;
            case "Tran Utilization":
                return TranUtil;
            case "Tran Throughput":
                return TranTput;
            case "Tran Response Time":
                return TranRespT;
            case "Tardiness":
                return Tard;
            case "System Tardiness":
                return SysTard;
            default:
                return null;
        }
    }
}
