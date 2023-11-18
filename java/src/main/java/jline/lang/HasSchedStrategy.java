package jline.lang;

import jline.lang.constant.SchedStrategy;
import jline.lang.distributions.Distribution;

public interface HasSchedStrategy {
    SchedStrategy getSchedStrategy();
    Distribution getServiceProcess(JobClass jobClass);
    double minRate();
    double maxRate();
    double avgRate();
    int rateCt();
}
