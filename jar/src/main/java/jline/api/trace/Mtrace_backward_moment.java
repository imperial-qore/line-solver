package jline.api.trace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public final class Mtrace_backward_moment {
    private Mtrace_backward_moment() {}

    /**
     * Computes backward moments of a multi-class trace.
     * Backward moments characterize the time until the previous arrival.
     *
     * Returns a Double if classIndex specified, double[] if all classes.
     */
    public static Object mtrace_backward_moment(double[] interArrivalTimes, int[] classLabels,
                                                 int order, int classIndex) {
        if (interArrivalTimes.length != classLabels.length) {
            throw new IllegalArgumentException("Inter-arrival times and class labels must have same length");
        }
        if (order < 1) {
            throw new IllegalArgumentException("Moment order must be >= 1");
        }

        int max = 0;
        for (int i = 0; i < classLabels.length; i++) {
            if (classLabels[i] > max) max = classLabels[i];
        }
        int numClasses = max + 1;

        if (classIndex >= 0) {
            return Double.valueOf(computeClassBackwardMoment(interArrivalTimes, classLabels, classIndex, order));
        } else {
            double[] result = new double[numClasses];
            for (int c = 0; c < numClasses; c++) {
                result[c] = computeClassBackwardMoment(interArrivalTimes, classLabels, c, order);
            }
            return result;
        }
    }

    public static Object mtrace_backward_moment(double[] interArrivalTimes, int[] classLabels, int order) {
        return mtrace_backward_moment(interArrivalTimes, classLabels, order, -1);
    }

    /**
     * Compute backward moment for a specific class
     */
    private static double computeClassBackwardMoment(double[] interArrivalTimes, int[] classLabels,
                                                      int targetClass, int order) {
        List<Double> backwardTimes = new ArrayList<Double>();
        double cumulativeTime = 0.0;
        Map<Integer, Double> lastClassTime = new HashMap<Integer, Double>();

        for (int i = 0; i < interArrivalTimes.length; i++) {
            cumulativeTime += interArrivalTimes[i];
            int currentClass = classLabels[i];

            if (currentClass == targetClass) {
                double backwardTime;
                if (lastClassTime.containsKey(Integer.valueOf(targetClass))) {
                    backwardTime = cumulativeTime - lastClassTime.get(Integer.valueOf(targetClass)).doubleValue();
                } else {
                    backwardTime = cumulativeTime;
                }
                backwardTimes.add(Double.valueOf(backwardTime));
            }

            lastClassTime.put(Integer.valueOf(currentClass), Double.valueOf(cumulativeTime));
        }

        if (backwardTimes.isEmpty()) {
            return 0.0;
        }

        if (order == 1) {
            double sum = 0.0;
            for (int i = 0; i < backwardTimes.size(); i++) {
                sum += backwardTimes.get(i).doubleValue();
            }
            return sum / backwardTimes.size();
        } else {
            double sum = 0.0;
            for (int i = 0; i < backwardTimes.size(); i++) {
                sum += Math.pow(backwardTimes.get(i).doubleValue(), (double) order);
            }
            return sum / backwardTimes.size();
        }
    }

    /**
     * Computes conditional backward moments given the forward recurrence time.
     */
    public static double mtrace_backward_moment_conditional(double[] interArrivalTimes, int[] classLabels,
                                                             int conditioningClass, int order) {
        List<Double> conditionalBackwardTimes = new ArrayList<Double>();
        double cumulativeTime = 0.0;
        Map<Integer, List<Double>> classArrivalTimes = new HashMap<Integer, List<Double>>();

        for (int i = 0; i < interArrivalTimes.length; i++) {
            cumulativeTime += interArrivalTimes[i];
            int currentClass = classLabels[i];

            List<Double> existing = classArrivalTimes.get(Integer.valueOf(currentClass));
            if (existing == null) {
                existing = new ArrayList<Double>();
                classArrivalTimes.put(Integer.valueOf(currentClass), existing);
            }
            existing.add(Double.valueOf(cumulativeTime));
        }

        List<Double> conditioningArrivals = classArrivalTimes.get(Integer.valueOf(conditioningClass));
        if (conditioningArrivals == null) {
            return 0.0;
        }

        for (int ci = 0; ci < conditioningArrivals.size(); ci++) {
            double conditioningTime = conditioningArrivals.get(ci).doubleValue();
            for (Map.Entry<Integer, List<Double>> entry : classArrivalTimes.entrySet()) {
                int classIdx = entry.getKey().intValue();
                List<Double> arrivalTimes = entry.getValue();
                if (classIdx != conditioningClass) {
                    Double recentArrival = null;
                    for (int k = 0; k < arrivalTimes.size(); k++) {
                        double t = arrivalTimes.get(k).doubleValue();
                        if (t < conditioningTime) {
                            if (recentArrival == null || t > recentArrival.doubleValue()) {
                                recentArrival = Double.valueOf(t);
                            }
                        }
                    }
                    if (recentArrival != null) {
                        conditionalBackwardTimes.add(Double.valueOf(conditioningTime - recentArrival.doubleValue()));
                    }
                }
            }
        }

        if (conditionalBackwardTimes.isEmpty()) {
            return 0.0;
        }

        if (order == 1) {
            double sum = 0.0;
            for (int i = 0; i < conditionalBackwardTimes.size(); i++) {
                sum += conditionalBackwardTimes.get(i).doubleValue();
            }
            return sum / conditionalBackwardTimes.size();
        } else {
            double sum = 0.0;
            for (int i = 0; i < conditionalBackwardTimes.size(); i++) {
                sum += Math.pow(conditionalBackwardTimes.get(i).doubleValue(), (double) order);
            }
            return sum / conditionalBackwardTimes.size();
        }
    }
}
