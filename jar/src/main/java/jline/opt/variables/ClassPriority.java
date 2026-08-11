package jline.opt.variables;

import jline.lang.JobClass;
import jline.lang.Network;

import java.util.ArrayList;
import java.util.List;

/**
 * Optimize the priority of job classes. In 'levels' mode each class gets an
 * integer priority in [minPriority, maxPriority] (one encoded dimension per
 * class); in 'permutation' mode the encoded random keys induce a priority
 * ordering (n-1 dimensions). Mirrors native-Python {@code ClassPriority}.
 */
public class ClassPriority extends DecisionVariable {

    private final List<JobClass> jobclasses;
    private final String mode;   // "levels" or "permutation"
    private final int minPriority;
    private final int maxPriority;

    public ClassPriority(List<JobClass> jobclasses) {
        this(jobclasses, "levels", 1, 10, "class_priorities");
    }

    public ClassPriority(List<JobClass> jobclasses, String mode, int minPriority,
                         int maxPriority, String name) {
        super(name);
        this.jobclasses = jobclasses;
        this.mode = mode;
        this.minPriority = minPriority;
        this.maxPriority = maxPriority;
        if ("levels".equals(mode)) {
            this.dimension = jobclasses.size();
        } else {
            this.dimension = Math.max(1, jobclasses.size() - 1);
        }
    }

    public List<JobClass> getJobClasses() {
        return jobclasses;
    }

    public String getMode() {
        return mode;
    }

    public double[][] getBounds() {
        return unitBounds(dimension);
    }

    public Object decode(double[] x) {
        if ("levels".equals(mode)) {
            int[] priorities = new int[x.length];
            for (int i = 0; i < x.length; i++) {
                double p = minPriority + x[i] * (maxPriority - minPriority);
                priorities[i] = (int) Math.round(p);
            }
            return priorities;
        }
        // permutation from random keys, descending order (argsort of -keys)
        int n = jobclasses.size();
        if (n == 1) {
            return new int[]{0};
        }
        double[] keys = new double[n];
        for (int i = 0; i < x.length; i++) {
            keys[i] = x[i];
        }
        keys[n - 1] = 0.0;
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) {
            order[i] = i;
        }
        // stable descending sort by key (numpy argsort of -keys is stable ascending
        // of -keys == descending of keys, ties by original index)
        java.util.Arrays.sort(order, new java.util.Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                int c = Double.compare(keys[b], keys[a]);
                if (c != 0) {
                    return c;
                }
                return Integer.compare(a, b);
            }
        });
        int[] out = new int[n];
        for (int i = 0; i < n; i++) {
            out[i] = order[i];
        }
        return out;
    }

    public void apply(Network model, Object value) {
        int[] v = (int[]) value;
        if ("levels".equals(mode)) {
            for (int i = 0; i < jobclasses.size(); i++) {
                JobClass cls = resolveOrSelf(model, jobclasses.get(i));
                cls.setPriority(v[i]);
            }
        } else {
            for (int rank = 0; rank < v.length; rank++) {
                int classIdx = v[rank];
                int priority = v.length - rank;
                JobClass cls = resolveOrSelf(model, jobclasses.get(classIdx));
                cls.setPriority(priority);
            }
        }
    }

    private JobClass resolveOrSelf(Network model, JobClass jc) {
        JobClass r = resolveClass(model, jc);
        return r != null ? r : jc;
    }

    public String getVariableType() {
        return "class_priority";
    }

    // reference to keep import used when list ops needed
    static List<JobClass> asList(JobClass... cs) {
        List<JobClass> l = new ArrayList<JobClass>();
        for (JobClass c : cs) {
            l.add(c);
        }
        return l;
    }
}
