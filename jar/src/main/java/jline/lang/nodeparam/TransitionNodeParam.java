package jline.lang.nodeparam;

import jline.lang.Mode;
import jline.lang.NodeParam;
import jline.lang.constant.ProcessType;
import jline.lang.constant.TimingStrategy;
import jline.util.SerializableFunction;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import java.util.List;
import java.util.Map;

public class TransitionNodeParam extends ServiceNodeParam {
    public List<Matrix> enabling;
    public List<Matrix> inhibiting;
    public List<String> modenames;
    public Matrix nmodeservers;
    public List<TimingStrategy> timing;
    public Map<Mode, ProcessType> firingprocid;
    public Map<Mode, MatrixCell> firingproc;
    public Map<Mode, Matrix> firingpie;
    public Matrix firingphases;
    public List<Matrix> firing;
    public Matrix firingprio;
    public Matrix fireweight;
    // Marking-dependent firing-rate multiplier g_mode(marking) per mode; a null
    // entry == unit multiplier. Argument is the node-indexed input-place marking.
    public List<SerializableFunction<Matrix, Double>> firingdep;
    public int nmodes;

    @Override
    public boolean isEmpty() {
        return super.isEmpty() &&
                (enabling == null || enabling.isEmpty()) &&
                (inhibiting == null || inhibiting.isEmpty()) &&
                (modenames == null || modenames.isEmpty()) &&
                nmodeservers == null &&
                (timing == null || timing.isEmpty()) &&
                (firingprocid == null || firingprocid.isEmpty()) &&
                (firingproc == null || firingproc.isEmpty()) &&
                (firingpie == null || firingpie.isEmpty()) &&
                firingphases == null &&
                (firing == null || firing.isEmpty()) &&
                firingprio == null &&
                fireweight == null &&
                (firingdep == null || firingdep.isEmpty()) &&
                nmodes == 0;
    }
}
