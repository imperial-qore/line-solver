import jline.io.CLI;
import jline.lang.constant.GlobalConstants;

import static jline.io.CLI.solverLN;

public class JarMain {
    public static void main(String[] args) throws Exception {
        System.out.format("LINE Java Solver "+ GlobalConstants.Version +" - Command Line Interface\n");
        System.out.format("Copyright (c) 2012-2024, QORE group, Imperial College London\n");
        System.out.format("-------------------------\n");

        CLI.solverLN(args);
    }
}
