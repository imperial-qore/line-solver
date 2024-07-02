package jline.io;

import jline.lang.constant.GlobalConstants;
import jline.lang.constant.VerboseLevel;

import java.util.logging.Level;
import java.util.logging.Logger;

public class InputOutput {

    private static String lastWarning = "";
    private static boolean suppressedWarnings = false;
    private static long suppressedWarningTic = System.currentTimeMillis();

    private static void line_printf(String message) {
        Logger.getLogger(InputOutput.class.getName()).log(Level.INFO, message);
    }

    public static void line_warning(String caller, String msg, Object... args) {
        Logger logger = Logger.getLogger(InputOutput.class.getName());

        // Check if global verbosity level allows warnings
        if (GlobalConstants.Verbose == VerboseLevel.SILENT) {
            return;
        }

        String errmsg = String.format(msg, args);
        String finalmsg = String.format("[%s] %s", caller, errmsg);

        try {
            long currentTime = System.currentTimeMillis();

            if (finalmsg.compareTo(lastWarning)!=0 || (currentTime - suppressedWarningTic) > 60000) {
                line_printf(finalmsg, Level.WARNING);
                lastWarning = finalmsg;
                suppressedWarnings = false;
                suppressedWarningTic = currentTime;
            } else {
                if (!suppressedWarnings) {
                    line_printf(String.format("Warning [%s]: %s", caller, "Warning message casted more than once, repetitions will not be printed on screen for 60 seconds.\n"), Level.WARNING);
                    suppressedWarnings = true;
                    suppressedWarningTic = currentTime;
                }
            }
        } catch (Exception e) {
            logger.log(Level.SEVERE, "Exception in line_warning", e);
        }
    }

    private static void line_printf(String message, Level level) {
        Logger.getLogger(InputOutput.class.getName()).log(level, message);
    }

    public static void line_error(String caller, String msg) {
        Logger logger = Logger.getLogger(InputOutput.class.getName());
        String finalmsg = String.format("[%s] %s", caller, msg);
        line_printf(finalmsg, Level.SEVERE);
        System.exit(1);
    }

    public static String mfilename(Object obj) {
        return obj.getClass().getEnclosingMethod().getName();
    }

    public static void main(String[] args) throws InterruptedException {
        // Example usage
        String name = mfilename(new Object() {
        });
        line_printf("This is an info message\nWith new line\n" + name);
        for (int i=0; i<70; i++ ){
            Thread.sleep(1000);
            line_warning(name, "This is an warning message\nWith new line");
        }
        line_error(name, "This is an error message\nWith new line");
    }
}
