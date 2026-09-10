package jline.unified;

import java.util.ArrayList;
import java.util.List;

/**
 * Data class representing test results.
 */
public class TestResults {
    public int passed = 0;
    public int failed = 0;
    public int skipped = 0;
    public final List<String> errors = new ArrayList<String>();

    public TestResults() {}
}
