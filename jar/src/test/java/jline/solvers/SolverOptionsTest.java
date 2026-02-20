package jline.solvers;

import jline.lang.constant.SolverType;
import jline.VerboseLevel;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;
import static jline.TestTools.*;

/**
 * Comprehensive test to verify all SolverOptions parameters can be set and retrieved correctly.
 * This test ensures that the MATLAB parseCommonOptions function can configure all Java options.
 */
public class SolverOptionsTest {
    
    
    @Test
    public void testAllSolverOptionsParameters() {
        // Create a SolverOptions instance and set all parameters
        SolverOptions options = new SolverOptions(SolverType.MVA);
        
        // Basic options
        options.cache = false;
        options.cutoff = new Matrix(2, 3);
        options.cutoff.fill(100.0);
        options.force = true;
        options.hide_immediate = false;
        options.init_sol = new Matrix(3, 1);
        options.init_sol.fill(1.0);
        options.iter_max = 500;
        options.iter_tol = 1e-8;
        options.tol = 1e-6;
        options.keep = true;
        options.lang = "matlab";
        options.method = "exact";
        options.remote = false;
        options.remote_endpoint = "192.168.1.100";
        options.samples = 50000;
        options.seed = 12345;
        options.stiff = false;
        options.timespan = new double[]{0.0, 100.0};
        options.verbose = VerboseLevel.STD;
        
        // ODE solver options
        options.setODEMinStep(0.0001);
        options.setODEMaxStep(1.0);
        
        // Config options
        options.config.highvar = "interp";
        options.config.multiserver = "ps";
        options.config.np_priority = "fifo";
        
        List<Double> pstar = new ArrayList<>();
        pstar.add(1.5);
        pstar.add(2.0);
        pstar.add(2.5);
        options.config.pstar = pstar;
        
        options.config.fork_join = "enabled";
        options.config.merge = "aggressive";
        options.config.compress = "bdd";
        options.config.space_max = 1000000;
        options.config.interlocking = true;
        options.config.eventcache = false;
        options.config.hide_immediate = false;
        options.config.state_space_gen = "full";
        
        // Verify all parameters are set correctly
        assertFalse(options.cache);
        assertNotNull(options.cutoff);
        assertEquals(2, options.cutoff.getNumRows());
        assertEquals(3, options.cutoff.getNumCols());
        assertEquals(100.0, options.cutoff.get(0, 0), ZERO_TOL);

        assertTrue(options.force);
        assertFalse(options.hide_immediate);
        assertNotNull(options.init_sol);
        assertEquals(3, options.init_sol.getNumRows());
        assertEquals(1.0, options.init_sol.get(0, 0), ZERO_TOL);
        
        assertEquals(500, options.iter_max);
        assertEquals(1e-8, options.iter_tol, ZERO_TOL);
        assertEquals(1e-6, options.tol, ZERO_TOL);
        assertTrue(options.keep);
        assertEquals("matlab", options.lang);
        assertEquals("exact", options.method);
        assertFalse(options.remote);
        assertEquals("192.168.1.100", options.remote_endpoint);
        assertEquals(50000, options.samples);
        assertEquals(12345, options.seed);
        assertFalse(options.stiff);
        assertArrayEquals(new double[]{0.0, 100.0}, options.timespan, ZERO_TOL);
        assertEquals(VerboseLevel.STD, options.verbose);
        
        // Verify ODE solver options
        assertEquals(0.0001, options.odesolvers.odeminstep, ZERO_TOL);
        assertEquals(1.0, options.odesolvers.odemaxstep, ZERO_TOL);
        
        // Verify config options
        assertEquals("interp", options.config.highvar);
        assertEquals("ps", options.config.multiserver);
        assertEquals("fifo", options.config.np_priority);
        assertEquals(3, options.config.pstar.size());
        assertEquals(1.5, options.config.pstar.get(0), ZERO_TOL);
        assertEquals(2.0, options.config.pstar.get(1), ZERO_TOL);
        assertEquals(2.5, options.config.pstar.get(2), ZERO_TOL);
        assertEquals("enabled", options.config.fork_join);
        assertEquals("aggressive", options.config.merge);
        assertEquals("bdd", options.config.compress);
        assertEquals(1000000, options.config.space_max);
        assertTrue(options.config.interlocking);
        assertFalse(options.config.eventcache);
        assertFalse(options.config.hide_immediate);
        assertEquals("full", options.config.state_space_gen);
    }
    
    @Test
    public void testSolverOptionsDefaults() {
        // Test that default values are set correctly
        SolverOptions options = new SolverOptions();
        
        // Verify some key defaults
        assertTrue(options.cache);
        assertFalse(options.force);
        assertTrue(options.hide_immediate);
        assertEquals(100, options.iter_max);
        assertEquals(0.0001, options.iter_tol, ZERO_TOL);
        assertEquals(0.0001, options.tol, ZERO_TOL);
        assertTrue(options.keep);
        assertEquals("java", options.lang);
        assertEquals("default", options.method);
        assertFalse(options.remote);
        assertEquals("127.0.0.1", options.remote_endpoint);
        assertEquals(10000, options.samples);
        assertTrue(options.stiff);
        assertEquals(VerboseLevel.STD, options.verbose);
        
        // Verify config defaults
        assertEquals("default", options.config.highvar);
        assertEquals("default", options.config.multiserver);
        assertEquals("default", options.config.np_priority);
        assertEquals("default", options.config.fork_join);
        assertFalse(options.config.interlocking);
        assertTrue(options.config.eventcache);
        assertFalse(options.config.hide_immediate);
    }
    
    @Test
    public void testSolverOptionsCopy() {
        // Test that cloning works correctly for all parameters
        SolverOptions original = new SolverOptions(SolverType.FLUID);
        
        // Set some non-default values
        original.cache = false;
        original.iter_max = 999;
        original.tol = 1e-9;
        original.config.highvar = "test";
        original.config.space_max = 500000;
        
        SolverOptions cloned = original.copy();
        
        // Verify cloned values match
        assertEquals(original.cache, cloned.cache);
        assertEquals(original.iter_max, cloned.iter_max);
        assertEquals(original.tol, cloned.tol, ZERO_TOL);
        assertEquals(original.config.highvar, cloned.config.highvar);
        assertEquals(original.config.space_max, cloned.config.space_max);
        
        // Verify it's a deep copy by modifying original
        original.cache = true;
        original.config.highvar = "modified";
        
        // Cloned should retain original values
        assertFalse(cloned.cache);
        assertEquals("test", cloned.config.highvar);
    }
    
    @Test
    public void testCutoffParameterParsingWithSingleNumericArgument() {
        // Test parsing cutoff as single numeric argument (Integer)
        SolverOptions options1 = Solver.parseOptions(new SolverOptions(), 100);
        assertNotNull(options1.cutoff);
        assertEquals(100.0, options1.cutoff.get(0, 0), ZERO_TOL);
        
        // Test parsing cutoff as single numeric argument (Double)
        SolverOptions options2 = Solver.parseOptions(new SolverOptions(), 50.5);
        assertNotNull(options2.cutoff);
        assertEquals(50.5, options2.cutoff.get(0, 0), ZERO_TOL);
        
        // Test parsing cutoff with key-value pair (should still work)
        SolverOptions options3 = Solver.parseOptions(new SolverOptions(), "cutoff", 75);
        assertNotNull(options3.cutoff);
        assertEquals(75.0, options3.cutoff.get(0, 0), ZERO_TOL);
        
        // Test that it works with SolverOptions(SolverType) constructor too
        SolverOptions options4 = Solver.parseOptions(new SolverOptions(SolverType.CTMC), 200);
        assertNotNull(options4.cutoff);
        assertEquals(200.0, options4.cutoff.get(0, 0), ZERO_TOL);
    }
}
