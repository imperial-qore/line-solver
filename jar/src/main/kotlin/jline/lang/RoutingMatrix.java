/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.lang;

import jline.GlobalConstants;
import jline.lang.constant.RoutingStrategy;
import jline.lang.nodes.Cache;
import jline.lang.nodes.ClassSwitch;
import jline.lang.nodes.Node;
import jline.lang.sections.ClassSwitcher;
import jline.util.matrix.Matrix;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class representing a probabilistic routing matrix
 */
public class RoutingMatrix implements Serializable {
    private final List<JobClass> jobClasses;
    private final List<Node> nodes;
    //In order not to modify most part of the code, we will use list. But when set it to NetworkStruct, it will be transferred to map.
    private final List<List<Matrix>> routings;
    private final Matrix csMatrix;
    private final Network model;
    private Map<JobClass, Integer> classIndexMap;
    private Map<Node, Integer> nodeIndexMap;
    private boolean hasUnappliedConnections;
    private boolean hasClassSwitches;

    public RoutingMatrix() {
        this.jobClasses = new ArrayList<JobClass>();
        this.nodes = new ArrayList<Node>();
        this.hasUnappliedConnections = false;

        int I = this.nodes.size();
        int K = this.jobClasses.size();
        this.csMatrix = new Matrix(K, K, K * K);
        for (int i = 0; i < K; i++)
            this.csMatrix.set(i, i, 1.0);
        this.routings = new ArrayList<List<Matrix>>();
        this.model = new Network("");
        this.hasClassSwitches = false;
    }

    public RoutingMatrix(Network model, List<JobClass> jobClasses, List<Node> nodes) {
        int nJobClasses = jobClasses.size();
        int nNodes = nodes.size();
        this.jobClasses = new ArrayList<JobClass>(jobClasses);
        this.nodes = new ArrayList<Node>(nodes);
        this.classIndexMap = new HashMap<JobClass, Integer>();
        this.nodeIndexMap = new HashMap<Node, Integer>();
        this.hasUnappliedConnections = false;

        for (int j = 0; j < this.nodes.size(); j++) {
            this.nodeIndexMap.put(this.nodes.get(j), j);
        }

        this.model = model;
        this.csMatrix = new Matrix(nJobClasses, nJobClasses, nJobClasses * nJobClasses);
        for (int i = 0; i < nJobClasses; i++)
            this.csMatrix.set(i, i, 1.0);
        this.hasClassSwitches = false;

        routings = new ArrayList<List<Matrix>>(nJobClasses);
        for (int i = 0; i < nJobClasses; i++) {
            List<Matrix> frame = new ArrayList<Matrix>(nJobClasses);
            for (int j = 0; j < nJobClasses; j++)
                frame.add(generateEmptyNodeOrClassRouting(nNodes));

            this.routings.add(frame);
            this.classIndexMap.put(this.jobClasses.get(i), i);
        }
    }

    public void addClass(JobClass jobClass) {
        if (this.jobClasses.contains(jobClass)) {
            // idempotent
            return;
        }

        int classIdx = this.jobClasses.size();
        this.jobClasses.add(jobClass);
        this.classIndexMap.put(jobClass, classIdx);

        int nJobClasses = jobClasses.size();
        int nNodes = nodes.size();

        this.csMatrix.expandMatrix(nJobClasses, nJobClasses, nJobClasses * nJobClasses);
        this.csMatrix.set(nJobClasses - 1, nJobClasses - 1, 1);

        List<Matrix> frame = new ArrayList<Matrix>();
        for (int i = 0; i < nJobClasses - 1; i++) {
            this.routings.get(i).add(this.generateEmptyNodeOrClassRouting(nNodes)); // Old class to the new class
            frame.add(this.generateEmptyNodeOrClassRouting(nNodes)); // New class to old class
        }

        frame.add(this.generateEmptyNodeOrClassRouting(nNodes));  // New class to New class
        routings.add(frame);
    }

    public void addConnection(Node sourceNode, Node destNode) {
        for (JobClass jobClass : this.jobClasses) {
            this.addConnection(sourceNode, destNode, jobClass, jobClass);
        }
    }

    public void addConnection(Node sourceNode, Node destNode, JobClass jobClass) {
        if (sourceNode.getRoutingStrategy(jobClass) == RoutingStrategy.DISABLED) {
            return;
        }

        this.hasUnappliedConnections = true;
        this.addConnection(jobClass, jobClass, sourceNode, destNode, Double.NaN);
    }

    public void addConnection(Node sourceNode, Node destNode, double probability) {
        for (JobClass jobClass : this.jobClasses) {
            this.addConnection(sourceNode, destNode, jobClass, probability);
        }
    }

    public void addConnection(Node sourceNode, Node destNode, JobClass originClass, JobClass targetClass) {
        if (sourceNode.getRoutingStrategy(originClass) == RoutingStrategy.DISABLED) {
            return;
        }

        this.hasUnappliedConnections = true;
        this.addConnection(originClass, targetClass, sourceNode, destNode, Double.NaN);
    }

    public void addConnection(Node sourceNode, Node destNode, JobClass jobClass, double probability) {
        // Only skip DISABLED routing when probability is not explicitly specified (NaN).
        // When probability is explicit (e.g., from serialRouting with 1.0), the caller is
        // intentionally overriding the default DISABLED state set by Dispatcher initialization.
        if (Double.isNaN(probability) && sourceNode.getRoutingStrategy(jobClass) == RoutingStrategy.DISABLED) {
            return;
        }

        if (Double.isNaN(probability)) {
            this.hasUnappliedConnections = true;
        }

        this.addConnection(jobClass, jobClass, sourceNode, destNode, probability);
    }

    public void addConnection(Node sourceNode, Node destNode, JobClass originClass, JobClass targetClass, double probability) {
        addConnection(originClass, targetClass, sourceNode, destNode, probability);
    }

    public void addConnection(JobClass originClass, JobClass targetClass, Node sourceNode, Node destNode, double probability) {

        int originClassIdx = getClassIndex(originClass);
        int targetClassIdx = getClassIndex(targetClass);
        int sourceNodeIdx = getNodeIndex(sourceNode);
        int destNodeIdx = getNodeIndex(destNode);

        if (probability <= GlobalConstants.Zero) {
            // Still need to set to 0 to clear any existing routing
            routings.get(originClassIdx).get(targetClassIdx).unsafeSet(sourceNodeIdx, destNodeIdx, 0.0);
            return;
        }

        if (!originClass.equals(targetClass) || sourceNode instanceof Cache || destNode instanceof Cache) {
            this.hasClassSwitches = true;
            csMatrix.set(originClassIdx, targetClassIdx, 1);
        }

        routings.get(originClassIdx).get(targetClassIdx).unsafeSet(sourceNodeIdx, destNodeIdx, probability);

        if (sourceNode instanceof ClassSwitch) {
            ClassSwitcher server = ((ClassSwitcher) sourceNode.getServer());
            for (int r = 0; r < this.jobClasses.size(); r++) {
                for (int s = 0; s < this.jobClasses.size(); s++) {
                    if (server.applyCsFun(r, s) > 0)
                        csMatrix.set(r, s, 1.0);
                }
            }
        }
    }

    public void addNode(Node node) {
        if (this.nodes.contains(node)) {
            return;
        }

        int nodeIdx = this.nodes.size();

        this.nodes.add(node);
        this.nodeIndexMap.put(node, nodeIdx);

        int I = this.nodes.size();

        for (List<Matrix> classArray : this.routings) {
            for (int i = 0; i < classArray.size(); i++)
                classArray.get(i).expandMatrix(I, I, I * I);
        }
    }

    private Matrix generateEmptyNodeOrClassRouting(int size) {
        return new Matrix(size, size, size * size);
    }

    public Matrix get(JobClass jobclass1, JobClass jobclass2) {
        return routings.get(jobclass1.getIndex() - 1).get(jobclass2.getIndex() - 1);
    }

    public double get(JobClass jobclass1, JobClass jobclass2, Node node1, Node node2) {
        return routings.get(jobclass1.getIndex() - 1).get(jobclass2.getIndex() - 1).get(node1.getNodeIndex(), node2.getNodeIndex());
    }

    public Matrix get(int jobclass1, int jobclass2) {
        return routings.get(jobclass1 - 1).get(jobclass2 - 1);
    }

    public double get(int jobclass1, int jobclass2, int node1, int node2) {
        return routings.get(jobclass1 - 1).get(jobclass2 - 1).get(node1 - 1, node2 - 1);
    }

    public JobClass getClass(int classIndex) {
        return this.jobClasses.get(classIndex);
    }

    public int getClassIndex(JobClass jobClass) {
        return this.classIndexMap.get(jobClass);
    }

    public int getNodeIndex(Node node) {
        return this.nodeIndexMap.get(node);
    }

    public void print() {
        for (int r = 0; r < this.jobClasses.size(); r++) {
            for (int s = 0; s < this.jobClasses.size(); s++) {
                System.out.println(System.out.format("Class %d -> Class %d:", r, s));
                this.routings.get(r).get(s).print();
            }
        }
    }

    public void resolveClassSwitches() {
        //line 163 - 168
        int nNodes = nodes.size();
        int nClasses = jobClasses.size();
        List<List<Matrix>> csnodematrix = new ArrayList<List<Matrix>>(nNodes);
        for (int i = 0; i < nNodes; i++) {
            List<Matrix> classMatrix = new ArrayList<Matrix>(nNodes);
            for (int j = 0; j < nNodes; j++)
                classMatrix.add(generateEmptyNodeOrClassRouting(nClasses));
            csnodematrix.add(classMatrix);
        }

        //line 170 - 179
        for (int row = 0; row < nClasses; row++) {
            for (int col = 0; col < nClasses; col++) {
                Matrix nodeRouting = routings.get(row).get(col);
                if (nodeRouting.getNonZeroLength() > 0) {
                    int[] col_idx = nodeRouting.getColIndexes();
                    int[] nz_rows = nodeRouting.getNonZeroRows();
                    double[] nz_values = nodeRouting.getNonZeroValues();

                    for (int colIdx = 0; colIdx < nodeRouting.getNumCols(); colIdx++) {
                        int col1 = col_idx[colIdx];
                        int col2 = col_idx[colIdx + 1];

                        for (int i = col1; i < col2; i++) {
                            int rowIdx = nz_rows[i];
                            double values = nz_values[i];
                            csnodematrix.get(rowIdx).get(colIdx).set(row, col, values);
                        }
                    }
                }
            }
        }

        //line 196 - 207
        for (int i = 0; i < nNodes; i++) {
            for (int j = 0; j < nNodes; j++) {
                Matrix classRouting = csnodematrix.get(i).get(j);
                Matrix res = classRouting.sumRows();
                classRouting.divideRows(res.getNonZeroValues(), 0);
                for (int r = 0; r < nClasses; r++) {
                    if (res.get(r) == 0)
                        classRouting.set(r, r, 1.0);
                }
            }
        }

        //line 209 - 220
        int[][] csid = new int[nNodes][nNodes];
        for (int i = 0; i < nNodes; i++) {
            for (int j = 0; j < nNodes; j++) {
                Matrix classRouting = csnodematrix.get(i).get(j);
                if (!classRouting.isDiag()) {
                    String csname = "CS_" + nodes.get(i).getName() + "_to_" + nodes.get(j).getName();
                    ClassSwitch csnode = new ClassSwitch(model, csname, classRouting);
                    csnode.autoAdded = true;
                    this.addNode(csnode); //line 239 - 243
                    csid[i][j] = model.getNumberOfNodes() - 1;
                }
            }
        }

        // lines 222-233
        for (int i = 0; i < nNodes; i++) {
            // This is to ensure that also stateful cs like caches are accounted
            if (this.nodes.get(i) instanceof Cache) {
                Cache cache = (Cache) this.nodes.get(i);
                for (int r = 0; r < cache.getHitClass().getNumCols(); r++) {
                    if (cache.getHitClass().get(r) != -1) {
                        csMatrix.set(r, (int) cache.getHitClass().get(r), 1);
                    }
                }
                for (int r = 0; r < cache.getMissClass().getNumCols(); r++) {
                    if (cache.getMissClass().get(r) != -1) {
                        csMatrix.set(r, (int) cache.getMissClass().get(r), 1);
                    }
                }
            }
        }

        //line 245 - 260
        for (int i = 0; i < nNodes; i++) {
            for (int j = 0; j < nNodes; j++) {
                if (csid[i][j] > 0) {
                    for (int r = 0; r < nClasses; r++) {
                        for (int s = 0; s < nClasses; s++) {
                            Matrix nodeRouting = routings.get(r).get(s);
                            if (nodeRouting.get(i, j) > 0) {
                                Matrix from = routings.get(r).get(r);
                                from.set(i, csid[i][j], from.get(i, csid[i][j]) + nodeRouting.get(i, j));
                                nodeRouting.remove(i, j);
                                Matrix to = routings.get(s).get(s);
                                to.set(csid[i][j], j, 1.0);
                            }
                        }
                    }
                }
            }
        }
        this.hasClassSwitches = false;
    }

    public void resolveUnappliedConnections() {

        int I = nodes.size();
        for (List<Matrix> jobClassRouting : this.routings) {
            for (Matrix nodeRouting : jobClassRouting) {
                for (int row = 0; row < I; row++) {
                    double residProb = 1;
                    int nUnapplied = 0;
                    for (int col = 0; col < I; col++) {
                        Double routingAmount = nodeRouting.get(row, col);
                        if (routingAmount.isNaN()) {
                            nUnapplied++;
                        } else {
                            residProb -= routingAmount;
                        }
                    }
                    if (nUnapplied == 0) {
                        continue;
                    }
                    double unitProb = residProb / nUnapplied;
                    for (int col = 0; col < I; col++) {
                        if (Double.isNaN(nodeRouting.get(row, col)))
                            nodeRouting.set(row, col, unitProb);
                    }
                }
            }
        }
        this.hasUnappliedConnections = false;
    }

    private Map<JobClass, Map<JobClass, Matrix>> routingListToMap() {
        Map<JobClass, Map<JobClass, Matrix>> routingMap = new HashMap<JobClass, Map<JobClass, Matrix>>();
        for (int i = 0; i < routings.size(); i++) {
            JobClass jobClass = jobClasses.get(i);
            List<Matrix> routingList = routings.get(i);
            Map<JobClass, Matrix> map = new HashMap<JobClass, Matrix>();
            for (int j = 0; j < routingList.size(); j++) {
                map.put(jobClasses.get(j), routingList.get(j).copy());
            }
            routingMap.put(jobClass, map);
        }
        return routingMap;
    }

    private Map<JobClass, Map<JobClass, Double>> routingMatrixToMap(Matrix classRouting) {
        Map<JobClass, Map<JobClass, Double>> csm = new HashMap<JobClass, Map<JobClass, Double>>();
        int[] col_idx = classRouting.getColIndexes();
        int[] nz_rows = classRouting.getNonZeroRows();
        double[] nz_values = classRouting.getNonZeroValues();
        for (int colIdx = 0; colIdx < classRouting.getNumCols(); colIdx++) {
            int col1 = col_idx[colIdx];
            int col2 = col_idx[colIdx + 1];

            for (int m = col1; m < col2; m++) {
                int rowIdx = nz_rows[m];
                double values = nz_values[m];
                Map<JobClass, Double> map = csm.getOrDefault(jobClasses.get(rowIdx), new HashMap<JobClass, Double>());
                map.put(jobClasses.get(colIdx), values);
                csm.put(jobClasses.get(rowIdx), map);
            }
        }
        return csm;
    }

    public void set(JobClass jobclass1, JobClass jobclass2, Node srcNode, Node destNode, double probability) {
        this.addConnection(jobclass1, jobclass2, srcNode, destNode, probability);
    }

    public void set(JobClass jobclass1, Node srcNode, Node destNode, double probability) {
        this.addConnection(srcNode, destNode, jobclass1, probability);
    }

    public void set(Node srcNode, Node destNode, double probability) {
        this.addConnection(srcNode, destNode, probability);
    }

    public void set(JobClass jobclass, Matrix rt) {
        for (int i = 0; i < this.nodes.size(); i++)
            for (int k = 0; k < this.nodes.size(); k++)
                set(jobclass, this.nodes.get(i), this.nodes.get(k), rt.get(i, k));
    }

    public void set(JobClass jobclass1, JobClass jobclass2, Matrix rt) {
        for (int i = 0; i < this.nodes.size(); i++)
            for (int k = 0; k < this.nodes.size(); k++)
                set(jobclass1, jobclass2, this.nodes.get(i), this.nodes.get(k), rt.get(i, k));
    }

    public void set(JobClass jobclass, RoutingMatrix rt) {
        this.set(jobclass, jobclass, rt.get(jobclass, jobclass));
    }

    public void set(JobClass jobclass1, JobClass jobclass2, RoutingMatrix rt) {
        this.set(jobclass1, jobclass2, rt.get(jobclass1, jobclass2));
    }
    
    public void set(int classIndex1, int classIndex2, Matrix rt) {
        if (classIndex1 >= 0 && classIndex1 < this.jobClasses.size() && 
            classIndex2 >= 0 && classIndex2 < this.jobClasses.size()) {
            
            JobClass localFromClass = this.jobClasses.get(classIndex1);
            JobClass localToClass = this.jobClasses.get(classIndex2);
            
            for (int i = 0; i < this.nodes.size(); i++) {
                for (int k = 0; k < this.nodes.size(); k++) {
                    set(localFromClass, localToClass, this.nodes.get(i), this.nodes.get(k), rt.get(i, k));
                }
            }
        }
    }

    public void set(Node srcNode, Node destNode) {
        this.addConnection(srcNode, destNode);
    }

    public void setRouting(Network model) {
        if (this.hasUnappliedConnections) {
            this.resolveUnappliedConnections();
        }
        Map<JobClass, Map<JobClass, Matrix>> rtorig = routingListToMap();

        if (this.hasClassSwitches) {
            this.resolveClassSwitches();
        }

        //line 262-273
        for (int r = 0; r < this.jobClasses.size(); r++) {
            Matrix routing = routings.get(r).get(r);

            int[] col_idx = routing.getColIndexes();
            int[] nz_rows = routing.getNonZeroRows();
            double[] nz_values = routing.getNonZeroValues();

            for (int colIdx = 0; colIdx < routing.getNumCols(); colIdx++) {
                int col1 = col_idx[colIdx];
                int col2 = col_idx[colIdx + 1];

                for (int i = col1; i < col2; i++) {
                    int rowIdx = nz_rows[i];
                    double value = nz_values[i];
                    if (value > GlobalConstants.Zero) {
                        model.addLink(rowIdx, colIdx);
                        this.nodes.get(rowIdx).setRouting(this.jobClasses.get(r), RoutingStrategy.PROB, this.nodes.get(colIdx), value);
                    }
                }
            }

        }
        // Update rtorig without replacing the entire NetworkStruct.
        // This preserves computed fields (rates, isstatedep, etc.) during relink.
        model.updateRtorig(rtorig);
        model.setCsMatrix(this.csMatrix);
        if (this.model != model) {
            this.model.updateRtorig(rtorig);
            this.model.setCsMatrix(this.csMatrix);
        }
    }

    /*
     * ===================================================================================
     * MISSING METHODS FROM MATLAB ROUTINGMATRIX IMPLEMENTATION - NOT YET MIGRATED
     * ===================================================================================
     *
     * Based on analysis of /matlab/src/lang/RoutingMatrix.m
     */

    // =================== MATLAB OPERATOR OVERLOADING METHODS ===================
    // public Object subsref(SubscriptReference s)  // Custom subscripted reference behavior
    // public RoutingMatrix subsasgn(SubscriptReference s, Object... values)  // Custom subscripted assignment behavior

    // =================== ENHANCED SETTING METHODS ===================
    // public void set(int jobclass1, Object mat)  // Set routing matrix for single job class
    // public void set(int jobclass1, int jobclass2, Matrix mat)  // Set routing matrix between job classes
    // public void set(int jobclass1, int jobclass2, int node1, int node2, double val)  // Set single routing probability

    // =================== ENHANCED ACCESS METHODS ===================
    // public Matrix[][] getCell()  // Get routing matrix as cell array (equivalent to MATLAB cell)

    // =================== STATIC UTILITY METHODS ===================
    // public static RoutingMatrixReturn rtnodes2rtorig(NetworkStruct sn)  // Convert node routing to original format
    // 
    // Helper class for static method return:
    // public static class RoutingMatrixReturn {
    //     public Matrix[][][] rtorigcell;  // Original routing cell array
    //     public Matrix rtorig;            // Original routing matrix
    // }

}
