package jline.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;

/**
 * An undirected graph used for routing and class-switching matrices.
 */
public class UndirectedGraph {
	
	private final int V;
	private final List<Set<Integer>> adj;	//undirected graph
	private final Set<Set<Integer>> wcc;
	private final Set<Integer> colsToIgnore;
	
	public UndirectedGraph(Matrix param, Set<Integer> colsToIgnore) {
		this.V = param.getNumCols();
		
		this.adj = new ArrayList<Set<Integer>>();
		for(int i = 0; i < V; i++)
			this.adj.add(new HashSet<Integer>());
		
		this.colsToIgnore = colsToIgnore;
		
		this.wcc = new TreeSet<Set<Integer>>(new Comparator<Set<Integer>>(){
			@Override
			public int compare(Set<Integer> arg0, Set<Integer> arg1) {
				Iterator<Integer> it0 = arg0.iterator();
				Iterator<Integer> it1 = arg1.iterator();
				while(it0.hasNext()) {
					if(it1.hasNext()) {
						int val0 = it0.next();
						int val1 = it1.next();
						if (val0 != val1)
							return val0 - val1;
					} else {
						return 1;
					}
				}
				
				if (it1.hasNext()) 
					return -1;
				else
					return 0;
			}
		});
		
		int[] col_idx = param.col_idx;
		int[] nz_rows = param.nz_rows;
		for(int colIdx = 0; colIdx < param.getNumCols(); colIdx++) {
			int col1 = col_idx[colIdx];
			int col2 = col_idx[colIdx+1];
			
			for(int i = col1; i < col2; i++) {
				int rowIdx = nz_rows[i];
				this.addEdge(rowIdx, colIdx);
				this.addEdge(colIdx, rowIdx);
			}
		}
	}
	
	public UndirectedGraph(Matrix param) {
		this(param, null);
	}
	
	public void addEdge(int s, int d) {
		if (s >= V || d >= V) 
			throw new RuntimeException("The index of row or column out of order");
		adj.get(s).add(d);
	}
	
	public void computeWeaklyConnectedComponents() {
		boolean[] visitedVertices = new boolean[V];
		
		for(int i = 0; i < V; i++) {
			if (!visitedVertices[i]) {
				Set<Integer> res = new TreeSet<Integer>();
				findCC(i, visitedVertices, res);
				if (colsToIgnore != null) {
					Iterator<Integer> it = colsToIgnore.iterator();
					while(it.hasNext()) {
						int num = it.next();
						res.remove(num);
					}
				}
				wcc.add(res);
			}
		}
	}
	
	public void findCC(int i, boolean[] visited, Set<Integer> chain) {
		visited[i] = true;
		chain.add(i);
		
		for(Integer v : adj.get(i)) {
			if (!visited[v])
				findCC(v, visited, chain);
		}
	}
	
	public Set<Set<Integer>> getWCC(){
		return this.wcc;
	}
}
