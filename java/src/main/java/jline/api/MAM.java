package jline.api;

import jline.lang.JobClass;
import jline.lang.NetworkStruct;
import jline.lang.nodes.Station;
import jline.solvers.SolverOptions;
import jline.util.Matrix;

import java.util.*;

import static jline.lib.KPCToolbox.*;
import static jline.lib.M3A.*;


/**
* APIs for Matrix-Analytic Methods and Markovian Arrival Processes (MAPs).
*/
public class MAM {

	/**
	 * Compress a Marked Markovian Arrival processes.
	 * @param MMAP - representation of the MMAP
	 * @param config - configuration options
	 * @return compressed representation   
	 */	
	public static Map<Integer,Matrix> mmap_compress(Map<Integer,Matrix> MMAP, SolverOptions config){
		Map<Integer,Matrix> mmap = new HashMap<>();
		int K = MMAP.size()-2;
		if(Objects.equals(config.method, "default") ||Objects.equals(config.method, "mixture")||Objects.equals(config.method, "mixture.order1")){
			Matrix lambda = mmap_lambda(MMAP);
			Map<Integer,Map<Integer,Matrix>> AMAPs = mmap_maps(MMAP);
				for(int k=0;k<K;k++){
					AMAPs.put(k,mmpp2_fit1(map_mean(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_scv(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_skew(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_idc(AMAPs.get(k).get(0),AMAPs.get(k).get(1))));
					lambda.scale(1/lambda.elementSum());
					mmap = mmap_mixture(lambda,AMAPs);
				}

		}else if(Objects.equals(config.method,"mixture.order2")){
			//mmap = mmap_mixture_fit_mmap(MMAP).MMAP;
		}
		//TODO: more methods
		return mmap;
	}

	/**
	 * Compress a Marked Markovian Arrival processes as a mixture of 2-state acyclic Markovian Arrival Processes
	 * @param MMAP - representation of the MMAP
	 * @return compressed representation   
	 */
	public static Map<Integer,Matrix> mmap_compress(Map<Integer,Matrix> MMAP){
		Map<Integer,Matrix> mmap = new HashMap<>();
		int K = MMAP.size()-2;

		Matrix lambda = mmap_lambda(MMAP);
		Map<Integer,Map<Integer,Matrix>> AMAPs = mmap_maps(MMAP);
		for(int k=0;k<K;k++){
			AMAPs.put(k,mmpp2_fit1(map_mean(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_scv(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_skew(AMAPs.get(k).get(0),AMAPs.get(k).get(1)),map_idc(AMAPs.get(k).get(0),AMAPs.get(k).get(1))));
			lambda.scale(1/lambda.elementSum());
			mmap = mmap_mixture(lambda,AMAPs);
		}

		return mmap_normalize(mmap);
	}

	/**
	 * Re-indexes a collection of phase-type distributions (PHs)
	 * @param PHs - a collection of PHs indexed by station and class objects
	 * @return a collection of PHs indexed by station and class indexes
	 */
	public static Map<Integer,Map<Integer,Map<Integer,Matrix>>> ph_reindex(Map<Station,Map<JobClass,Map<Integer,Matrix>>> PHs, NetworkStruct sn){
		Map<Integer,Map<Integer,Map<Integer,Matrix>>> result = new HashMap<>();

		for (int i=0;i<sn.nstations;i++){
			for(int j=0; j<sn.nclasses;j++){
				if(j==0){
					result.put(i,new HashMap<>());
				}
				result.get(i).put(j,PHs.get(sn.stations.get(i)).get(sn.jobclasses.get(j)));
			}
		}
		return result;
	}

	/**
	 * Convert the M3A MMAP representation into the shorter BUTools representation.
	 * In M3A a MMAP is represented as MAP followed by D1 markings: D0,D1,D1a,D1b,..
	 * In BUTools a MMAP is represented as D0 followed by D1 markings: D0,D1a,D1b,..
	 * @param mmap - representation of the MMAP into M3A format
	 * @return MMAP representation in BUTools format
	 */
	public static Map<Integer,Matrix> mmap_shorten(Map<Integer,Matrix> mmap){
		Map<Integer,Matrix> result = new HashMap<>();
		for(int i=0;i<mmap.size();i++){
			if(i==0){
				result.put(0,mmap.get(0));
			}else if(i!=1){
				result.put(i-1,mmap.get(i));
			}
		}
		return result;
	}


}
