package networkCalcPackage;

import java.io.PrintWriter;

import net.sf.javaml.clustering.Clusterer;
import net.sf.javaml.clustering.IterativeKMeans;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.Instance;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.clustering.evaluation.ClusterEvaluation;
import net.sf.javaml.clustering.evaluation.SumOfAveragePairwiseSimilarities;
import net.sf.javaml.clustering.evaluation.SumOfSquaredErrors;

class Clustering {
	private Dataset DataSet = new DefaultDataset();
	private Dataset[] clusters;
	
	
	public Clustering (GCNMatrix Distances) {
		/* Once data is moved into the clustering stage, there is no recomputing;
		 * Load data instantly, cluster later
		 */
	    for (int i = 0; i < Distances.getNumRows() ; i++) {
	    	double[] Row = Distances.getRowByIndexDbl(i);
	    	String Locus = Distances.getRowName(i);
	        Instance this_Instance = new DenseInstance(Row,Locus); 
	        DataSet.add(this_Instance);
	        
	    }
	}
	
	public Dataset[] runKMeans (int k) {
		ClusterEvaluation ce = new SumOfAveragePairwiseSimilarities();
		ClusterEvaluation sse = new SumOfSquaredErrors();
		Clusterer ikm = new IterativeKMeans(1,k,ce);
		Clusterer km = new KMeans(k);
		Dataset[] clusters = km.cluster(DataSet);
		_clustersToFile(clusters,"test");
		return clusters;
	}
	private static void _clustersToFile (Dataset[] clusters, String path){
		try {	
				for(int c=0;c<clusters.length;c++){
					String nPath = "Test.cluster."+c+".txt";
					PrintWriter writer = new PrintWriter(nPath,"UTF-8");
					for(int I=0;I<clusters[c].size();I++){
						Instance this_Instance = clusters[c].get(I);
						String name = (String) this_Instance.classValue();
						writer.println(name);
					}
					writer.close();
				}
				
			} catch (Exception e){
				// 
			}
			
	}
}