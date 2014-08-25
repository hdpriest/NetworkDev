/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package networkCalcPackage;

import com.apporiented.algorithm.clustering.AverageLinkageStrategy;
import com.apporiented.algorithm.clustering.Cluster;
import com.apporiented.algorithm.clustering.ClusteringAlgorithm;
import com.apporiented.algorithm.clustering.DefaultClusteringAlgorithm;
import com.apporiented.algorithm.clustering.visualization.DendrogramPanel;

/**
 *
 * @author Henry Priest, extending and adapting 
 * https://github.com/lbehnke/hierarchical-clustering-java.git
 */
public class Clustering {
            public static double[][] getDistances (GCNMatrix Similarities) {
            int Rows = Similarities.getNumRows();
            int Cols = Similarities.getNumColumns();
            double[][] distances = new double[Rows][Cols];
            for(int i=0;i<Rows;i++){
                for(int j=0;j<Cols;j++){
                    distances[i][j]=1 - Similarities.getValueByEntry(i,j);
                }
            }
            return distances;
        }
    
        public static void getClusters (GCNMatrix Similarities,String Method){
            ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
            // COULD run into serious problems with non-symmetric matricies
            double[][] distances = Clustering.getDistances(Similarities);
            String[] names = Similarities.getRowNames();
            Cluster cluster = alg.performClustering(distances,names,new AverageLinkageStrategy());
            DendrogramPanel dp = new DendrogramPanel();
            dp.setModel(cluster);
        }
    
}
