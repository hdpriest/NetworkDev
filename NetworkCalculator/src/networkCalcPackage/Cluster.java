package networkCalcPackage;

import java.util.Iterator;
import java.util.List;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;


class Cluster {
        private GCNMatrix DISS;
        private float INF = Float.POSITIVE_INFINITY;
        private Dendrogram Dendrogram;
        private int Criteria;
        
        
	public Cluster (GCNMatrix Similarities,int Crit) {
            /*
            Adapted from Langfelder & Yau's adaptation of original Fortran code of Fionn Murtagh
            http://cran.r-project.org/web/packages/flashClust/index.html
            http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/flashClust/

            Published here:
            http://www.jstatsoft.org/v46/i11/paper

            No changes, just a java implementation with perhaps some multithreading support ... we'll see.
            */
            /*
            Init variables, get Distances, etc...
            */

            
            int N = Similarities.getNumRows();
            Criteria = Crit;
            DISS = new GCNMatrix(N,N);
            for (int i = 0; i < N ; i++) {
                float[] Row = Similarities.getRowByIndexAsDistance(i);
                DISS.addRow(Row);
	    }
            Dendrogram = new Dendrogram(DISS);
            Dendrogram.getDendrogram(Criteria);
	}
        
        public void dynamicTreeCut (int MinSize) {
            // implementation of dynamic tree cut and flashclust from horvath and langfelder et al
            // This will be pretty ham-handed.
            //float cutoff = 0.99f * (1.0f-(CRIT[1]));
            float[] Dist = Dendrogram.getHeights();
            int N = Dendrogram.getNumberOfBranches();
            float cutoff = 0.99f * Dist[N-2];
            System.err.println("Cutoff is 99% of dendrogram height: " + cutoff);
            int[][] Clusters = Dendrogram.staticCut(cutoff,MinSize);
         
            for(int i=0;i<N;i++){
                if(Clusters[i] == null) continue;
                int thisN = Clusters[i].length;
                if(thisN < MinSize) continue;
                GCNMatrix thisDist = new GCNMatrix(thisN,thisN);
                System.out.println("Main Cluster size: " + thisN);
                for(int j=0;j<thisN;j++){
                    int Ind_i = Clusters[i][j];
                    float[] Row = new float[thisN];
                    for(int J=0;J<thisN;J++){
                        int Ind_j = Clusters[i][J];
                        float d = DISS.getValueByEntry(Ind_i,Ind_j);
                        Row[J]=d;
                    }
                    thisDist.addRow(Row);        
                }
                _adaptiveTreeCut(thisDist);
            }
            System.exit(0);
        }
        
        private void _adaptiveTreeCut(GCNMatrix thisDist) {
            Dendrogram this_Dendro = new Dendrogram(thisDist);
            this_Dendro.getDendrogram(Criteria);
// thisDendro now has the unique merge order and heights of the rows of thisDist
            int[] this_IA = this_Dendro.getMergeRoots();
            int[] this_IB = this_Dendro.getMergeLeaves();
            float[] this_Dist  = this_Dendro.getHeights();
            float L=_getMean(this_Dist);
            int[] this_Order = this_Dendro.getDendroOrder();
            /*    
            for (int Z =0;Z<this_IB.length;Z++){
               System.out.println(this_IA[Z]+ "\t" + this_IB[Z] + "\t" + this_Dist[Z] + "\t" + this_Order[Z]);
            }
            System.exit(0);
            */
            int[][] clusters = _treeCutCore(this_Order,this_Dist,L,25);
            for(int c=0;c<clusters.length;c++){
                System.out.println("\t Sub Cluster size: " + clusters[c].length);
            }
        }
        
        private float[] _getS (int[] this_Order, float L,float[] this_Dist) { 
            float[] S = new float[this_Order.length];
            for(int o=0;o<this_Order.length;o++){
            //for(int h=0;h<H.length;h++){
                S[o] = this_Dist[this_Order[o]]-L; // go get the heights for the branches in O, and calculate S based on those
                System.out.println(o + "\t" + this_Order[o] +"\t" + S[o] + "\t" + L);
            }
            return S;
        }
        
        private int[][] _treeCutCore (int[] this_Order, float[] this_Dist, float L_o, int T) {
            boolean[] TP = new boolean[this_Order.length];
            float[] S = _getS(this_Order,L_o,this_Dist);
            int last=0;
            ArrayList<Integer[]> Clusters = new ArrayList<Integer[]>();
            for(int s=0;s<S.length-1;s++){
                TP[s] = (S[s] * S[s+1] <= 0.0f);
                if(TP[s] == true){ 
                    int C = s-last+1;
/*
I think the correct thing to do here is to implement the multiple height cutter
(AdaptiveTreeCut) and remove the C<T statement. Then, implement the merger
that is based on heights. If a height doesn't return any cluters (because
they're all too small), then a higher or lower cutoff might do it correctly.
That, combined with merging the small modules into larger modules
might be the path forward.
*/                        
                    Integer[] cluster = {last,C};
                    Clusters.add(cluster);
                    last = s;
                }
            }
            // Build ordered set of clusters // AND filter at the same time
            /*
            int c = 0;
            while(c < Clusters.size()-1){
                Integer[] cluster = Clusters.get(c);
                if(cluster.length < T) continue;
                if(Clusters.get(c+1).length < T){
                    float meanMain = 
                }
                    
                
                c++;
            }
                 */   
            Iterator<Integer[]> I = Clusters.iterator();
            int[][] Cs = new int[Clusters.size()][];
            int m = 0;
            while(I.hasNext()){
                Integer[] Cluster = I.next();
                int breakpoint = Cluster[0];
                int forwardrun = Cluster[1];
                int[] cluster = new int[forwardrun];
                for(int i = breakpoint;i<forwardrun;i++){
                    int ind = i-breakpoint;
                    cluster[ind] = i;
                }
                Cs[m]=cluster;
                m++;
            }
            return Cs;
               // Either clean-up here, or not. Probably best here.
        }
        private float _getMean (float[] Distances){
            float mean =0.0f;
            for(int i=0;i<Distances.length;i++){
                mean = mean + Distances[i];
            }
            mean = mean / Distances.length;
            return mean;
        }
        
  
	private float[] _subsetArray (int[] index, float[] TA){
            float[] subset = new float[index.length];
            for(int i=0;i<index.length;i++){
                subset[i] = TA[index[i]];
            }
            return subset;
        }
        
        private int[] _subsetArray (int[] index, int[] TA){
            int[] subset = new int[index.length];
            for(int i=0;i<index.length;i++){
                subset[i] = TA[index[i]];
            }
            return subset;
        }
        
	private static void _clustersToFile (int[] Cluster, int M){
            try {	
		//Iterator<Integer> Node = Cluster.iterator();
                String nPath = "Cluster." + M + ".txt"; 
                PrintWriter writer = new PrintWriter(nPath,"UTF-8");
                //while(Node.hasNext()){
                //	int node = Node.next();
                for(int i=0;i<Cluster.length;i++){
                	int node = Cluster[i];
                	writer.println(node);
                }
                writer.close();	
            } catch (Exception e){
				// 
            }
			
	}
     
        
}