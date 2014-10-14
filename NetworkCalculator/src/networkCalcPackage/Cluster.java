package networkCalcPackage;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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
            int[] Order = Dendrogram.getDendroOrder();
            float cutoff = 0.99f * Dist[N-2];
            System.err.println("Cutoff is 99% of dendrogram height: " + cutoff);
            
            ArrayList<int[]> Clusters = Dendrogram.staticCut(cutoff,MinSize);

            int c = 0;
            while(c < Clusters.size()){
                int[] this_cluster = Clusters.get(c); 
                /*
                 * The content of this_cluster is important it is a set of branch IDs which are 
                 * not necessarily contiguous
                 * not necessarily origin-0
                 * not necessarily sorted
                 * this_cluster[i]=branch_id
                 * Dist[branch_id]=height
                 */
                int this_N = this_cluster.length;
                if(Clusters.size() > N){ // checks because Clusters changes many times
                    System.err.println("More clusters than genes. Something has gone horribly wrong.\n\n");
                    System.exit(0);
                }
                //System.out.println("Processing cluster " + c + ", of size: " + this_N);
                if(this_N<MinSize){
                    c++;
                    continue;
                }

                ArrayList<int[]> new_Clusters = _adaptiveTreeCut(this_cluster);
/*
 * The content of new_Clusters is a set of new branch IDs which:
 * contain all integers from i=0 to i<this_cluster.length
 * are origin-0
 * not necessarily sorted
 * no missing values
 * These branch IDs refer to the IDs of this_cluster.                
 */
                if(new_Clusters.size() == 0){
                        System.err.println("ZERO new clusters... whooops\n\n");
                        System.exit(0);
                }
                if(new_Clusters.size() == 1){
                // if new_Clusters size ==1 (i.e., no new clusters here) then.... proceed to next cluster.
                    System.out.println("Found no new clusters");
                    c++;
                    continue;
                }else{
                // if new_Clusters size ne 0 && ne 1, insert new clusters starting at c, and do not iterate (i.e., start @ c again)
                    int n = 0;
                    Clusters.remove(c);
                    while(n < new_Clusters.size()){
                        int nc = c+n; // starts @ c
                        System.out.println("Inserting cluster at " + nc);
                        Clusters.add(nc,new_Clusters.get(n));
                        n++;
                    }
                    // no iteration
                    continue;// continues anyway; this just makes you look like a tool;
                }
            }
            System.out.println("Done.");
            System.out.println("Obtained " + Clusters.size() + " Clusters.");
            Iterator<int[]> it = Clusters.iterator();
            while(it.hasNext()){
                int[] cluster = it.next();
                System.out.println("Final cluster size: " + cluster.length);
            }
/*
            }
            for(int i=0;i<N;i++){
                if(Clusters[i] == null) continue;
                int thisN = Clusters[i].length;
                if(thisN < MinSize) continue;
                System.out.println("Main Cluster size: " + thisN);
                _adaptiveTreeCut(Clusters[i]);
            }
*/
            System.exit(0);
        }
        
        private GCNMatrix _getDistForCluster (int[] cluster){
            int thisN = cluster.length;
            GCNMatrix thisDist = new GCNMatrix(thisN,thisN);
            for(int j=0;j<thisN;j++){
                int Ind_i = cluster[j]; // branch id at j
                float[] Row = new float[thisN];
                for(int J=0;J<thisN;J++){
                    int Ind_j = cluster[J]; // branch id at J
                    float d = DISS.getValueByEntry(Ind_i,Ind_j); 
                    Row[J] = d;
                }
                thisDist.addRow(Row);
            }
            return thisDist;
        }
        
        private float _getMax (float[] array){
            float max=0.0f;
            for(int i=0;i<array.length;i++){
                if(array[i] > max) max = array[i];
            }
            return max;
        }
        
        private float _getMin (float[] array){
            float min=INF;
            for(int i=0;i<array.length;i++){
                if(array[i] < min) min = array[i];
            }
            return min;
        }
        
        private ArrayList<int[]> _adaptiveTreeCut(int[] Cluster) {
            GCNMatrix thisDist = _getDistForCluster(Cluster);
            Dendrogram this_Dendro = new Dendrogram(thisDist);
            this_Dendro.getDendrogram(Criteria);
// thisDendro now has the unique merge order and heights of the rows of thisDist
            float[] this_Dist  = this_Dendro.getHeights();
            int[] this_Order = this_Dendro.getDendroOrder();
/*
but now all of the branch IDs change!            
here is the problem. As you descend levels of this algorithm:
1) static cut
2) dynamic tree cut (per cluster)
3) adaptive tree cut (new dendro)
4) tree cut core (clusters)

The branch IDs of 4 don't travel back to 1... this means that at the final step, when you re-ass
all of the cluster branches to the main cluster tree, everything is FUBAR'd (doesn't conjugate like that...)

Solution: this_Order (above), holds the ordering of the branches of this_Dendro
*/            
            /*
             * The content of Cluster is important it is a set of branch IDs which are 
             * not necessarily contiguous
             * not necessarily origin-0
             * not necessarily sorted
             * Cluster[i]=branch_id
             * Dist[branch_id]=height
             */
            float L_naught=_getMean(this_Dist);
            float L_max   = 0.5f * (L_naught + _getMax(this_Dist));
            float L_min   = 0.5f * (L_naught + _getMin(this_Dist));
            // call l_naught, then l_min, then l_max, 
            //int[][] clusters = _treeCutCore(this_Order,this_Dist,L_naught,25);
            ArrayList<int[]> clusters = new ArrayList<int[]>();
            if(Cluster.length < 10){
            	// can't cluster things that are too small.
            	clusters.add(Cluster);
            	return clusters;
            }
            System.out.println("working on cluster of size "+ Cluster.length);
            clusters =	_treeCutCore(this_Order,this_Dist,L_naught,25);
            System.out.println("Found " + clusters.size() + " clusters on this iteration (naught) " + L_naught);
            //if(Arrays.equals(Cluster,clusters[0])){
            if(clusters.size() <= 1){
                clusters = _treeCutCore(this_Order,this_Dist,L_min,25);
                System.out.println("Found " + clusters.size() + " clusters on this iteration (lower) " + L_min );
            }
            if(clusters.size() <= 1){
                clusters = _treeCutCore(this_Order,this_Dist,L_max,25);
                System.out.println("Found " + clusters.size() + " clusters on this iteration (upper) " + L_max);
            }
            
            if(clusters.size() == 0){ // No new clusters are found at any cut height
            	clusters.add(Cluster); // Send back the original.
            	return clusters;
            }else if(clusters.size() == 1){
            	clusters.remove(0); // treeCutCore can drop nodes from clusters and return only one.
            	// this might not be correct. will have to look.
            	clusters.add(Cluster);
            	return clusters;
            }else{
            	for(int c=0;c<clusters.size();c++){
            		System.out.println("\t Sub Cluster size: " + clusters.get(c).length);
            		int[] cluster = clusters.get(c);
            		for(int i=0;i<cluster.length;i++){
            			cluster[i]=Cluster[cluster[i]]; // map back to original branch IDs
            		}
            		clusters.set(c, cluster); // maybe a more elegant way to do this.
            		// will ONLY work if treeCutCore doesn't further mangle IDs.
            	}
            }
            return clusters;
        }
        
        private float[] _getS (int[] this_Order, float L,float[] this_Dist) { 
            float[] S = new float[this_Order.length];
            for(int o=0;o<this_Order.length;o++){
                S[o] = this_Dist[this_Order[o]]-L; // go get the heights for the branches in O, and calculate S based on those
                //System.out.println(o + "\t" + this_Order[o] +"\t" + S[o] + "\t" + L);
            }
            return S;
        }
        private float _getMeanForCluster (float[] this_Dist, int[] this_Order, int[] this_Cluster){
            float S = 0.0f;
            for(int c=0;c<this_Cluster.length;c++){
                S = S + this_Dist[this_Order[this_Cluster[c]]];
            }
            S = S / this_Cluster.length;
            return S;
        }
        private ArrayList<int[]> _treeCutCore (int[] this_Order, float[] this_Dist, float L_o, int T) {
            boolean[] TP = new boolean[this_Order.length];
            float[] S = _getS(this_Order,L_o,this_Dist);
            int last=0;
            ArrayList<int[]> Clusters = new ArrayList<int[]>();
            // Anything based on S is based on this_Order 
            // below, branch s corresponds to this_Order[s];
            ArrayList<Integer> Breakpoints = new ArrayList<Integer>();
            ArrayList<Integer> ForwardRuns = new ArrayList<Integer>();
            Breakpoints.add(0);
            //find all breakpoints
            for(int s=0;s<S.length-1;s++){
            	boolean tp = (S[s] * S[s+1] <= 0.0f); // true if sign of s and s+1 are non-equal
            	if(tp == true){
            		int length = s - last;
            		ForwardRuns.add(length); // refers to the forward run of the previous breakpoint
            		Breakpoints.add(s); // adds this breakpoint
            		last = s;
            	}else{
            		
            	}
            }
            ForwardRuns.add(S.length-last); // caps off the forward run array
            int tau = 10;
            // identify significant breakpoints with forward run length > Tau
            int b = 1; // always retain the first breakpoint, at zero...
            System.err.println("iterating through breakpoints...");
            while(b < Breakpoints.size()-1){
            	System.err.println("working on " + b);
            	if(ForwardRuns.get(b) < tau){
            		ForwardRuns.remove(b);
            		Breakpoints.remove(b);
            		int newForward = Breakpoints.get(b) - Breakpoints.get(b-1);
            		ForwardRuns.set(b-1, newForward);
            		// b is gone. there is a new b.
            		// no need to iterate.
            	}else{
            		// leave in place. is a significant run
            		b++;
            	}
            	// any remaining breakpoints likely have invalid run lengths
            }
            System.err.println("Done.\nAdding Clusters");
            //
            b=0;
            while(b<Breakpoints.size()-2){
            	int this_s = Breakpoints.get(b);
            	int this_fr= ForwardRuns.get(b);
            	int end = this_s+this_fr-1;
            	int[] cluster = new int[this_fr];
            	for(int i = this_s;i<=end;i++){
            		int index=i-this_s;
            		cluster[index]=i;
            	}
            	Clusters.add(cluster);
            	b++;
            }
            System.err.println("done.");
            /*
            for(int s=0;s<S.length-1;s++){
                TP[s] = (S[s] * S[s+1] <= 0.0f);
                int C = s-last+1;
                if(TP[s] == true){ 
                    //System.err.println("Found: C: " + C + " and last: " + last + " and s " + s);
                    if (C<10) continue;
                    int[] cluster = new int[C];
                    for(int i=last;i<=s;i++){
                        int index=i-last;
                        cluster[index]=i; // add branches to cluster.
                    }
                    Clusters.add(cluster);
                    last = s;
                }else if((s==S.length-1) && (last ==0)){ // if whole thing is one cluster
                	int[] cluster = new int[C];
                    for(int i=last;i<=s;i++){
                        int index=i-last;
                        cluster[index]=i; // add branches to cluster.
                    }
                    Clusters.add(cluster);
                }
            }
            */
            // Build ordered set of clusters // AND filter at the same time
          
            int c = 0;
            while(c < Clusters.size()-1){
                int[] this_Cluster = Clusters.get(c);
                int[] next_Cluster = Clusters.get(c+1);
                int this_n = this_Cluster.length;
                int next_n = next_Cluster.length;
                /*if(this_Cluster.length < T){
                    c++;
                    continue;
                }*/
                float main_Mean = _getMeanForCluster(this_Dist,this_Order,this_Cluster);
                float next_Mean = _getMeanForCluster(this_Dist,this_Order,next_Cluster);
                System.out.println("main: " + main_Mean + "("+this_n+") next: " + next_Mean +"(" + next_n + ")");
                //if(Clusters.get(c+1).length < T){
                if((this_n > T) && (next_n < T)){
                    if(next_Mean <= main_Mean){
                        int[] new_Cluster = ArrayUtils.addAll(this_Cluster,next_Cluster);
                        Clusters.set(c, new_Cluster);
                        Clusters.remove(c+1);
                        System.out.println("merging");
                    }else{
                        c++;
                    }
                }else if((this_n < T) && (next_n < T)){
                    int[] new_Cluster = ArrayUtils.addAll(this_Cluster,next_Cluster);
                    Clusters.set(c, new_Cluster);
                    Clusters.remove(c+1);
                    System.out.println("merging");
                }else if((this_n < T) && (next_n > T)){
                    if(next_Mean <= main_Mean){
                        int[] new_Cluster = ArrayUtils.addAll(this_Cluster,next_Cluster);
                        Clusters.set(c+1, new_Cluster);
                        Clusters.remove(c);
                        System.out.println("merging");
                    }else{
                        c++;
                    }
                }else{
                    c++;
                }
            }
            /*      
            Iterator<Integer[]> I = Clusters.iterator();
            int[][] Cs = new int[Clusters.size()][];
            int m = 0;
            while(I.hasNext()){
                Integer[] CLUSTER = I.next();
                int[] cluster = new int[CLUSTER.length];
                for(int i=0;i<CLUSTER.length;i++){
                    cluster[i]=(int) CLUSTER[i];
                }
                Cs[m]=cluster;
                m++;
            }
            return Cs;
            */
            return Clusters;
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