package networkCalcPackage;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

class Cluster {

    private static GCNMatrix DISS;
    private float INF = Float.POSITIVE_INFINITY;
    private Dendrogram Dendrogram;
    private int Criteria;

    public Cluster(GCNMatrix Similarities, int Crit) {
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
        DISS = new GCNMatrix(N, N);
        DISS = Operations.copyNames(Similarities.getRowNames(), DISS);
        for (int i = 0; i < N; i++) {
            float[] Row = Similarities.getRowByIndexAsDistance(i);
            DISS.addRow(Row);
        }

        Dendrogram = new Dendrogram(DISS);
        Dendrogram.getDendrogram(Criteria);
    }

    public ArrayList<int[]> dynamicTreeCut(int MinSize) {
            // implementation of dynamic tree cut and flashclust from horvath and langfelder et al
        // This will be pretty ham-handed.
        float[] Dist = Dendrogram.getHeights();
        int N = Dendrogram.getNumberOfBranches();
        int[] Order = Dendrogram.getDendroOrder();
        float cutoff = 0.99f * Dist[N - 2];
        System.err.println("Cutoff is 99% of dendrogram height: " + cutoff);

        ArrayList<int[]> Clusters = Dendrogram.staticCut(cutoff, MinSize);

        int c = 0;
        for (c = 0; c < Clusters.size(); c++) {
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
            if (Clusters.size() > N) {
                System.err.println("More clusters than genes. Something has gone horribly wrong.\n\n");
                System.exit(0);
            }
            //System.out.println("Processing cluster " + c + ", of size: " + this_N);
            if (this_N < MinSize) {
                c++;
                continue;
            }
            ArrayList<int[]> new_Clusters = _adaptiveTreeCut(this_cluster, MinSize);
            /*
             * The content of new_Clusters is a set of new branch IDs which:
             * contain all integers from i=0 to i<this_cluster.length
             * are origin-0
             * not necessarily sorted
             * no missing values
             * These branch IDs refer to the IDs of this_cluster.                
             */
            if (new_Clusters.size() == 0) {
                System.err.println("ZERO new clusters... whooops\n\n");
                System.exit(0);
            }
            if (new_Clusters.size() == 1) {
                //System.out.println("Found no new clusters");
                c++;
                continue;
            } else {
                int n = 0;
                Clusters.remove(c);
                while (n < new_Clusters.size()) {
                    int nc = c + n;
                    //  System.out.println("Inserting cluster at " + nc);
                    Clusters.add(nc, new_Clusters.get(n));
                    n++;
                }
                    //c++;
                //c=c+n+1;
            }
        }
        System.out.println("Done.");
            //System.out.println("Obtained " + Clusters.size() + " Clusters.");

        return Clusters;
    }

    private GCNMatrix _getDistForCluster(int[] cluster) {
        int thisN = cluster.length;
        GCNMatrix thisDist = new GCNMatrix(thisN, thisN);
        for (int j = 0; j < thisN; j++) {
            int Ind_i = cluster[j];
            float[] Row = new float[thisN];
            for (int J = 0; J < thisN; J++) {
                int Ind_j = cluster[J];
                float d = DISS.getValueByEntry(Ind_i, Ind_j);
                Row[J] = d;
            }
            thisDist.addRow(Row);
        }
        return thisDist;
    }

    private float _getMax(float[] array) {
        float max = 0.0f;
        for (int i = 0; i < array.length; i++) {
            if (array[i] > max) {
                max = array[i];
            }
        }
        return max;
    }

    private float _getMin(float[] array) {
        float min = INF;
        for (int i = 0; i < array.length; i++) {
            if (array[i] < min) {
                min = array[i];
            }
        }
        return min;
    }

    private ArrayList<int[]> _adaptiveTreeCut(int[] Cluster, int Tau) {
        GCNMatrix thisDist = _getDistForCluster(Cluster);
        Dendrogram this_Dendro = new Dendrogram(thisDist);
        this_Dendro.getDendrogram(Criteria);
        float[] this_Dist = this_Dendro.getHeights();
        int[] this_Order = this_Dendro.getDendroOrder();

        float L_naught = _getMean(this_Dist);
        float L_max = 0.5f * (L_naught + _getMax(this_Dist));
        float L_min = 0.5f * (L_naught + _getMin(this_Dist));
        ArrayList<int[]> clusters = new ArrayList<int[]>();
        if (Cluster.length < Tau) {
            clusters.add(Cluster);
            return clusters;
        }
        //System.out.println("working on cluster of size "+ Cluster.length);
        // max -> naught -> min continues to cut based on the distances until more than 1 cluster is obtained.
        // Basically, it biases for larger clusters.
        clusters = _treeCutCore(this_Order, this_Dist, L_max, Tau);
        //   System.out.println("Found " + clusters.size() + " clusters on this iteration (naught) " + L_naught);
        if (clusters.size() <= 1) {
            clusters = _treeCutCore(this_Order, this_Dist, L_naught, Tau);
            //     System.out.println("Found " + clusters.size() + " clusters on this iteration (lower) " + L_min );
        }
        if (clusters.size() <= 1) {
            clusters = _treeCutCore(this_Order, this_Dist, L_min, Tau);
            //     System.out.println("Found " + clusters.size() + " clusters on this iteration (upper) " + L_max);
        }

        if (clusters.size() == 0) {
            clusters.add(Cluster);
            return clusters;
        } else if (clusters.size() == 1) {
            clusters.remove(0);
            clusters.add(Cluster);
            return clusters;
        } else {
            for (int c = 0; c < clusters.size(); c++) {
                //	System.out.println("\t Sub Cluster size: " + clusters.get(c).length);
                int[] cluster = clusters.get(c);
                for (int i = 0; i < cluster.length; i++) {
                    cluster[i] = Cluster[cluster[i]];
                }
                clusters.set(c, cluster);
            }
        }
        return clusters;
    }

    private float[] _getS(int[] this_Order, float L, float[] this_Dist) {
        float[] S = new float[this_Order.length];
        for (int o = 0; o < this_Order.length; o++) {
            S[o] = this_Dist[this_Order[o]] - L;
            //     System.out.println(o + "\t" + this_Order[o] +"\t" + S[o] + "\t" + L);
        }
        return S;
    }

    private float _getMeanForCluster(float[] this_Dist, int[] this_Order, int[] this_Cluster) {
        float S = 0.0f;
        for (int c = 0; c < this_Cluster.length; c++) {
            S = S + this_Dist[this_Order[this_Cluster[c]]];
        }
        S = S / this_Cluster.length;
        return S;
    }

    private ArrayList<int[]> _treeCutCore(int[] this_Order, float[] this_Dist, float L_o, int T) {
        boolean[] TP = new boolean[this_Order.length];
        float[] S = _getS(this_Order, L_o, this_Dist);
        int last = 0;
        ArrayList<int[]> Clusters = new ArrayList<int[]>();
        ArrayList<Integer> Breakpoints = new ArrayList<Integer>();
        ArrayList<Integer> ForwardRuns = new ArrayList<Integer>();
        ArrayList<Integer> ReverseRuns = new ArrayList<Integer>();
        int tau = T;
        for (int s = 0; s < S.length - 1; s++) {
            boolean tp = (S[s] * S[s + 1] <= 0.0f);
            if (tp == true) {
                int RR = s - last;
                //if (RR <= tau) continue;
                int FR = S.length - s;
                for (int f = s + 1; f < S.length - 1; f++) {
                    if (S[f] * S[f + 1] <= 0.0f) {
                        FR = f - s;
                        break;
                    }
                }
                ForwardRuns.add(FR);
                Breakpoints.add(s);
                ReverseRuns.add(RR);
                //            System.err.println("adding " + s + " as a breakpoint, and " + FR + " as the forward run and " + RR + " as the reverse run");
                last = s;
            } else if (s == S.length - 2) {
                int RR = S.length - last - 1;
                //if (RR <= tau) continue;
                int FR = 0;
                int t = S.length - 1;
                ForwardRuns.add(FR);
                Breakpoints.add(t);
                ReverseRuns.add(RR);
                //      System.err.println("adding " + t + " as a breakpoint, and " + FR + " as the forward run and " + RR + " as the reverse run");
            } else {

            }
        }
        ForwardRuns.add(S.length - last);

        int b = 0;
        //System.err.println("iterating through breakpoints...");
        while (b < Breakpoints.size()) {
            	//System.err.println("working on " + b);
            //if((ForwardRuns.get(b) > tau) && (ReverseRuns.get(b) > tau)){
            if (ReverseRuns.get(b) > tau) {
                //  System.err.println("Significant breakpoint at " + Breakpoints.get(b) + " with FR: " + ForwardRuns.get(b) + " And RR: " + ReverseRuns.get(b));
                b++;
            } else {
                //System.err.println("Removing breakpoint at " + Breakpoints.get(b) + " with FR: " + ForwardRuns.get(b) + " And RR: " + ReverseRuns.get(b));
                Breakpoints.remove(b);
                ForwardRuns.remove(b);
                ReverseRuns.remove(b);
            }
        }
          //  System.err.println("Done.\nAdding Clusters");

        last = 0;
        b = 0;
        while (b < Breakpoints.size()) {
            int this_s = Breakpoints.get(b);
            int begin = last;
            int this_rr = this_s - last;
            int[] cluster = new int[this_rr + 1];
            for (int i = begin; i <= this_s; i++) {
                int index = i - begin;
                cluster[index] = i;
            }
            //    System.err.println("adding cluster of size " + cluster.length);
            Clusters.add(cluster);
            b++;
            last = this_s;
        }
            //System.err.println("done.");

        int c = 0;
        while (c < Clusters.size() - 1) {
            int[] this_Cluster = Clusters.get(c);
            int[] next_Cluster = Clusters.get(c + 1);
            int this_n = this_Cluster.length;
            int next_n = next_Cluster.length;
            float main_Mean = _getMeanForCluster(this_Dist, this_Order, this_Cluster);
            float next_Mean = _getMeanForCluster(this_Dist, this_Order, next_Cluster);

            if ((this_n > T) && (next_n < T)) {
                if (next_Mean <= main_Mean) {
                    int[] new_Cluster = ArrayUtils.addAll(this_Cluster, next_Cluster);
                    Clusters.set(c, new_Cluster);
                    Clusters.remove(c + 1);
                    //          System.out.println("merging");
                } else {
                    c++;
                }
            } else if ((this_n < T) && (next_n < T)) {
                int[] new_Cluster = ArrayUtils.addAll(this_Cluster, next_Cluster);
                Clusters.set(c, new_Cluster);
                Clusters.remove(c + 1);
                //    System.out.println("merging");
            } else if ((this_n < T) && (next_n > T)) {
                if (next_Mean <= main_Mean) {
                    int[] new_Cluster = ArrayUtils.addAll(this_Cluster, next_Cluster);
                    Clusters.set(c + 1, new_Cluster);
                    Clusters.remove(c);
                    //      System.out.println("merging");
                } else {
                    c++;
                }
            } else {
                c++;
            }

        }
        c = 0;
        while (c < Clusters.size() - 1) {
            int[] this_Cluster = Clusters.get(c);
            if (this_Cluster.length < T) {
                Clusters.remove(c);
                //    System.err.println("Removing cluster of size " + this_Cluster.length);
            } else {
                c++;
            }
        }
        return Clusters;
    }

    private float _getMean(float[] Distances) {
        float mean = 0.0f;
        for (int i = 0; i < Distances.length; i++) {
            mean = mean + Distances[i];
        }
        mean = mean / Distances.length;
        return mean;
    }

    private float[] _subsetArray(int[] index, float[] TA) {
        float[] subset = new float[index.length];
        for (int i = 0; i < index.length; i++) {
            subset[i] = TA[index[i]];
        }
        return subset;
    }

    private int[] _subsetArray(int[] index, int[] TA) {
        int[] subset = new int[index.length];
        for (int i = 0; i < index.length; i++) {
            subset[i] = TA[index[i]];
        }
        return subset;
    }

    private static void _clustersToFile(int[] Cluster, int M) {
        try {
            String nPath = "Cluster." + M + ".txt";
            PrintWriter writer = new PrintWriter(nPath, "UTF-8");
            for (int i = 0; i < Cluster.length; i++) {
                int node = Cluster[i];
                String name = DISS.getRowName(node);
                writer.println(name);
            }
            writer.close();
        } catch (Exception e) {
            // 
        }

    }

}
