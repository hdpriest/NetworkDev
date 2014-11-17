package networkCalcPackage;

import java.util.Iterator;
import java.util.List;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.apache.commons.lang3.ArrayUtils;

class Dendrogram {

    private GCNMatrix DISS;
    private int Critereon;
    private int N;
    private int[] IA;
    private int[] IB;
    private int NCL;
    private float[] CRIT;
    private int[] ORDER;
    private int[] MEMBR;
    private int[] NN;
    private float[] DISNN;
    private boolean[] FLAG;
    private float INF = Float.POSITIVE_INFINITY;

    public Dendrogram(GCNMatrix Distances) {
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

        DISS = Distances;
        N = DISS.getNumRows();
        NCL = N;
        FLAG = new boolean[N];
        NN = new int[N];
        DISNN = new float[N];
        IA = new int[N];
        IB = new int[N];
        CRIT = new float[N];
        MEMBR = new int[N];
        for (int i = 0; i < N; i++) {
            FLAG[i] = true;
            MEMBR[i] = 1;
        }
    }

    public void getDendrogram(int Crit) {
        Critereon = Crit;

        for (int i = 0; i < N - 1; i++) {
            int J = -1;
            float minD = INF;
            for (int j = i + 1; j < N; j++) {
                if (DISS.getValueByEntry(i, j) > minD) {
                    continue;
                }
                minD = DISS.getValueByEntry(i, j);
                J = j;

            }
            NN[i] = J;
            DISNN[i] = minD;
        }

        while (NCL > 1) {
            float minD = INF;
            int IM = 0;
            int JM = 0;
            for (int i = 0; i < N - 1; i++) {
                if (FLAG[i] != true) {
                    continue;
                }
                if (DISNN[i] > minD) {
                    continue;
                }
                minD = DISNN[i];
                IM = i;
                JM = NN[i];
            }
            NCL--;
            int I2 = (IM < JM) ? IM : JM;
            int J2 = (IM > JM) ? IM : JM;

            IA[N - NCL - 1] = I2;
            IB[N - NCL - 1] = J2;
            CRIT[N - NCL - 1] = minD;

            FLAG[J2] = false;

            minD = INF;
            int JJ = -1;
            for (int k = 0; k < N; k++) {
                if (FLAG[k] != true) {
                    continue;
                }
                if (k == I2) {
                    continue;
                }
                int X = MEMBR[I2] + MEMBR[J2] + MEMBR[k];
                float XX = DISS.getValueByEntry(I2, J2);
                float IK = DISS.getValueByEntry(I2, k);
                float JK = DISS.getValueByEntry(J2, k);
                float newDiss = 1.0f;
                switch (Critereon) { // works on String or int (and others) TODO encode cases...
                    //case "ward":
                    case 1:
                        newDiss = ((MEMBR[I2] + MEMBR[k]) * IK) + ((MEMBR[J2] + MEMBR[k]) * JK) - (MEMBR[k] * XX);
                        newDiss = newDiss / X;
                        break;
                    //case "single":
                    case 2:
                        newDiss = (IK <= JK) ? IK : JK;
                        break;
                    //case "complete":
                    case 3:
                        newDiss = (IK >= JK) ? IK : JK;
                        break;
                    //case "average":
                    case 4:
                        newDiss = (MEMBR[I2] * IK + MEMBR[J2] * JK) / (MEMBR[I2] + MEMBR[J2]);
                        break;
                    //case "mcquitty":
                    case 5:
                        newDiss = 0.5f * IK + 0.5f * JK;
                        break;
                    //case "median":
                    case 6:
                        // also Gower's
                        newDiss = 0.5f * IK + 0.5f * JK - 0.25f * XX;
                        break;
                    //case "centroid":
                    case 7:
                        newDiss = (MEMBR[I2] * IK + MEMBR[J2] * JK - MEMBR[I2] * MEMBR[J2] * XX / (MEMBR[I2] + MEMBR[J2])) / (MEMBR[I2] + MEMBR[J2]);
                        break;
                    //etc
                    }
                DISS.setValueByEntry(newDiss, I2, k);
                if (I2 > k) {
                    continue;
                }
                if (newDiss > minD) {
                    continue;
                }
                minD = newDiss;
                JJ = k;
            }
            MEMBR[I2] = MEMBR[I2] + MEMBR[J2];
            DISNN[I2] = minD;
            NN[I2] = JJ;

            if (Critereon > 5) {
                for (int i = 0; i < N - 1; i++) {
                    if (FLAG[i] != true) {
                        continue;
                    }
                    if ((i == I2) || (NN[i] == I2) || (NN[i] == J2)) {
                        minD = INF;
                        for (int j = i + 1; j < N; j++) {
                            if (FLAG[j] != true) {
                                continue;
                            }
                            if (i == j) {
                                continue;
                            }
                            float D = DISS.getValueByEntry(i, j);
                            if (D >= minD) {
                                continue;
                            }
                            minD = D;
                            JJ = j;
                        }
                        NN[i] = JJ;
                        DISNN[i] = minD;
                    } else {
                        float D = DISS.getValueByEntry(i, I2);
                        if (D >= DISNN[i]) {
                            continue;
                        }
                        DISNN[i] = D;
                        NN[i] = I2;
                    }
                }
            } else {
                for (int i = 0; i < N - 1; i++) {
                    if (FLAG[i] != true) {
                        continue;
                    }
                    if ((NN[i] == I2) || (NN[i] == J2)) {
                        minD = INF;
                        for (int j = i + 1; j < N; j++) {
                            if (FLAG[j] != true) {
                                continue;
                            }
                            if (i == j) {
                                continue;
                            }
                            float D = DISS.getValueByEntry(i, j);
                            if (D > minD) {
                                continue;
                            }
                            minD = D;
                            JJ = j;
                        }
                        NN[i] = JJ;
                        DISNN[i] = minD;
                    }
                }
            }
        }
    }

    public ArrayList<int[]> staticCut(float cutoff, int MinSize) {
        int[][] I_Clusters = new int[N][];
        for (int i = 0; i < N; i++) {
            int[] clust = new int[1];
            int branch_id = i;
            clust[0] = branch_id;
            I_Clusters[i] = clust;
        }
        for (int i = 0; i < N - 1; i++) {
            if (CRIT[i] > cutoff) {
                continue;
            }
            if (I_Clusters[IB[i]] == null) {
                continue;
            }
            int[] i_clust = ArrayUtils.addAll(I_Clusters[IA[i]], I_Clusters[IB[i]]);
            I_Clusters[IA[i]] = i_clust;
            I_Clusters[IB[i]] = null;

        }
        ArrayList<int[]> Clusters = new ArrayList<int[]>();
        for (int i = 0; i < N - 1; i++) {
            if (I_Clusters[i] == null) {
                continue;
            }
            if (I_Clusters[i].length == 1) {
                continue;
            }
            Clusters.add(I_Clusters[i]);
        }
        return Clusters;
    }

    public int[] getMergeLeaves() { // pretty sure thats not a leaf.
        return IB;
    }

    public int[] getMergeRoots() {
        return IA;
    }

    public float[] getHeights() {
        return CRIT;
    }

    public int getNumberOfBranches() {
        return N;
    }

    private float _getMean(float[] Distances) {
        float mean = 0.0f;
        for (int i = 0; i < Distances.length; i++) {
            mean = mean + Distances[i];
        }
        mean = mean / Distances.length;
        return mean;
    }

    public int[] getDendroOrder() {
        /*            
         C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
         C                                                               C
         C  Given a HIERARCHIC CLUSTERING, described as a sequence of    C
         C  agglomerations, prepare the seq. of aggloms. and "horiz."    C
         C  order of objects for plotting the dendrogram using S routine C
         C  'plclust'.                                                   C
         C                                                               C
         C  Parameters:                                                  C
         C                                                               C
         C  IA, IB:       vectors of dimension N defining the agglomer-  C
         C                 ations.                                       C
         C  IIA, IIB:     used to store IA and IB values differently     C
         C                (in form needed for S command 'plclust'        C
         C  IORDER:       "horiz." order of objects for dendrogram       C
         C                                                               C
         C  F. Murtagh, ESA/ESO/STECF, Garching, June 1991               C
         C                                                               C
         C  HISTORY                                                      C
         C                                                               C
         C  Adapted from routine HCASS, which additionally determines    C
         C   cluster assignments at all levels, at extra comput. expense C
         C                                                               C
         C  This routine copied by Peter Langfelder from the source      C
         C  of R package stats.                                          C
         C                                                               C
         C---------------------------------------------------------------C
         SUBROUTINE HCASS2(N,IA,IB,IORDER,IIA,IIB)
         c Args
         INTEGER N,IA(N),IB(N),IORDER(N),IIA(N),IIB(N)
         c Var
         INTEGER I, J, K, K1, K2, LOC
         C
         C     Following bit is to get seq. of merges into format acceptable to plclust
         C     I coded clusters as lowest seq. no. of constituents; S's 'hclust' codes
         C     singletons as -ve numbers, and non-singletons with their seq. nos.
         C
         */
        int N = IA.length;
        int[] IIA = new int[N];
        int[] IIB = new int[N];
        int[] IORDER = new int[N];
        for (int i = 0; i < N; i++) {
            IIA[i] = IA[i];
            IIB[i] = IB[i];
        }

        for (int i = 0; i <= N - 1; i++) {
            int k = (IA[i] < IB[i]) ? IA[i] : IB[i];
            for (int j = i + 1; j <= N - 1; j++) {
                if (IA[j] == k) {
                    IIA[j] = (-1 * i);
                }
                if (IB[j] == k) {
                    IIB[j] = (-1 * i);
                }
            }
        }

        for (int i = 0; i <= N - 1; i++) {
            IIA[i] = -1 * IIA[i];
            IIB[i] = -1 * IIB[i];
        }

        for (int i = 0; i <= N - 1; i++) {
            if ((IIA[i] > 0) && (IIB[i] < 0)) {
                int K = IIA[i];
                IIA[i] = IIB[i];
                IIB[i] = K;
            }
            if ((IIA[i] > 0) && (IIB[i] > 0)) {
                int K1 = (IIA[i] < IIB[i]) ? IIA[i] : IIB[i];
                int K2 = (IIA[i] > IIB[i]) ? IIA[i] : IIB[i];
                IIA[i] = K1;
                IIB[i] = K2;
            }
        }
        IORDER[0] = IIA[N - 2];
        IORDER[1] = IIB[N - 2];

        int LOC = 1;
        for (int i = N - 3; i >= 0; i--) {
            for (int j = 0; j <= LOC; j++) {
                if (IORDER[j] == i) {
                    IORDER[j] = IIA[i];
                    if (j == LOC) {
                        LOC = LOC + 1;
                        IORDER[LOC] = IIB[i];
                    } else {
                        LOC = LOC + 1;
                        for (int k = LOC; k >= j + 2; k--) {
                            IORDER[k] = IORDER[k - 1];
                        }
                        IORDER[j + 1] = IIB[i];
                    }
                    break;
                }
            }
        }
        for (int i = 0; i < N; i++) {
            IORDER[i] = -1 * IORDER[i];
        }
        ORDER = IORDER;
        return IORDER;
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
                writer.println(node);
            }
            writer.close();
        } catch (Exception e) {
        }

    }

}
