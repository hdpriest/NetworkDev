package networkCalcPackage;

import java.io.PrintWriter;


class Clustering {
        private GCNMatrix DISS;
        //private float[][] DISS;
        private int Critereon;
        private int N;
        private int[] IA;
        private int[] IB;
        private int NCL;
        private float[] CRIT;
        private int[] MEMBR;
        private int[] NN;
        private float[] DISNN;
        private boolean[] FLAG;
        private float INF = Float.POSITIVE_INFINITY;
        
        
	
	
	public Clustering (GCNMatrix Similarities) {
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
            N = Similarities.getNumRows();
            IA = new int[N];
            IB = new int[N];
            CRIT = new float[N];
            MEMBR = new int[N];
            NN = new int[N];
            DISNN = new float[N];
            FLAG = new boolean[N];
            DISS = new GCNMatrix(N,N);
            NCL = N;
            for (int i = 0; i < N ; i++) {
                FLAG[i] = true; // all objects start agglomerable
                MEMBR[i] = 1; // all clusters start with 1 entry
                float[] Row = Similarities.getRowByIndexAsDistance(i);
                DISS.addRow(Row);
	    }
	}
        
        public void Cluster (int Crit){
            Critereon = Crit;
                        
            // Find NN for all I
            
            for(int i=0;i<N-1;i++){
                int J=-1;
                float minD = INF;
                for(int j=i+1;j<N;j++){
                    if(DISS.getValueByEntry(i,j) > minD) continue;
                    minD = DISS.getValueByEntry(i,j);
                    J=j;
                    
                }
                NN[i]=J;
                DISNN[i]=minD;
            }
            
            while(NCL>1){
                float minD = INF;
                int IM=0;
                int JM=0;
                for(int i=0;i<N-1;i++){
                    if(FLAG[i] != true) continue;
                    if(DISNN[i] > minD) continue;
                    minD = DISNN[i];
                    IM=i;
                    JM=NN[i];
                }
                NCL--;
                int I2 = (IM < JM) ? IM : JM;
                int J2 = (IM > JM) ? IM : JM; // I2 < J2
                // So, at each step, IA holds the root node, IB holds the node that was merged into it, and CRIT holds the dissimilarity.
                // These are in order, meaning that the first entry is that pair of NNs which was closest. The 2nd entry is the 2nd closest.
                // So, we can obtain the order of merges which will tell us the horizontal order of the dendrogram
                // There's no gty that IA[0] and IA[1] are the same root...
                
                IA[N-NCL]=I2;
                IB[N-NCL]=J2;
                CRIT[N-NCL]=minD;

                FLAG[J2] = false;

                minD = INF;
                int JJ=-1;
                // UPDATE DISTANCES
                for(int k=0;k<N;k++){
                    if(FLAG[k] != true) continue;
                    if(k==I2) continue;
                    int X = MEMBR[I2] + MEMBR[J2] + MEMBR[k];
                    float XX = DISS.getValueByEntry(I2,J2);
                    // IND1 in fortran = K,I2
                    // IND2 in fortran = K,J2 
                    float IK = DISS.getValueByEntry(I2,k); // IND 1
                    float JK = DISS.getValueByEntry(J2,k); // IND 2
                    float newDiss=1.0f;
                    switch (Critereon) { // works on String or int (and others) TODO encode cases...
                        //case "ward":
                        case 1:
                            newDiss = ((MEMBR[I2]+MEMBR[k]) * IK) + ((MEMBR[J2]+MEMBR[k]) * JK) - (MEMBR[k]*XX);
                            newDiss = newDiss/X;
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
                            newDiss = (MEMBR[I2]*IK+MEMBR[J2]*JK)/(MEMBR[I2]+MEMBR[J2]);
                            break;
                        //case "mcquitty":
                        case 5:
                            newDiss = 0.5f*IK+0.5f*JK; 
                            break;
                        //case "median":
                        case 6:
                            // also Gower's
                            newDiss = 0.5f*IK + 0.5f*JK - 0.25f*XX;
                            break;     
                        //case "centroid":
                        case 7:
                            newDiss = (MEMBR[I2]*IK+MEMBR[J2]*JK-MEMBR[I2]*MEMBR[J2]*XX/(MEMBR[I2]+MEMBR[J2]))/(MEMBR[I2]+MEMBR[J2]);
                            break;                        
                        //etc
                    }
                    DISS.setValueByEntry(newDiss, I2, k);
                    if(I2 > k) continue; // what is this for??
                    if(newDiss > minD) continue;
                    minD = newDiss;
                    JJ=k; // only designate k as the NN if k > I2 .... can only merge left??
                }
                MEMBR[I2] = MEMBR[I2] + MEMBR[J2];
                DISNN[I2] = minD;
                NN[I2]=JJ;

                /*
                 OK. Update NNs:
                if NN(i) == I2, update DISNN(i)
                if NN(i) == J2, update DISNN, update NN(i)=I2
                 Yau & PL have a better updater for 6 & 7
                */ 
                if(Critereon > 5){
                    for(int i=0;i<N-1;i++){
                        if(FLAG[i] != true) continue; // FLAG[i] == false >> i has been clustered
                        if((i == I2) || (NN[i] == I2) || (NN[i] == J2)){ 
                            // If i != I2, i's NN isn't I2 or J2, it's NN didn't change, unless the new I2 Distance is < DISNN
                            minD = INF;
                            for(int j=i+1;j<N;j++){
                                if(FLAG[j] != true) continue;
                                if (i == j) continue; // Why is this here??
                                float D = DISS.getValueByEntry(i,j);
                                if(D >= minD) continue;
                                minD = D;
                                JJ = j;
                            }
                            NN[i]=JJ;
                            DISNN[i]=minD;
                        }else{
                            float D = DISS.getValueByEntry(i,I2);
                            if(D >= DISNN[i]) continue; // Skip this if the updated i-I2 distance is GE the current NN of i
                            DISNN[i] = D;
                            NN[i] = I2;
                        }
                    }
                }else{
                    for(int i=0;i<N-1;i++){
                        if(FLAG[i] != true) continue; // omits J2
                        if((NN[i] == I2) || (NN[i] == J2)){ // Otherwise, NN didn't change
                            minD = INF;
                            for(int j=i+1;j<N;j++){
                                if(FLAG[j] != true) continue;
                                if(i == j) continue; // again. seemed useless.
                                float D = DISS.getValueByEntry(i,j);
                                if(D > minD) continue;
                                minD = D;
                                JJ = j;
                            }
                            NN[i] = JJ;
                            DISNN[i] = minD;
                        }
                    }
                }
            }
            System.err.println("Done clustering...\n");
        }
        
        public void getClusters () {
            /*
        private GCNMatrix DISS;
        //private float[][] DISS;
        private int Critereon;
        private int N;
        private int[] IA;
        private int[] IB;
        private int NCL;
        private float[] CRIT;
        private int[] MEMBR;
        private int[] NN;
        private float[] DISNN;
        private boolean[] FLAG;
        private float INF = Float.POSITIVE_INFINITY;
            */
            float cutoff = 0.99f * (CRIT[N-1]);
            System.out.println("Cutoff is 99% of dendrogram height: " +cutoff);
            for(int i=1;i<N;i++){
                System.out.println(IA[i] + " " + IB[i] + " " + CRIT[i]);
            }
            System.exit(0);
            /*
            I would say the desired return from this is... a list of clusters, and the clusters are lists of nodes.
            So, it could really be ArrayList<array> or something similar...
            But, keep in mind the R package ALSO has dynamicTreeCut - which calculates the diff threshold
            Also, it has a min module size, which would keep you from reporting a bunch of singletons as
            single gene modules... 
            */
            
        }
        
        public void getDendroOrder () {
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
            int[] IIA = new int[N];
            int[] IIB = new int[N];
            int[] IORDER = new int [N];
            for(int i=0;i<N;i++){
                IIA[i]=IA[i];
                IIB[i]=IB[i];
            }
            
            for(int i=0;i<N-2;i++){
                int k = (IA[i] < IB[i]) ? IA[i] : IB[i];
                for(int j=i+1;j<N-1;j++){
                    if (IA[j] == k) IIA[j]=(-1 * i);
                    if (IB[j] == k) IIB[j]=(-1 * i);
                }
            }
            /// I don't understand the point of the above and below stanzas yet...
            for(int i=0;i<N-1;i++){
                IIA[i] = -1 * IIA[i];
                IIB[i] = -1 * IIB[i];
            }
            
            for(int i=0;i<N-1;i++){
                if((IIA[i] > 0) && (IIB[i] < 0)){
                   int K = IIA[i];
                   IIA[i] = IIB[i];
                   IIB[i] = K;
                }
                if((IIA[i] > 0) && (IIB[i] > 0)){
                    int K1 = (IIA[i] < IIB[i]) ? IIA[i] : IIB[i];
                    int K2 = (IIA[i] > IIB[i]) ? IIA[i] : IIB[i];
                    IIA[i] = K1;
                    IIB[i] = K2;
                }
            }
            IORDER[1] = IIA[N-1];
            IORDER[2] = IIB[N-1];
            int LOC = 2;
            for(int i=N-3;i>=0;i--){
                for(int j=0;j<LOC;j++){
                    if(IORDER[j] == i){
                        IORDER[j] = IIA[i];
                        if(j == LOC){
                            LOC=LOC+1;
                            IORDER[LOC]=IIB[i];
                        }else{
                            LOC=LOC+1;
                            for(int k=LOC;k>=j+2;k--){
                                IORDER[k] = IORDER[k-1];
                            }
                            IORDER[j+1]=IIB[i];
                        }
                        break;
                    }
                }
            }
            for(int i=0;i<N;i++){
                IORDER[i]=-1 * IORDER[i];
                System.out.println(i+ " "+ IORDER[i]);
            }
            
            
        }
                
	/*
	
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
        */
        
}