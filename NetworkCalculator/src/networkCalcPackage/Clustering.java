package networkCalcPackage;

import java.io.PrintWriter;
/*
Adapted from Langfelder & Yau's adaptation of original Fortran code of Fionn Murtagh
http://cran.r-project.org/web/packages/flashClust/index.html
http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/flashClust/

Published here:
http://www.jstatsoft.org/v46/i11/paper

No changes, just a java implementation with perhaps some multithreading support ... we'll see.
*/

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
            float minD = INF;
            for(int i=0;i<N-1;i++){
                int J=-1;
                for(int j=i+1;j<N;j++){
                    if(DISS.getValueByEntry(i,j) > minD) continue;
                    minD = DISS.getValueByEntry(i,j);
                    J=j;
                }
                NN[i]=J;
                DISNN[i]=minD;
            }
            
            while(NCL>1){
                minD = INF;
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
        }
                
	/*
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
        */
        /*
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  DATA(N,M)         input data matrix,                      C
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C  Modified by Peter Langfelder, implemented bug fix
C  by Chi Ming Yau
C                                                            C
C------------------------------------------------------------C
      SUBROUTINE HC(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,
     X                FLAG,DISS)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DOUBLE PRECISION MEMBR(N),DISS(LEN)
      INTEGER IA(N),IB(N)
      DOUBLE PRECISION CRIT(N)
      DIMENSION NN(N),DISNN(N)
      LOGICAL FLAG(N)
      DOUBLE PRECISION INF
c     was 1D+20
      DATA INF/1.D+300/
c
c     unnecessary initialization of im jj jm to keep g77 -Wall happy
c
      IM = 0
      JJ = 0
      JM = 0

C
C  Initializations
C
      DO I=1,N
c         MEMBR(I)=1.
         FLAG(I)=.TRUE.
      ENDDO
      NCL=N
C
C  Construct dissimilarity matrix
C
C      DO I=1,N-1
C         DO J=I+1,N
C            IND=IOFFSET(N,I,J)
C            DISS(IND)=0.
C            DO K=1,M
C               DISS(IND)=DISS(IND)+(DATA(I,K)-DATA(J,K))**2
C            ENDDO
C            IF (IOPT.EQ.1) DISS(IND)=DISS(IND)/2.
C           (Above is done for the case of the min. var. method
C            where merging criteria are defined in terms of variances
C            rather than distances.)
C          ENDDO
C       ENDDO

C
C  Carry out an agglomeration - first create list of NNs
C


      DO I=1,N-1
         DMIN=INF
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            IF (DISS(IND).GE.DMIN) GOTO 500
               DMIN=DISS(IND)
               JM=J
  500    CONTINUE
         ENDDO
         NN(I)=JM
         DISNN(I)=DMIN
      ENDDO
C
  400 CONTINUE
C     Next, determine least diss. using list of NNs
      DMIN=INF
      DO I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 600
         IF (DISNN(I).GE.DMIN) GOTO 600
            DMIN=DISNN(I)
            IM=I
            JM=NN(I)
  600    CONTINUE
      ENDDO
      NCL=NCL-1
C
C  This allows an agglomeration to be carried out.
C
      I2=MIN0(IM,JM)
      J2=MAX0(IM,JM)
      IA(N-NCL)=I2
      IB(N-NCL)=J2
      CRIT(N-NCL)=DMIN
C
C  Update dissimilarities from new cluster.
C
      FLAG(J2)=.FALSE.
      DMIN=INF
      DO K=1,N
         IF (.NOT.FLAG(K)) GOTO 800
         IF (K.EQ.I2) GOTO 800
         X=MEMBR(I2)+MEMBR(J2)+MEMBR(K)
         IF (I2.LT.K) THEN
                           IND1=IOFFSET(N,I2,K)
                      ELSE
                           IND1=IOFFSET(N,K,I2)
         ENDIF
         IF (J2.LT.K) THEN
                           IND2=IOFFSET(N,J2,K)
                      ELSE
                           IND2=IOFFSET(N,K,J2)
         ENDIF
         IND3=IOFFSET(N,I2,J2)
         XX=DISS(IND3)
C
C  WARD'S MINIMUM VARIANCE METHOD - IOPT=1.
C
         IF (IOPT.EQ.1) THEN
            DISS(IND1)=(MEMBR(I2)+MEMBR(K))*DISS(IND1)+
     X                 (MEMBR(J2)+MEMBR(K))*DISS(IND2)-
     X                 MEMBR(K)*XX
            DISS(IND1)=DISS(IND1)/X
         ENDIF
C
C  SINGLE LINK METHOD - IOPT=2.
C
         IF (IOPT.EQ.2) THEN
            DISS(IND1)=MIN(DISS(IND1),DISS(IND2))
         ENDIF
C
C  COMPLETE LINK METHOD - IOPT=3.
C
         IF (IOPT.EQ.3) THEN
            DISS(IND1)=MAX(DISS(IND1),DISS(IND2))
         ENDIF
C
C  AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4.
C
         IF (IOPT.EQ.4) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2))/
     X                 (MEMBR(I2)+MEMBR(J2))
         ENDIF
C
C  MCQUITTY'S METHOD - IOPT=5.
C
         IF (IOPT.EQ.5) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)
         ENDIF
C
C  MEDIAN (GOWER'S) METHOD - IOPT=6.
C
         IF (IOPT.EQ.6) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)-0.25*XX
         ENDIF
C
C  CENTROID METHOD - IOPT=7.
C
         IF (IOPT.EQ.7) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2)-
     X          MEMBR(I2)*MEMBR(J2)*XX/(MEMBR(I2)+MEMBR(J2)))/
     X          (MEMBR(I2)+MEMBR(J2))
            ENDIF
C
         IF (I2.GT.K) GOTO 800
         IF (DISS(IND1).GE.DMIN) GOTO 800
            DMIN=DISS(IND1)
            JJ=K
  800    CONTINUE
      ENDDO
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
      DISNN(I2)=DMIN
      NN(I2)=JJ
C
C  Update list of NNs insofar as this is required.
C  This part modified by Chi Ming Yau and PL. For methods IOPT=6 and 7
C  use modified updating of nearest neighbors that is a bit slower but
C  necessary.

      IF (IOPT.GT.5) THEN
        DO I=1,N-1
          IF (.NOT.FLAG(I)) GOTO 900
          IF (I.EQ.I2) GOTO 850
          IF (NN(I).EQ.I2) GOTO 850
          IF (NN(I).EQ.J2) GOTO 850
C  Compare DISNN(I) with updated DISS between I and I2
          IF (I2.LT.I) THEN
                            IND=IOFFSET(N,I2,I)
                       ELSE
                            IND=IOFFSET(N,I,I2)
          ENDIF
          DMIN=DISS(IND)
          IF (DMIN.GE.DISNN(I)) GOTO 900
             DISNN(I)=DMIN
             NN(I)=I2
          GOTO 900
 850      CONTINUE
C        (Redetermine NN of I:)
          DMIN=INF
          DO J=I+1,N
             IND=IOFFSET(N,I,J)
             IF (.NOT.FLAG(J)) GOTO 870
             IF (I.EQ.J) GOTO 870
             IF (DISS(IND).GE.DMIN) GOTO 870
                DMIN=DISS(IND)
                JJ=J
 870         CONTINUE
          ENDDO
          NN(I)=JJ
          DISNN(I)=DMIN
 900      CONTINUE
        ENDDO
      ELSE
C       For methods IOPT<6 use the original fast update.
        DO I=1,N-1
           IF (.NOT.FLAG(I)) GOTO 901
           IF (NN(I).EQ.I2) GOTO 851
           IF (NN(I).EQ.J2) GOTO 851
           GOTO 901
  851      CONTINUE
C          (Redetermine NN of I:)
           DMIN=INF
           DO J=I+1,N
              IND=IOFFSET(N,I,J)
              IF (.NOT.FLAG(J)) GOTO 871
              IF (I.EQ.J) GOTO 871
              IF (DISS(IND).GE.DMIN) GOTO 871
                 DMIN=DISS(IND)
                 JJ=J
  871         CONTINUE
           ENDDO
           NN(I)=JJ
           DISNN(I)=DMIN
  901      CONTINUE
        ENDDO
      ENDIF

C
C  Repeat previous steps until N-1 agglomerations carried out.
C
      IF (NCL.GT.1) GOTO 400
C
C
      RETURN
      END
C
C
      FUNCTION IOFFSET(N,I,J)
C  Map row I and column J of upper half diagonal symmetric matrix
C  onto vector.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C  Convert integer I to a double
C  This hopefully prevents overflow errors when I^2 is greater than
C  2^31.
      IF (N.GT.32768) THEN
         XI = DBLE(I)
         IOFFSET=J+NINT( (XI-1)*N - (XI*(XI+1))/2)
      ELSE
         IOFFSET=J+(I-1)*N-(I*(I+1))/2
      ENDIF
      RETURN
      END
        */

}