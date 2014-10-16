package networkCalcPackage;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;





public class Operations {
	
	public static float GINI (float[] array1,float[] array2){
		float GINI_coeff;
		float Numerator=0.0f;
		float Denominator=0.0f;
		int[] sortRanks1 = getIndicesInOrder(array1);
		int[] sortRanks2 = getIndicesInOrder(array2);
		for(int i=0;i<array1.length;i++){
			float v2=((2*(i+1))-array1.length-1) * array1[sortRanks1[i]];
			float v1=((2*(i+1))-array1.length-1) * array1[sortRanks2[i]];
			Denominator+=v2;
			Numerator += v1;
		}
		GINI_coeff=Numerator/Denominator;
		return GINI_coeff;
	}
	
        public static GCNMatrix copyNames (String[] names, GCNMatrix NetB){
            /// if you always copy Rows->rows and columns, you can't go wrong
            NetB.setRowNames(names);
            NetB.setColumnNames(names);
            return NetB;
        }
        
	public static int[] getIndicesInOrder(float[] array) {
	    Map<Integer, Float> map = new HashMap<Integer, Float>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);

	    List<Entry<Integer, Float>> l = 
	                           new ArrayList<Entry<Integer, Float>>(map.entrySet());

	    Collections.sort(l, new Comparator<Entry<?, Float>>() {
	            @Override
	            public int compare(Entry<?, Float> e1, Entry<?, Float> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });

	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();

	    return result;
	}
	/*
	public static GCNMatrix calculateGINIcoefficient (GCNMatrix InputFrame){
		int D = InputFrame.getNumRows();
		GCNMatrix Similarity = new GCNMatrix(D,D);
                Similarity = Operations.copyNames(InputFrame, Similarity);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				float GCC1=(float) 0.0;
				float GCC2=(float) 0.0;
				if(i==j){
					GCC1 = (float) 1.0;
					GCC2 = (float) 0.0;
				}else{
					float[] I_data = InputFrame.getRowByIndex(i);
					float[] J_data = InputFrame.getRowByIndex(j);
					GCC1 = GINI(I_data,J_data);
					GCC2 = GINI(J_data,I_data);
					//System.err.println();
				}
				if(Math.abs(GCC1)<Math.abs(GCC2)){
					Similarity.setValueByEntry(GCC1,i,j);
					Similarity.setValueByEntry(GCC1,j,i);
				}else if (Math.abs(GCC2)<Math.abs(GCC1)){
					Similarity.setValueByEntry(GCC2,i,j);
					Similarity.setValueByEntry(GCC2,j,i);
				}else{
					Similarity.setValueByEntry(GCC1,i,j);
					Similarity.setValueByEntry(GCC1,j,i);
				}
			}
		}
		return Similarity;
	}
	
	public static GCNMatrix calculateGINIcoefficient (GCNMatrix Expression,int Threads){ 
		//// it really seems as though all these concurrent methods
		//// could be condensed - the matrix handling is all fairly uniform, its the calls that differ.
		int D = Expression.getNumRows();
		GCNMatrix Similarity = new GCNMatrix(D,D);
                Similarity = Operations.copyNames(Expression, Similarity);
		ExecutorService pool = Executors.newFixedThreadPool(Threads);
		ExecutorCompletionService<HashMap<String,Float>> completionService = new ExecutorCompletionService<>(pool);
		List<Future<HashMap<String,Float>>> taskList = new ArrayList<Future<HashMap<String,Float>>>();
		ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				String S = i+"-"+j;
				queue.add(S);
			}
		}
		for ( int i = 0; i < Threads; i++ ) {
			Callable<HashMap<String,Float>> worker = new ConcurrentProcessing(Expression,queue,"gini");
			Future<HashMap<String,Float>> submit = completionService.submit(worker);
			taskList.add(submit);  
		}
		
		
		
		for(int t=0;t<Threads;t++){
			try{
				HashMap<String,Float> hm = completionService.take().get();
				for(Map.Entry<String,Float> entry : hm.entrySet()){
					String s = entry.getKey();
					String[] S = s.split("-");
					Float d = entry.getValue();
					int i = Integer.parseInt(S[0]);
					int j = Integer.parseInt(S[1]);
					Similarity.setValueByEntry((float) d,i,j);
					Similarity.setValueByEntry((float) d,j,i);
				}
			}catch(InterruptedException e){
				e.printStackTrace();
			}catch (ExecutionException e){
				e.printStackTrace();
			}
			
		}
		pool.shutdown();
		return Similarity;
	}
	*/
	public static GCNMatrix compareNetworksViaTOM (GCNMatrix Net1, GCNMatrix Net2){
		int D = Net1.getNumRows();
		GCNMatrix ReturnFrame = new GCNMatrix(D,D);
                ReturnFrame = Operations.copyNames(Net1.getRowNames(), ReturnFrame);
		for(int i=0;i<D;i++){
			
			for(int j=0;j<D;j++){
				float T=0;
				if(i==j){
					float product=0;
					float i_k = Net1.findK(i, i);
					float j_k = Net2.findK(j,j);
					for(int u=0;u<D;u++){
						if((u != i) && (u != j) && (Net1.testValue(i, u)) && (Net2.testValue(j, u))){
							float i_v = Net1.getValueByEntry(i,u);
							float j_v = Net2.getValueByEntry(j,u);
							float max = Math.max(i_v,j_v);
							product += i_v * j_v / max;
							/// if node is not connected to anything, all products are zero
						}
					}
					float k_min = Math.min(i_k, j_k);
					float DFIJ=0;
					T=(product+DFIJ)/(k_min + 1 - DFIJ); // if one node unconnected, = 0+0/0+1-0
					// if IJ are totally connected, all products > 0, but < kmin
					//T=(product+DFIJ)/(k_min + 1);
				}else{
				}
				ReturnFrame.setValueByEntry(T, i, j);
			}
		}
		return ReturnFrame;
	}

	public static GCNMatrix calculateTOM (GCNMatrix Adjacency, int Threads){
		int D = Adjacency.getNumRows();
		GCNMatrix ReturnMatrix = new GCNMatrix(D,D);
                ReturnMatrix = Operations.copyNames(Adjacency.getRowNames(), ReturnMatrix);
		ExecutorService pool = Executors.newFixedThreadPool(Threads);
		ExecutorCompletionService<HashMap<String,float[]>> completionService = new ExecutorCompletionService<>(pool);
		List<Future<HashMap<String,float[]>>> taskList = new ArrayList<Future<HashMap<String,float[]>>>();
		ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		//System.err.println("Processing topological overlap using " + Threads + " threads.");
		for(int i=0;i<D;i++){
			String S = String.valueOf(i);
			queue.add(S);
		}
		for ( int i = 0; i < Threads; i++ ) {
			Callable<HashMap<String,float[]>> worker = new ConcurrentProcessing(Adjacency,queue,"tom");
			Future<HashMap<String,float[]>> submit = completionService.submit(worker);
			taskList.add(submit);  
                        
		}
		
		
		
		for(int t=0;t<Threads;t++){
			try{
				HashMap<String,float[]> hm = completionService.take().get();
				//System.err.println("obtained result for thread " + t);
				int r = 0;
				for(Map.Entry<String,float[]> entry : hm.entrySet()){
                                        String s = entry.getKey();
                                        int i = Integer.valueOf(s);
                                        int size = D-i;
                                        float[] d = new float[size];
                                        d = entry.getValue();
                                        for(int j=0;j<d.length;j++){
                                            int coord = j+i;
                                            ReturnMatrix.setValueByEntry(d[j],i,coord);
                                            r++;
                                        }
				}
				//System.err.println("Processed "+ r + " records");
			}catch(InterruptedException e){
				e.printStackTrace();
			}catch (ExecutionException e){
				e.printStackTrace();
			}
			//System.err.println("Thread " + t + " complete.");
		}
		//System.err.println("Done.");
		pool.shutdownNow();
		return ReturnMatrix;
	}
	
	private static void explode (){
		System.err.println("No correlation got to this point!?");
		System.exit(0);
	}
	
	public static GCNMatrix calculateAdjacency (ExpressionFrame Expression,String corr, String adj, float mu, float alpha, int Threads){
		
		// prepwork
		int D = Expression.getNumRows();
		GCNMatrix Adjacency = new GCNMatrix(D,D);
                 Adjacency = Operations.copyNames(Expression.getRowNames(), Adjacency);
                switch (corr) {
                    case "gini": Expression.calculateGiniSums(); // half of every gini coeff is a pre-calculable value
                    break;
                    case "pcc": Expression.calculateMeans(); // most of a pcc is pre-calculateable;
                    break;
                    default: explode();
                    break;
                }
        
        // thread prep
		ExecutorService pool = Executors.newFixedThreadPool(Threads);
		ExecutorCompletionService<HashMap<String,float[]>> completionService = new ExecutorCompletionService<>(pool);
		List<Future<HashMap<String,float[]>>> taskList = new ArrayList<Future<HashMap<String,float[]>>>();
		ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		//System.err.println("Processing adjacency using " + Threads + " threads.");
		
		//queue prep
		for(int i=0;i<D;i++){
			String S = String.valueOf(i);
			queue.add(S);
		}
		
		//add the tasks
		for ( int i = 0; i < Threads; i++ ) {
			Callable<HashMap<String,float[]>> worker = new ConcurrentProcessing(Expression,queue,corr,"sigmoid",mu,alpha);
			Future<HashMap<String,float[]>> submit = completionService.submit(worker);
			taskList.add(submit);  
		}
		 
		// collect results
		for(int t=0;t<Threads;t++){
			try{
				HashMap<String,float[]> hm = completionService.take().get();
				//System.err.println("obtained result for thread " + t);
				int r=0;
				for(Map.Entry<String,float[]> entry : hm.entrySet()){
										String s = entry.getKey();
//					String[] S = s.split("-");
										int i = Integer.valueOf(s);
										int size = D-i;
                                        float[] d = new float[size];
                                        d = entry.getValue();
                                        for(int j=0;j<d.length;j++){
                                            int coord = j+i;
                                            Adjacency.setValueByEntry(d[j],i,coord);
					//System.out.println(i+"\t"+j+"\t"+d);
					//Adjacency.setValueByEntry( d,i,j);
					//Adjacency.setValueByEntry( d,j,i);
                                            r++;
                                        }
					
				}
				//System.err.println("Processed "+ r + " records");
			}catch(InterruptedException e){
				e.printStackTrace();
			}catch (ExecutionException e){
				e.printStackTrace();
			}
			//System.err.println("Thread " + t + " complete.");
		}
		//System.err.println("Done.");
		pool.shutdownNow();
		return Adjacency;
	}
	
	/*
	public static GCNMatrix calculateSigmoidAdjacency (GCNMatrix Similarity,float mu, float alpha, int Threads){
		int D = Similarity.getNumRows();
		GCNMatrix Adjacency = new GCNMatrix(D,D);
                Adjacency = Operations.copyNames(Similarity, Adjacency);
		ExecutorService pool = Executors.newFixedThreadPool(Threads);
		ExecutorCompletionService<HashMap<String,Float>> completionService = new ExecutorCompletionService<>(pool);
		List<Future<HashMap<String,Float>>> taskList = new ArrayList<Future<HashMap<String,Float>>>();
		ConcurrentLinkedQueue<String> queue = new ConcurrentLinkedQueue<String>();
		System.err.println("Processing adjacency using " + Threads + " threads.");
		
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				String S = i+"-"+j;
				queue.add(S);
			}
		}
		
		for ( int i = 0; i < Threads; i++ ) {
			Callable<HashMap<String,Float>> worker = new ConcurrentProcessing(Similarity,queue,"sigmoid",mu,alpha);
			Future<HashMap<String,Float>> submit = completionService.submit(worker);
			taskList.add(submit);  
		}
		
	
		for(int t=0;t<Threads;t++){
			try{
				HashMap<String,Float> hm = completionService.take().get();
				System.err.println("obtained result for thread " + t);
				int r=0;
				for(Map.Entry<String,Float> entry : hm.entrySet()){
					String s = entry.getKey();
					String[] S = s.split("-");
					Float d = entry.getValue();
					int i = Integer.parseInt(S[0]);
					int j = Integer.parseInt(S[1]);
					//System.out.println(i+"\t"+j+"\t"+d);
					Adjacency.setValueByEntry( d,i,j);
					Adjacency.setValueByEntry( d,j,i);
					r++;
				}
				System.err.println("Processed "+ r + " records");
			}catch(InterruptedException e){
				e.printStackTrace();
			}catch (ExecutionException e){
				e.printStackTrace();
			}
			System.err.println("Thread " + t + " complete.");
		}
		System.err.println("Done.");
		pool.shutdown();
		return Adjacency;
	}
	*/
	public static GCNMatrix calculateDifference (GCNMatrix mat1, GCNMatrix mat2){
		int D = mat1.getNumRows();
		GCNMatrix Difference = new GCNMatrix(D,D);
                Difference = Operations.copyNames(mat1.getRowNames(), Difference);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				float v1=mat1.getValueByEntry(i, j);
				float v2=mat2.getValueByEntry(i, j);
				//if((v1 != 0) & (v2 != 0)){
				float d1 =v1-v2;
//                                System.out.println("Val1: " + v1 +" Val2: " + v2 + " diff " + d1);
				Difference.setValueByEntry(d1,i,j);
				//}
			}
		}
		return Difference;
	}
		
      
public static void generateHistogramHM (GCNMatrix DataFrame, String pathOut, String Title,String Xlab, String Ylab,boolean print) {
	int H = DataFrame.getNumRows();
	int W = DataFrame.getNumColumns();
	DecimalFormat df = new DecimalFormat("#.##");
	df.setRoundingMode(RoundingMode.HALF_UP);
	TreeMap<Float,Integer> HMHistogram = new TreeMap<Float,Integer>();
	for(int i=0;i<H;i++){
		for(int j=0;j<W;j++){
			//System.out.println("Val: "+DataFrame[i][j]+"\n");
			if(DataFrame.getValueByEntry(i,j) != 0){
				try{
					Float V = (Float.valueOf(df.format(DataFrame.getValueByEntry(i,j))));
                                        
					if(HMHistogram.containsKey(V)){
						Integer I = HMHistogram.get(V);
                                               // System.out.println("Putting " + V+ " and " + I);
						HMHistogram.put(V,I+1);
					}else{
						HMHistogram.put(V,1);
                                               // System.out.println("Putting " + V+ " and 1");
					}
				}catch(NumberFormatException ex){
					System.out.println("Obtain " + DataFrame.getValueByEntry(i,j) +" from matrix.");
					System.exit(1);
				}
			}
		}
	}
	XYSeriesCollection dataset = new XYSeriesCollection();
	XYSeries series = new XYSeries("Values");
	for(Map.Entry<Float,Integer> entry : HMHistogram.entrySet()) {
			  Float A;
                          A = entry.getKey();
			  Integer value;
                          value = entry.getValue();
                          //System.err.println("Got " + A + " and " + value);
                          if(print == true){
                            System.out.println(A +","+value);
                          }
			  series.add((double) A, (double) value,true);
	}
	dataset.addSeries(series);
	JFreeChart chart = ChartFactory.createXYLineChart(Title,Xlab,Ylab,dataset, PlotOrientation.VERTICAL, 
			 false, true, false);
	chart.setBackgroundPaint(Color.white);
	chart.setAntiAlias(true);	
	final Plot plot = chart.getPlot();
	plot.setBackgroundPaint(Color.white);
	plot.setOutlinePaint(Color.black);
	try {
		ChartUtilities.saveChartAsJPEG(new File(pathOut), chart, 500, 300);
	} catch (IOException e) {
		System.err.println("Problem occurred creating chart.");
	}
}

	public static float permuteData(ExpressionFrame expF1, ExpressionFrame expF2, int P,String out,int threads) {
		int s1 = expF1.getNumColumns();
		int s2 = expF2.getNumColumns();
		int R = expF1.getNumRows();
		int M = Math.min(expF1.getNumColumns(), expF2.getNumColumns());
		int S = s1+s2;
		int s = M * 2;
		Integer[][] Sets = new Integer[P][];
		Sets = _getPermutations(s1,s2,P);
		//System.out.println("size1: "+s1+"\nsize2: "+s2);
		// for all i in Sets p, obtain values of i<s from exp1, values of i>s
                
                float CUTOFF = 1.0f;
                
                ArrayList<TreeMap<Float,Integer>> Perms = new ArrayList<>();
		for(int p=0;p<P;p++){
			System.out.println("Permutation: " + p);
			ExpressionFrame pF1 = new ExpressionFrame(R,M);
			ExpressionFrame pF2 = new ExpressionFrame(R,M);
			for(int r=0;r<R;r++){
				float[] rF1 = expF1.getRowByIndex(r);
				float[] rF2 = expF2.getRowByIndex(r);
				float[] nR1 = new float[M];
				float[] nR2 = new float[M];
				for(int m=0;m<M;m++){
					int ind = Sets[p][m];
					//System.out.println("getting "+ ind + " ("+ m +")");
					if(ind<s1){
						//System.out.println("getting "+ind+" from 1 s1 is " + s1);
						nR1[m]=rF1[ind];
					}else if(ind>=s1){
						ind = ind - s1;
					//	System.out.println("getting "+ind+" from 2 s1 is " + s1);
						nR1[m]=rF2[ind];
					}
				}
				pF1.addRow(nR1);
				for(int m=M;m<s;m++){
					int Mind = m-M;
					int ind = Sets[p][m];
					//System.out.println("getting "+ ind + " ("+ m +")");
					if(ind<s1){
						//System.out.println("getting "+ind+" from 1 s1 is " + s1);
						nR2[Mind]=rF1[ind];
					}else if(ind>=s1){
						ind = ind - s1;
						//System.out.println("getting "+ind+" from 2 s1 is " + s1);
						nR2[Mind]=rF2[ind];
					}
				}
				pF2.addRow(nR2);
			}	
			GCNMatrix CurrentMatrix1 = Operations.calculateAdjacency(pF1,"pcc","sigmoid",0.8f,20.0f,threads);
			GCNMatrix CurrentMatrix2 = Operations.calculateAdjacency(pF2,"pcc","sigmoid",0.8f,20.0f,threads);
			CurrentMatrix1.calculateKs();
			CurrentMatrix2.calculateKs();
                        CurrentMatrix1 = Operations.calculateTOM(CurrentMatrix1,16);
                        CurrentMatrix2 = Operations.calculateTOM(CurrentMatrix2,16);
			GCNMatrix Difference = Operations.calculateDifference(CurrentMatrix1,CurrentMatrix2);
                        TreeMap<Float,Integer> Distribution = Difference.generateDistribution();
                        Perms.add(Distribution);
                        Difference.maskMatrix(0.02f);
			String O2 = out + "/Selfwise."+ p + ".testing.jpeg";
                        Operations.generateHistogramHM(Difference, O2, "Cross-network Adjacency diffs Zm vs Sv", "selfwise TOM", "Count", false);
                        /*CurrentMatrix1 = Operations.calculateTOM(CurrentMatrix1, threads);
                        CurrentMatrix2 = Operations.calculateTOM(CurrentMatrix2, threads);
                        Difference = Operations.calculateDifference(CurrentMatrix1,CurrentMatrix2);
                        String O1 = out + "/Pairwise." + p + ".jpeg";
                        Difference.maskMatrix(0.02f);
                        Operations.generateHistogramHM(Difference, O1, "Pairwise Adjacency Differences Zm vs Sv", "cross-pair Delta-Adj", "Count", false);
                        */
		}
                GCNMatrix NetworkA = Operations.calculateAdjacency(expF1,"pcc","sigmoid",0.8f,20.0f,16);
                GCNMatrix NetworkB = Operations.calculateAdjacency(expF2,"pcc","sigmoid",0.8f,20.0f,16);
                NetworkA.calculateKs();
                NetworkB.calculateKs();
                NetworkA = Operations.calculateTOM(NetworkA,16);
                NetworkB = Operations.calculateTOM(NetworkB,16);
                GCNMatrix rDiff = Operations.calculateDifference(NetworkA, NetworkB);
                TreeMap<Float,Integer> Real = rDiff.generateDistribution();
                
                System.err.println("Cutoff\tAverage False\tTrue\tFDR");
                for(float c=0.0f;c<1.0f;c+=0.01f){
                    float C = c;
                    Double Total=0.0d;
                    for(int a=0;a<Perms.size();a++){
                        for(Map.Entry<Float,Integer> entry : Perms.get(a).entrySet()) {
                            Float A;
                            A = entry.getKey();
                            Double value;
                            value = Double.valueOf(entry.getValue());
                            if(Math.abs(A)>=C){
                                Total += value;
                            }
                        }
                    }
                    // Total holds all instances of adj value > C across all perms
                    Double Average = Total/Perms.size();
                    Double RealHits = 0.0d;
                    for(Map.Entry<Float,Integer> entry : Real.entrySet()) {
                        Float A;
                        A = entry.getKey();
                        Double value;
                        value = Double.valueOf(entry.getValue());
                        if(Math.abs(A)>=C){
                            RealHits += value;
                        }
                    }
                    double FDR = Average/RealHits;
                    if(FDR<=0.15){
                        System.err.println(C + "\t" + Average + "\t" + RealHits + "\t" + FDR);
                    }
                    
                    if(FDR <= 0.05){
                        if(CUTOFF==1){
                            CUTOFF = C;
                        }else{
                            
                        }
                    }
                    
                }
		// Now have permuted expression frames?
		// NEEDS TESTING.
                return CUTOFF;
	}
        
	public static void permuteDataHalf(ExpressionFrame expF1, ExpressionFrame expF2, int P,String out,int threads) {
		int s1 = expF1.getNumColumns();
		int s2 = expF2.getNumColumns();
		int R = expF1.getNumRows();
		int M = Math.min(expF1.getNumColumns(), expF2.getNumColumns());
		int S = s1+s2;
		int s = M * 2;
		Integer[][] Sets = new Integer[P][];
		Sets = _getPermutations(s1,s2,P);
		//System.out.println("size1: "+s1+"\nsize2: "+s2);
		// for all i in Sets p, obtain values of i<s from exp1, values of i>s
                GCNMatrix RealAdj = Operations.calculateAdjacency(expF1,"pcc","sigmoid",0.8f,16.0f,threads);
                RealAdj.calculateKs();
                GCNMatrix RealTom = Operations.calculateTOM(RealAdj, threads);
		for(int p=0;p<P;p++){
			System.out.println("Permutation: " + p);
			ExpressionFrame pF1 = new ExpressionFrame(R,M);
			for(int r=0;r<R;r++){
				float[] rF1 = expF1.getRowByIndex(r);
				float[] rF2 = expF2.getRowByIndex(r);
				float[] nR1 = new float[M];
				float[] nR2 = new float[M];
				for(int m=0;m<M;m++){
					int ind = Sets[p][m];
					//System.out.println("getting "+ ind + " ("+ m +")");
					if(ind<s1){
						//System.out.println("getting "+ind+" from 1 s1 is " + s1);
						nR1[m]=rF1[ind];
					}else if(ind>=s1){
						ind = ind - s1;
					//	System.out.println("getting "+ind+" from 2 s1 is " + s1);
						nR1[m]=rF2[ind];
					}
				}
				pF1.addRow(nR1);
			}
                       	GCNMatrix CurrentMatrix = Operations.calculateAdjacency(pF1,"pcc","sigmoid",0.8f,16.0f,threads);
                        
			CurrentMatrix.calculateKs();
                        
			GCNMatrix Difference = Operations.calculateDifference(RealAdj,CurrentMatrix);
                        Difference.maskMatrix(0.02f);
			String O2 = out + "/Selfwise."+ p + ".testing.jpeg";
                        Operations.generateHistogramHM(Difference, O2, "Cross-network Selfwise Topological Overlap Zm vs Sv", "selfwise TOM", "Count", false);
                        
                        /*
                        CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);
                        
                        Difference = Operations.calculateDifference(RealTom,CurrentMatrix);
                        String O1 = out + "/Pairwise." + p + ".jpeg";
                        Difference.maskMatrix(0.02f);
                        Operations.generateHistogramHM(Difference, O1, "Pairwise Adjacency Differences Zm vs Sv", "cross-pair Delta-Adj", "Count", false);
                        */
		}
		// Now have permuted expression frames?
		// NEEDS TESTING.
	}
        
	private static Integer[][] _getPermutations(int s1,int s2,int p) {
		Integer[][] Sets = new Integer[p][];
		int m = Math.min(s1,s2);
		int S = s1+s2;
		int s = m * 2;
		Integer ind[] = new Integer[S];
		for(int i=0;i<S;i++){
			ind[i]=i;
		}
		for(int i=0;i<p;i++){
			
			// obtain a array of random order indicies, i-> [0,123,2,5,3,1]
			Integer shuffled[] = fisherYates(ind);
			Integer[] perms = Arrays.copyOfRange(shuffled, 0, s);
			/*for(int x=0;x<perms.length;x++){
				System.out.print(perms[x]+",");
			}
			System.out.print("\n");*/
			Sets[i]=perms;
		}
		//System.exit(0);
		return Sets;
	}
	
	/// shuffle some primitives
	private static Integer[] fisherYates (Integer[] array){
	    Random rnd = new Random();
	    for (int i = array.length - 1; i > 0; i--){
	      int index = rnd.nextInt(i + 1);
	      // Simple swap
	      int a = array[index];
	      array[index] = array[i];
	      array[i] = a;
	    }
	    return array;
	  }
	public static void createDirectory (String dir){
		File theDir = new File(dir);
        if (!theDir.exists()) {
            System.out.println("creating directory: " + dir);
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                //TODO handle it
            }
        }
	}
	
}



