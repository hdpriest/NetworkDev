package networkCalcPackage;

import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

//public class SimilarityConcurrent implements Runnable {
public class ConcurrentProcessing implements Callable<HashMap<String,Float>> {
	private ConcurrentLinkedQueue<String> queue;
	private String m;
	private GCNMatrix Exp;
	private int D;
	private double M;
	private double A;
	//private ThreadLocal<HashMap<String, Double>> hm = new ThreadLocal<HashMap<String, Double>>();
	
	private void noCall () {
		System.err.println("Threading called with no method. Should not happen");
		System.exit(1);
	}
	
	public HashMap<String, Float> call() {
			HashMap<String, Float> hm = new HashMap<String, Float>();
		 	//try {
			String s = null;
	        	while ( ( s=queue.poll() ) != null ) {
	            //while ( true ) {
	        	float value = 0.0f;
	                //String s = queue.remove();
	                switch (m) {
	                	case "gini": value = doWork_gini(s);
	                	break;
	                	case "pcc": value = doWork_pcc(s);
	                	break;
	                	case "sigmoid" : value = doWork_sigmoid(s);
	                	break;
	                	case "tom" : value = doWork_tom(s);
	                	break;
	                	default: noCall();
	                	break;
	                }
	                hm.put(s, value);
	            }
	       // }
	        /*catch ( InterruptedException ie ) { 
	            // just terminate
	        }*/
		 	return hm;
	}
	
	public ConcurrentProcessing (GCNMatrix Expression,ConcurrentLinkedQueue<String> queue, String Method) {
        this.queue = queue;
        this.m = Method;
        this.Exp = Expression;
        this.M=0.0; // little hackery
        this.A=0.0;
        this.D=Exp.getNumRows();
    }
	
	public ConcurrentProcessing (GCNMatrix Similarity,ConcurrentLinkedQueue<String> queue, String Method,double mu, double a) {
        this.queue = queue;
        this.m = Method;
        this.Exp = Similarity;
        this.M = mu;
        this.A = a;
        this.D=Exp.getNumRows();
       
    }
	
	public float doWork_gini (String s){
		String[] S = s.split("-");
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		float GCC1=(float) 0.0;
		float GCC2=(float) 0.0;
		float gcc=(float) 0.0;
		if(i==j){
			GCC1 = (float) 1.0;
			GCC2 = (float) 0.0;
		}else{
			float[] I_data = Exp.getRowByIndex(i);
			float[] J_data = Exp.getRowByIndex(j);
			GCC1 = Operations.GINI(I_data,J_data);
			GCC2 = Operations.GINI(J_data,I_data);
		}
		if(Math.abs(GCC1)>Math.abs(GCC2)){
			gcc=GCC1;
		}else if (Math.abs(GCC2)>Math.abs(GCC1)){
			gcc=GCC2;
		}else{
			gcc=GCC1;
		}
		return gcc;
	}
	
	public float doWork_pcc (String s){
		String[] S = s.split("-");
		float correlation = (float) 0.0;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		double[] I_data = Exp.getRowByIndexDbl(i);
		double[] J_data = Exp.getRowByIndexDbl(j);
		if(i==j){
			correlation = 1.0f;
		}else{
			PearsonsCorrelation corr = new PearsonsCorrelation();
			correlation = (float) corr.correlation(I_data, J_data);
		}
		return correlation;
	}
	
	public float doWork_sigmoid (String s){
		String[] S = s.split("-");
		float adjacency = (float) 0.0;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		if(i==j){
			adjacency = (float) 1.0;
		}else{
			adjacency = (float) (1/(1+Math.exp(A*-1*(Math.abs(Exp.getValueByEntry(i,j))-M))));
		}
		return adjacency;
	}
	
	public float doWork_tom (String s){
		String[] S = s.split("-");
		float tom = 0.0f;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		if(i==j){
			tom = 1.0f;
		}else{
			float product=0;
			float i_k = Exp.findK(i, i);
			float j_k = Exp.findK(j,j);
			for(int u=0;u<D;u++){
				if((u != i) && (u != j) && (Exp.testValue(i, u)) && (Exp.testValue(j, u))){
					product += Exp.getValueByEntry(i,u) * Exp.getValueByEntry(j,u);
				}
			}
			float k_min = Math.min(i_k, j_k);
			float DFIJ=Exp.getValueByEntry(i,j);
			tom = ((product+DFIJ)/(k_min + 1 - DFIJ));
		}
		return tom;
	}
	
}


