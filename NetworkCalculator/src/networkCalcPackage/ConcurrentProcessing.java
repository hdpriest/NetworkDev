package networkCalcPackage;

import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

//public class SimilarityConcurrent implements Runnable {
public class ConcurrentProcessing implements Callable<HashMap<String,Double>> {
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
	
	public HashMap<String, Double> call() {
			HashMap<String, Double> hm = new HashMap<String, Double>();
		 	//try {
			String s = null;
	        	while ( ( s=queue.poll() ) != null ) {
	            //while ( true ) {
	        		double value = 0.0;
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
	
	public double doWork_gini (String s){
		String[] S = s.split("-");
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		double GCC1=0.0;
		double GCC2=0.0;
		double gcc=0.0;
		if(i==j){
			GCC1 = 1.0;
			GCC2 = 0.0;
		}else{
			double[] I_data = Exp.getRowByIndex(i);
			double[] J_data = Exp.getRowByIndex(j);
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
	
	public double doWork_pcc (String s){
		String[] S = s.split("-");
		double correlation = 0.0;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		double[] I_data = Exp.getRowByIndex(i);
		double[] J_data = Exp.getRowByIndex(j);
		if(i==j){
			correlation = 1.0;
		}else{
			PearsonsCorrelation corr = new PearsonsCorrelation();
			correlation = corr.correlation(I_data,J_data);
		}
		return correlation;
	}
	
	public double doWork_sigmoid (String s){
		String[] S = s.split("-");
		double adjacency = 0.0;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		if(i==j){
			adjacency = 1.0;
		}else{
			adjacency = 1/(1+Math.exp(A*-1*(Math.abs(Exp.getValueByEntry(i,j))-M)));
		}
		return adjacency;
	}
	
	public double doWork_tom (String s){
		String[] S = s.split("-");
		double tom = 0.0;
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		if(i==j){
			tom = 1.0;
		}else{
			double product=0;
			double i_k = Exp.findK(i, i);
			double j_k = Exp.findK(j,j);
			for(int u=0;u<D;u++){
				if((u != i) && (u != j) && (Exp.testValue(i, u)) && (Exp.testValue(j, u))){
					product += Exp.getValueByEntry(i,u) * Exp.getValueByEntry(j,u);
				}
			}
			double k_min = Math.min(i_k, j_k);
			double DFIJ=Exp.getValueByEntry(i,j);
			tom=(product+DFIJ)/(k_min + 1 - DFIJ);
		}
		return tom;
	}
	
}


