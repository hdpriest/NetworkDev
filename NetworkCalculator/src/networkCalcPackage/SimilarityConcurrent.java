package networkCalcPackage;

import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

//public class SimilarityConcurrent implements Runnable {
public class SimilarityConcurrent implements Callable<HashMap<String,Double>> {
	private final BlockingQueue<String> queue;
	private final String m;
	private final GCNMatrix Exp;
	//private ThreadLocal<HashMap<String, Double>> hm = new ThreadLocal<HashMap<String, Double>>();
	
	public HashMap<String, Double> call() {
			HashMap<String, Double> hm = new HashMap<String, Double>();
		 	try {
	        	while ( queue.peek() != null ) {
	            //while ( true ) {
	                String s = queue.take();
	                switch (m) {
	                	case "gini": hm = doWork_gini(s);
	                	break;
	                	case "pcc": hm = doWork_pcc(s);
	                	break;
	                	default: hm = doWork_pcc(s);
	                	break;
	                }
	            }
	        }
	        catch ( InterruptedException ie ) { 
	            // just terminate
	        }
		 	return hm;
	}
	
	public SimilarityConcurrent (GCNMatrix Expression,BlockingQueue<String> queue, String Method) {
        this.queue = queue;
        this.m = Method;
        this.Exp = Expression;
    }
	
	/*public void run() {
        try {
        	while ( queue.peek() != null ) {
            //while ( true ) {
                String s = queue.take();
                switch (m) {
                	case "gini": doWork_gini(s);
                	break;
                	case "pcc": doWork_pcc(s);
                	break;
                	default: doWork_pcc(s);
                	break;
                }
            }
        }
        catch ( InterruptedException ie ) { 
            // just terminate
        }
    }
	*/
	public HashMap<String,Double> doWork_gini (String s){
		String[] S = s.split("-");
		HashMap <String,Double> hm = new HashMap<String, Double>(); 
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		double GCC1=0.0;
		double GCC2=0.0;
		if(i==j){
			GCC1 = 1.0;
			GCC2 = 0.0;
		}else{
			double[] I_data = Exp.getRowByIndex(i);
			double[] J_data = Exp.getRowByIndex(j);
			GCC1 = Operations.GINI(I_data,J_data);
			GCC2 = Operations.GINI(J_data,I_data);
		}
		if(Math.abs(GCC1)<Math.abs(GCC2)){
			hm.put(s,(Double) GCC1);
		}else if (Math.abs(GCC2)<Math.abs(GCC1)){
			hm.put(s,(Double) GCC2);
		}else{
			hm.put(s,(Double) GCC1);
		}
		return hm;
	}
	
	public HashMap<String,Double> doWork_pcc (String s){
		String[] S = s.split("-");
		HashMap <String,Double> hm = new HashMap<String, Double>();
		int i = Integer.parseInt(S[0]);
		int j = Integer.parseInt(S[1]);
		double[] I_data = Exp.getRowByIndex(i);
		double[] J_data = Exp.getRowByIndex(j);
		if(i==j){
			hm.put(s,1.0);
		}else{
			PearsonsCorrelation corr = new PearsonsCorrelation();
			double correlation = corr.correlation(I_data,J_data);
			hm.put(s,(Double) correlation);
		}
		return hm;
	}
	
}


