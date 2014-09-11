package networkCalcPackage;

import java.util.HashMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

//public class SimilarityConcurrent implements Runnable {
public class ConcurrentProcessing implements Callable<HashMap<String,float[]>> {
	private ConcurrentLinkedQueue<String> queue;
	private String smethod;
	private String amethod;
	private ExpressionFrame Exp;
	private GCNMatrix Adj;
	private int D;
	private float M;
	private float A;
	//private ThreadLocal<HashMap<String, Double>> hm = new ThreadLocal<HashMap<String, Double>>();
	
	private void noCall () {
		System.err.println("Threading called with no method. Should not happen");
		System.exit(1);
	}
	
	public HashMap<String, float[]> call() {
		HashMap<String, float[]> hm = new HashMap<String, float[]>();
		 	//try {
		String s = null;
	        while ( ( s=queue.poll() ) != null ) {
	            //while ( true ) {
	        	int L = Integer.valueOf(s);
	        	int Size;
                        Size = D - L;
	        	float value[] = new float[Size];
	                //String s = queue.remove();
	                switch (smethod) {
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
	                /*for(int j=0;j<Size;j++){
	                	int coord = j+L;
	                	String k = L + "-" + coord;
	                	hm.put(k,value[j]);
	                }*/
	        }
	       // }
	        /*catch ( InterruptedException ie ) { 
	            // just terminate
	        }*/
		return hm;
	}
	public ConcurrentProcessing (GCNMatrix Adjacency,ConcurrentLinkedQueue<String> queue, String corr) {
        this.queue = queue;
        this.smethod = corr;
        this.Adj = Adjacency;
        this.D=Adj.getNumRows();
    }
	
	public ConcurrentProcessing (ExpressionFrame expression,ConcurrentLinkedQueue<String> queue, String corr,String adj, float mu, float a) {
        this.queue = queue;
        this.smethod = corr;
        this.amethod = adj;
        this.Exp = expression;
        this.M = mu;
        this.A = a;
        this.D=Exp.getNumRows();
    }
	
	public ConcurrentProcessing (GCNMatrix Similarity,ConcurrentLinkedQueue<String> queue, String Method,float mu, float a) {
        this.queue = queue;
        this.smethod = Method;
        this.Adj = Similarity;
        this.D=Adj.getNumRows();
       
    }
	
	public float[] doWork_gini (String s){
		int i = Integer.parseInt(s);
		int size = Exp.getNumRows() - i;
		float[] GINIS = new float[size];
		float[] I_data = Exp.getRowByIndex(i);
		int[] I_ranks = Operations.getIndicesInOrder(I_data);
		for(int j=i;j<Exp.getNumRows();j++){
			
			float gcc=(float) 0.0;
			int coord = j-i;
			
			if(i==j){
				gcc=1.0f;
			}else{
				float[] J_data = Exp.getRowByIndex(j);
				int[] J_ranks = Operations.getIndicesInOrder(J_data);
				float I_num=0.0f;
				float J_num=0.0f;
				float GCC1=0.0f;
				float GCC2=0.0f;
				for(int x=0;x<J_data.length;x++){
					I_num += ((2*(i+1))-I_data.length-1) * I_data[J_ranks[i]];
					J_num += ((2*(i+1))-J_data.length-1) * J_data[I_ranks[i]];
				}
				GCC1=I_num/Exp.getGiniDenom(j);
				GCC2=J_num/Exp.getGiniDenom(i);
				if(Math.abs(GCC1)>Math.abs(GCC2)){
					gcc=GCC1;
				}else if (Math.abs(GCC2)>Math.abs(GCC1)){
					gcc=GCC2;
				}else{
					gcc=GCC1;
				}
				
			}
			GINIS[coord]=_getSigmoid(gcc);
			
		}
		return GINIS;
	}
	
	public float[] doWork_pcc (String s){
		int i = Integer.parseInt(s);
		int size = Exp.getNumRows() - i;
		float[] means = new float[size];
		float[] I_data = Exp.getRowByIndex(i);
		float I_mean = Exp.getMean(i);
        	for(int j=i;j<Exp.getNumRows();j++){
			float correlation = 0.0f;
			int coord = j-i;
			if(i==j){
				correlation = 1.0f;
			}else{
				float[] J_data = Exp.getRowByIndex(j);
				float J_mean = Exp.getMean(j);
				double SQR1=0.0;
				double SQR2=0.0;
				double Na=0.0;
				for(int n=0;n<J_data.length;n++){
					double v1 = (double) I_data[n]-I_mean;
					double v2 = (double) J_data[n]-J_mean;
					Na += (v1 * v2);
					SQR1 += (v1 * v1);
					SQR2 += (v2 * v2);
				}
				correlation = (float) (Na/(Math.sqrt(SQR1) * Math.sqrt(SQR2)));
                                correlation = _getSigmoid(correlation);
			}
                        means[coord]=_getSigmoid(correlation);
		}
		return means;
	}
	
	private float _getSigmoid (float V){
		return (float) (1.0f/(1+Math.exp(A*-1*(Math.abs(V)-M))));
	}
	
	public float[] doWork_sigmoid (String s){
		int i = Integer.parseInt(s);
		int size = Exp.getNumRows() - i;
		float[] adjacency = new float[size];
		for(int j=i;j<Exp.getNumRows();j++){
			int coord = j-i;
			if(i==j){
				adjacency[coord] = 1.0f;
			}else{
				adjacency[coord] = _getSigmoid(Exp.getValueByEntry(i,j));
			}
		}
		return adjacency;
	}
	
	public float[] doWork_tom (String s){
		int i = Integer.parseInt(s);
		int size = D - i;
		float i_k = Adj.findK(i, i);
		float[] TOM = new float[size];
		for(int j=i;j<D;j++){
			int coord = j-i;
			float tom = 0.0f;
			if(i==j){
				tom = 1.0f;
			}else{
				float product=0;
				float j_k = Adj.findK(j,j);
				for(int u=0;u<D;u++){
					if((u != i) && (u != j) && (Adj.testValue(i, u)) && (Adj.testValue(j, u))){
						product += Adj.getValueByEntry(i,u) * Adj.getValueByEntry(j,u);
					}
				}
				float k_min = Math.min(i_k, j_k);
				float DFIJ=Adj.getValueByEntry(i,j);
				tom = ((product+DFIJ)/(k_min + 1 - DFIJ));
			}
			TOM[coord]=tom;
		}
		return TOM;
	}
	
}


