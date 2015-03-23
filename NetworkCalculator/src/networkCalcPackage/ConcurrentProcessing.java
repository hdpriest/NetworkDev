package networkCalcPackage;

import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;

public class ConcurrentProcessing implements Callable<HashMap<String, float[]>> {

    private ConcurrentLinkedQueue<String> queue;
    private String smethod;
    private String amethod;
    private ExpressionFrame Exp;
    private GCNMatrix Adj;
    private GCNMatrix Adj_b;
    private int D;
    private float M;
    private float A;

    private void noCall() {
        System.err.println("Threading called with no method. Something went very wrong.");
        System.exit(1);
    }

    public HashMap<String, float[]> call() {
        HashMap<String, float[]> hm = new HashMap<String, float[]>();
        String s = null;
        while ((s = queue.poll()) != null) {
            double proc = queue.size() * .001;
            double IntProc = java.lang.Math.ceil(proc);
            if((IntProc == proc) && (proc != 0.0f)){
                int rec = (int) proc * 1000 * D;
                System.err.println(rec+ " records remaining... ");
            }
            int L = Integer.valueOf(s);
            int Size;
            Size = D - L;
            float value[] = new float[Size];
            switch (smethod) {
                case "gini":
                    value = doWork_gini(s);
                    break;
                case "pcc":
                    value = doWork_pcc(s);
                    break;
                case "spearman":
                    value = doWork_spearman(s);
                    break;
                case "sigmoid":
                    value = doWork_sigmoid(s);
                    break;
                case "tom":
                    value = doWork_tom(s);
                    break;
                case "difference":
                    value = doWork_diff(s);
                    break;
                default:
                    noCall();
                    break;
            }
            hm.put(s, value);
        }
        return hm;
    }
    
    public ConcurrentProcessing(GCNMatrix Adjacency_a, GCNMatrix Adjacency_b, ConcurrentLinkedQueue<String> queue, String corr) {
        this.queue = queue;
        this.smethod = corr;
        this.Adj = Adjacency_a;
        this.Adj_b = Adjacency_b;
        this.D = Adj.getNumRows();
    }

    public ConcurrentProcessing(GCNMatrix Adjacency, ConcurrentLinkedQueue<String> queue, String corr) {
        this.queue = queue;
        this.smethod = corr;
        this.Adj = Adjacency;
        this.D = Adj.getNumRows();
    }

    public ConcurrentProcessing(ExpressionFrame expression, ConcurrentLinkedQueue<String> queue, String corr, String adj, float mu, float a) {
        this.queue = queue;
        this.smethod = corr;
        this.amethod = adj;
        this.Exp = expression;
        this.M = mu;
        this.A = a;
        this.D = Exp.getNumRows();
    }

    public ConcurrentProcessing(GCNMatrix Similarity, ConcurrentLinkedQueue<String> queue, String Method, float mu, float a) {
        this.queue = queue;
        this.smethod = Method;
        this.M = mu;
        this.A = a;
        this.Adj = Similarity;
        this.D = Adj.getNumRows();

    }
    private boolean _checkData (float[] Data){
        // We preclude loading of all null/zero data rows
        // However, permutation of datasets could result in a all-zero expression set
        // in this case, return 0, no correlation
        // Also need to handle no-variance dataset (all 1s, all anything, this will throw many correlations for a loop)
        float max = 0.0f;
        float min = Float.POSITIVE_INFINITY;
        for(int i=0;i<Data.length;i++){
            if(Data[i] > max){
                max = Data[i];    
            }
            if(Data[i] < min){
                min = Data[i];
            }
            
        }
        if(max == min) return false;
        return true;
    }
    
    public float[] doWork_spearman(String s) {
        int i = Integer.parseInt(s);
        int size = Exp.getNumRows() - i;
        float[] correlations = new float[size];
        float[] I_data = Exp.getRowByIndex(i);
        float DataSize = (float) I_data.length;
        int[] I_ranks = Operations.getRanks(I_data);
        if(_checkData(I_data)){
            for (int j = i; j < Exp.getNumRows(); j++) {

                float correlation = 0.0f;
                int coord = j - i;

                if (i == j) {
                    correlation = 1.0f;
                } else {
                    float[] J_data = Exp.getRowByIndex(j);
                    if(_checkData(J_data)){
                        int[] J_ranks = Operations.getRanks(J_data);
                        float sum=0.0f;
                        for(int X=0;X<J_ranks.length;X++){
                            float diff = (float) (I_ranks[X]-J_ranks[X]);
                            diff = diff*diff;
                            sum += diff;
                        }
                        float p = (float) 1.0f-( (6.0f*sum) / (DataSize* ( (DataSize*DataSize)-1 ) ) );
                        correlation=p;
                    }
                }
                correlations[coord] = _getSigmoid(correlation);
            }
        }else{
            for (int j = i; j < Exp.getNumRows(); j++) {
                int coord = j - i;
                if(j == i){
                    correlations[coord] = 1.0f;
                }else{
                    correlations[coord] = 0.0f;    
                }
            }
        }
        return correlations;
    }
    
    public float[] doWork_gini(String s) {
        int i = Integer.parseInt(s);
        int size = Exp.getNumRows() - i;
        float[] GINIS = new float[size];
        float[] I_data = Exp.getRowByIndex(i);
        int[] I_ranks = Operations.getIndicesInOrder(I_data);
        if(_checkData(I_data)){
            for (int j = i; j < Exp.getNumRows(); j++) {

                float gcc = 0.0f;
                int coord = j - i;

                if (i == j) {
                    gcc = 1.0f;
                } else {
                    float[] J_data = Exp.getRowByIndex(j);
                    if(_checkData(J_data)){
                        int[] J_ranks = Operations.getIndicesInOrder(J_data);
                        float I_num = 0.0f;
                        float J_num = 0.0f;
                        float GCC1 = 0.0f;
                        float GCC2 = 0.0f;
                        for (int x = 0; x < J_data.length; x++) {
                            I_num += ((2 * (x + 1)) - I_data.length - 1) * I_data[J_ranks[x]];
                            J_num += ((2 * (x + 1)) - J_data.length - 1) * J_data[I_ranks[x]];
                        }
                        GCC1 = I_num / Exp.getGiniDenom(j);
                        GCC2 = J_num / Exp.getGiniDenom(i);
                        if (Math.abs(GCC1) < Math.abs(GCC2)) {
                            gcc = GCC1;
                        } else if (Math.abs(GCC2) < Math.abs(GCC1)) {
                            gcc = GCC2;
                        } else {
                            gcc = GCC1;
                        }
                    }
                }
                GINIS[coord] = _getSigmoid(gcc);

            }
        }else{
            for (int j = i; j < Exp.getNumRows(); j++) {
                int coord = j - i;
                if(j == i){
                    GINIS[coord] = 1.0f;
                }else{
                    GINIS[coord] = 0.0f;    
                }
            }
        }
        return GINIS;
    }

    public float[] doWork_pcc(String s) {
        int i = Integer.parseInt(s);
        int size = Exp.getNumRows() - i;
        float[] correlations = new float[size];
        float[] I_data = Exp.getRowByIndex(i);
        if(_checkData(I_data)){
            float I_mean = Exp.getMean(i);
            for (int j = i; j < Exp.getNumRows(); j++) {
                float correlation = 0.0f;
                int coord = j - i;
                if (i == j) {
                    correlation = 1.0f;
                } else {
                    float[] J_data = Exp.getRowByIndex(j);
                    if(_checkData(J_data)){
                        float J_mean = Exp.getMean(j);
                        double SQR1 = 0.0;
                        double SQR2 = 0.0;
                        double Na = 0.0;
                        for (int n = 0; n < J_data.length; n++) {
                            double v1 = (double) I_data[n] - I_mean;
                            double v2 = (double) J_data[n] - J_mean;
                            Na += (v1 * v2);
                            SQR1 += (v1 * v1);
                            SQR2 += (v2 * v2);
                        }
                        correlation = (float) (Na / (Math.sqrt(SQR1) * Math.sqrt(SQR2)));
                        correlation = _getSigmoid(correlation);
                    }
                }
                correlations[coord] = correlation;
            }
        }else{
            for (int j = i; j < Exp.getNumRows(); j++) {;
                int coord = j - i;
                if (i == j) {
                    correlations[coord]=1.0f;
                }else{
                    correlations[coord]=0.0f;
                }
            }
        }
        return correlations;
    }

    private float _getSigmoid(float V) {
        if((A == 0.0f) && (M == 0.0f)){
            return V; // passthrough - must be notated in manual as A=0,M=0 =/= sigmoid(V)=V
        }else{
            float s = (float) (1.0f / (1 + Math.exp(A * -1 * (Math.abs(V) - M))));
            if(V<0.0f){
                s=s*-1;
            }
            return s;
        }
        
    }

    public float[] doWork_sigmoid(String s) {
        int i = Integer.parseInt(s);
        int size = Adj.getNumRows() - i;
        float[] adjacency = new float[size];
        float[] I_data = Adj.getRowByIndex(i);
        for (int j = i; j < Adj.getNumRows(); j++) {
            int coord = j - i;
            if (i == j) {
                adjacency[coord] = 1.0f;
            } else {
                float value = I_data[j];
                adjacency[coord] = _getSigmoid(value);
            }
        }
        return adjacency;
    }

    public float[] doWork_tom(String s) {
        int i = Integer.parseInt(s);
        int size = D - i;
        float i_k = Adj.findK(i, i);
        //System.err.println(i_k);
        float Approx_Cut = 0.05f; 
/* Approx_Cut: The level below which we assume all entries for i or j (as applied)
   are zero - in order for i_k <= Approx_Cut == true, the average adjacency edge
   strength must be equal to i_k/D, where D is the cardinality of the gene set
   So, this is approximate, but will save much time.     
*/
        float[] TOM = new float[size];
        if(i_k >= Approx_Cut){ // if i_k == 0, all a_ij == 0 all TOM(a_in==0)
            for (int j = i; j < D; j++) {
                int coord = j - i;
                float tom = 0.0f;
                if (i == j) {
                    tom = 1.0f;
                } else {
                    float product = 0;
                    float j_k = Adj.findK(j, j);
                    if(j_k <= Approx_Cut) continue; // if j_k == 0, j has no neighbors, can skip all iterations and checks
                    for (int u = 0; u < D; u++) {
                        if (u == i) continue;
                        if (u == j) continue;
                        if (Adj.testValue(i, u) == false) continue;
                        if (Adj.testValue(j, u) == false) continue;
                        float this_product = Math.abs(Adj.getValueByEntry(i, u) * Adj.getValueByEntry(j, u));                        
                        product += this_product;
                    }
                    float k_min = Math.min(i_k, j_k);
                    float DFIJ = Math.abs(Adj.getValueByEntry(i, j));
                    tom = ((product + DFIJ) / (k_min + 1 - DFIJ));
                }
                if(tom>1.0f){
                    tom=1.0f;
                }
                TOM[coord] = tom;
            }
        }
        return TOM;
    }
    public float[] doWork_diff(String s) {
        int i = Integer.parseInt(s);
        int size = D - i;
        float[] Diff = new float[size];
        for (int j = i; j < D; j++) {
            int coord = j - i;
            float diff = 0.0f;
            if (i == j) {
                diff = 0.0f;
            } else {
                diff = Adj_b.getValueByEntry(i,j) - Adj.getValueByEntry(i,j);
            }
            Diff[coord] = diff;
        }
        return Diff;
    }
}
