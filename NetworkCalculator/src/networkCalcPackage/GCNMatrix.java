package networkCalcPackage;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.math.RoundingMode;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;


class GCNMatrix {
	
	private float[][] DataFrame;
	private float[] k;
	private String[] X_lab;
	private String[] Y_lab;
	private int X_iterator;
	
        /* 
        We need to do some matrix-fu to make this memory-slimmer. 
        represent our matrix as a 1-d array of upper triangle matrix
        For a 1000x1000 matrix, number of elements that precedes row i is:
        1000+999+998 ... + (1000*(i-1)
        
        NumElementsPreceding ROW i, where N = scalar(DIM(matrix)), I is the 0-index row index:
        i( N-((i-1)/2) )
        
        And the number of elements that precede [i,j] in an upper-diagonal matrix is: 
        i( N-((i-1)/2) ) + 1) + (j-1) WHERE j>=i
        
        N=100, i = 0, j = 1 (first off-diagonal value)
        100(1) - 0 + 0 
        
        */
        
	public GCNMatrix (int Dim1, int Dim2) {
		DataFrame = new float[Dim1][Dim2];
		k = new float[Dim1];
		X_lab = new String[Dim1];
		Y_lab = new String[Dim2];
		X_iterator = -1;
	}
	
	public double[] generateHistogram (String pathOut, String Title,String Xlab, String Ylab) {
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		double[] Histogram=new double[201];
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.HALF_UP);
		for(int i=0;i<H;i++){
			for(int j=i;j<W;j++){
				//System.out.println("Val: "+DataFrame[i][j]+"\n");
				if(DataFrame[i][j] != 0.0){
					Double v= ((Double.valueOf(df.format(DataFrame[i][j])))+1)*100;
					int value = v.intValue();
					Histogram[value]++;
				}
			}
		}
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Values");
		for(int x=0;x<201;x++){
			double a = ((double)x/100)-1.0;
			Double A = Double.valueOf(df.format(a));
			//System.out.println("Index: " + x + " Value: "+ A +"\tObs: "+Histogram[x]);
			series.add((double) A, (double) Histogram[x],true);			
		}
		dataset.addSeries(series);
		JFreeChart chart = ChartFactory.createXYLineChart(Title,Xlab,Ylab,dataset, PlotOrientation.VERTICAL, 
				 false, true, false);
		final Plot plot = chart.getPlot();
		plot.setBackgroundPaint(Color.white);
		plot.setOutlinePaint(Color.black);
		chart.setBackgroundPaint(Color.white);
		chart.setAntiAlias(true);
		try {
			ChartUtilities.saveChartAsJPEG(new File(pathOut), chart, 500, 300);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
		return Histogram;
	}
	
	public void resetIterator () {
		X_iterator=-1;
	}
	
	public int getNumColumns () {
		int I = DataFrame[0].length;
		return I;
	}
	
	public int getNumRows () {
		int I = DataFrame.length;
		return I;
	}
	
	public boolean hasNext () {
		if(DataFrame[X_iterator+1] != null){
			return true;
		}else{
			return false;			
		}
	}
	public double[] getRowByIndexDbl (int I){
		float[] Row = new float[DataFrame[0].length];
                double[] dRow = new double[DataFrame[0].length];
		if(DataFrame[I] != null){
			Row =  DataFrame[I];
                        for(int i=0;i<DataFrame[0].length;i++){
                            dRow[i] = (double) Row[i];
                        }
			return dRow;
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return dRow;
	}
        
	public float[] getRowByIndex (int I){
		float[] Row = new float[DataFrame[0].length];
		if(DataFrame[I] != null){
			Row = DataFrame[I];
			return Row;
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return Row;
	}
	
	public float getValueByEntry (int I,int J){
		return DataFrame[I][J];
	}
	
	public void setValueByEntry (float Value,int I, int J){
		DataFrame[I][J]=Value;
	}
	public String[] getRowNames (){
            return X_lab;
        }
	public void setRowNames (String[] Rows) {
            X_lab = Rows;
	}
	
	public void setColumnNames (String[] Cols){
		Y_lab = Cols;
	}
	public float[][] getDataFrame () {
            return DataFrame;
        }
        
	public float[] getNextRow () {
		float[] thisRow = new float[DataFrame[0].length];
		thisRow=DataFrame[X_iterator+1];
		X_iterator++;
		return thisRow;
	}
	
	public void maskMatrix (float maskLevel) {
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		for(int i=0;i<H;i++){
			for(int j=0;j<W;j++){
				if(DataFrame[i][j]<maskLevel){
					DataFrame[i][j]=0;
				}
			}
		}
	}
	
	public boolean testValue (int i,int j){
		boolean res=true;
		if(DataFrame[i][j] == 0){
			res=false;
		}else{
			res=true;
		}
		return res;
	}
	
	public void calculateKs (){
		int H = DataFrame.length;
		for(int i=0;i<H;i++){
			float thisK=0f;
			for(int j=0;j<DataFrame[i].length;j++){
				if(i==j){
					
				}else{
					thisK+=DataFrame[i][j];
				}
			}
			k[i]=thisK;
		}	
	}
	
	public float findK (int R,int j){
		/*for(int i=0;i<DataFrame[R].length;i++){
			if(i==j){
				/// do not add DataFrame[R][j] to the k of R
			}else{
				K+=DataFrame[R][i];
			}
		}
		return K;*/
		float K=k[R];
		return K;
	}
	
	public void addRow (float[] Row){
		int I=X_iterator;
                System.arraycopy(Row, 0, DataFrame[I+1], 0, Row.length);
		X_iterator++;
	}
	
	public void changeRow (int I,float[] Row){
                System.arraycopy(Row, 0, DataFrame[I], 0, Row.length);
	}
	
	public void printMatrixToFile (String path,String Sep){
		try {
			PrintWriter writer = new PrintWriter(path,"UTF-8");
			writer.print(Sep);
			for(int y=0;y<Y_lab.length;y++){
				writer.print(Y_lab[y]);
				if(y != Y_lab.length-1){
					writer.print(Sep);
				}
			}
			writer.print("\n");
			for(int i=0;i<DataFrame.length;i++){
				float[] Row = new float[DataFrame[i].length];
				Row = DataFrame[i];
				writer.print(X_lab[i]);
				writer.print(Sep);
				for(int j=0;j<Row.length;j++){
					writer.print(Row[j]);
					if(j != Row.length-1){
						writer.print(Sep);
					}
				}
				writer.print("\n");
			}
			writer.close();
		} catch (Exception e){
			// 
		}
		
	}
}
