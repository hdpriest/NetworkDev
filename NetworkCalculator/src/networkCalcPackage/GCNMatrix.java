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
import org.jfree.chart.plot.PlotOrientation;

class GCNMatrix {
	
	private double[][] DataFrame;
	private String[] X_lab;
	private String[] Y_lab;
	private int X_iterator;
	
	public GCNMatrix (int Dim1, int Dim2) {
		DataFrame = new double[Dim1][Dim2];
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
			for(int j=0;j<W;j++){
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
	
	public double[] getRowByIndex (int I){
		double[] Row = new double[DataFrame[0].length];
		if(DataFrame[I] != null){
			Row = DataFrame[I];
			return Row;
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return Row;
	}
	
	public double getValueByEntry (int I,int J){
		return DataFrame[I][J];
	}
	
	public void setValueByEntry (double Value,int I, int J){
		DataFrame[I][J]=Value;
	}
	
	public void setRowNames (String[] Rows) {
		X_lab = Rows;
	}
	
	public void setColumnNames (String[] Cols){
		Y_lab = Cols;
	}
	
	public double[] getNextRow () {
		double[] thisRow = new double[DataFrame[0].length];
		thisRow=DataFrame[X_iterator+1];
		X_iterator++;
		return thisRow;
	}
	
	public void maskMatrix (double maskLevel) {
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
	
	public int findK (int R,int j){
		int K=0;
		for(int i=0;i<DataFrame[R].length;i++){
			if(i==j){
			}else{
				if(DataFrame[R][i] > 0){
					K++;
				}
			}
		}
		return K;
	}
	
	public void addRow (double[] Row){
		int I=X_iterator;
		for(int i=0;i<Row.length;i++){
			DataFrame[I+1][i]=Row[i];
		}
		X_iterator++;
	}
	
	public void changeRow (int I,double[] Row){
		for(int i=0;i<Row.length;i++){
			DataFrame[I][i]=Row[i];
		}
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
				double[] Row = new double[DataFrame[i].length];
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
