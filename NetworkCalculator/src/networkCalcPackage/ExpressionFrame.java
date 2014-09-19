package networkCalcPackage;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Scanner;
import java.math.RoundingMode;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotOrientation;


class ExpressionFrame {
	private int N;
	private int M;
	private float[][] DataFrame;
	private float[] k;
	private float[] means;
	private float[] gccSums;
	private String[] X_lab;
	private String[] Y_lab;
	private int X_iterator;
	    
	public ExpressionFrame (int Dim1, int Dim2) {
		DataFrame = new float[Dim1][Dim2];
		N = Dim1;
		M = Dim2;
		k = new float[Dim1];
		means = new float[Dim1];
		gccSums= new float[Dim1];
		X_lab = new String[Dim1];
		Y_lab = new String[Dim2];
		X_iterator = -1;
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
		float[] Row = _getRow(I);
		double[] dRow = new double[M];
		for(int i=0;i<M;i++){
			dRow[i] = (double) Row[i];
		}
		return dRow;
	}
        
	private float[] _getRow(int I){
		float[] Row = new float[M];
		for(int j=0;j<M;j++){
			Row[j]=_getValueByEntry(I,j);
		}
		return Row;
	}
	
	public float[] getRowByIndex (int I){
		if(DataFrame[I] != null){
			return _getRow(I);
		}else{
			System.err.println("Cannot get row "+I+" from matrix.\n\n");
			System.exit(0);
		}
		return null;
	}
	
	private float _getValueByEntry (int I,int J){
		return DataFrame[I][J];
	}
	
	public float getValueByEntry (int I,int J){
		return _getValueByEntry(I,J);
	}
	
	public void setValueByEntry (float Value,int I, int J){
		_setValueByEntry(Value,I,J);
	}
	
	private void _setValueByEntry (float Value,int I, int J){
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
                                float v = _getValueByEntry(i,j);
				if(v<maskLevel){
                                        _setValueByEntry(0.0f,i,j);
				}
			}
		}
	}
	
	public boolean testValue (int i,int j){
		boolean res=true;
		if(_getValueByEntry(i,j) == 0){
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
					thisK+=_getValueByEntry(i,j);
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

	public void calculateMeans() {
		for(int i=0;i<DataFrame.length;i++){
			float[] Row = _getRow(i);
			float s=0;
			for(int j=0;j<Row.length;j++){
				s+=Row[j];
			}
			means[i]=s/Row.length;
		}
		
	}

	public float getMean(int i){
		return means[i];
	}
	
	public float getGiniDenom(int i) {
		return gccSums[i];
	}
	
	public void calculateGiniSums() {
		for(int i=0;i<DataFrame.length;i++){
			float[] array1 = _getRow(i);
			float Denominator=0.0f;
			int[] sortRanks1 = Operations.getIndicesInOrder(array1);
			for(int j=0;j<array1.length;j++){
				float v2=((2*(j+1))-array1.length-1) * array1[sortRanks1[j]];
				Denominator+=v2;
			}
			gccSums[i]=Denominator;
		}
	}
}
