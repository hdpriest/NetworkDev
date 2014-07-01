package networkCalcPackage;

import java.io.PrintWriter;
import java.lang.Math;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class GCNMatrix {
	
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
				if(DataFrame[R][i] != 0){
					K++;
				}
			}
		}
		return K;
	}
	
	public void addRow (double[] Row){
		int I=X_iterator;
		for(int i=0;i<Row.length;i++){
			DataFrame[I][i]=Row[i];
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
			for(int i=0;i<DataFrame.length;i++){
				double[] Row = new double[DataFrame[i].length];
				Row = DataFrame[i];
				String row = StringUtils.join(Row,Sep);
				writer.println(row+"\n");
			}
			writer.close();
		} catch (Exception e){
			// 
		}
		
	}
}
