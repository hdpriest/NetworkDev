package networkCalcPackage;

import java.lang.Math;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class GCN_Matrix {
	
	private double[][] DataFrame;
	private String[] X_lab;
	private String[] Y_lab;
	private int X_iterator;
	
	public GCN_Matrix (int[] Dims) {
		DataFrame = new double[Dims[0]][Dims[1]];
		X_lab = new String[Dims[0]];
		Y_lab = new String[Dims[1]];
		X_iterator = -1;
	}
	
	public void resetIterator () {
		X_iterator=-1;
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
	
	public double[] getNextRow () {
		double[] thisRow = new double[DataFrame[0].length];
		thisRow=DataFrame[X_iterator+1];
		X_iterator++;
		return thisRow;
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

}
