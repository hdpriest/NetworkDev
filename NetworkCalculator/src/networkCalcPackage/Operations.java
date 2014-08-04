package networkCalcPackage;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class Operations {
	
	private static double GINI (double[] array1,double[] array2){
		double GINI_coeff;
		double Numerator=0.0;
		double Denominator=0.0;
		int[] sortRanks1 = getIndicesInOrder(array1);
		int[] sortRanks2 = getIndicesInOrder(array2);
		for(int i=0;i<array1.length;i++){
			double v2=((2*(i+1))-array1.length-1) * array1[sortRanks1[i]];
			double v1=((2*(i+1))-array1.length-1) * array1[sortRanks2[i]];
			Denominator+=v2;
			Numerator += v1;
		}
		GINI_coeff=Numerator/Denominator;
		//System.err.println("Obtained Num: "+ Numerator +"\nObtained Denom: " + Denominator + "\n" + "Gini: " + GINI_coeff);
		return GINI_coeff;
	}
	
	private static int[] getIndicesInOrder(double[] array) {
	    Map<Integer, Double> map = new HashMap<Integer, Double>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);

	    List<Entry<Integer, Double>> l = 
	                           new ArrayList<Entry<Integer, Double>>(map.entrySet());

	    Collections.sort(l, new Comparator<Entry<?, Double>>() {
	            @Override
	            public int compare(Entry<?, Double> e1, Entry<?, Double> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });

	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();

	    return result;
	}
	
	public static GCNMatrix calculateGINIcoefficient (GCNMatrix InputFrame){
		int D = InputFrame.getNumRows();
		GCNMatrix Similarity = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				double GCC1=0.0;
				double GCC2=0.0;
				if(i==j){
					GCC1 = 1.0;
					GCC2 = 0.0;
				}else{
					double[] I_data = InputFrame.getRowByIndex(i);
					double[] J_data = InputFrame.getRowByIndex(j);
					GCC1 = GINI(I_data,J_data);
					GCC2 = GINI(J_data,I_data);
					//System.err.println();
				}
				if(Math.abs(GCC1)<Math.abs(GCC2)){
					Similarity.setValueByEntry(GCC1,i,j);
					Similarity.setValueByEntry(GCC1,j,i);
				}else if (Math.abs(GCC2)<Math.abs(GCC1)){
					Similarity.setValueByEntry(GCC2,i,j);
					Similarity.setValueByEntry(GCC2,j,i);
				}else{
					Similarity.setValueByEntry(GCC1,i,j);
					Similarity.setValueByEntry(GCC1,j,i);
				}
			}
		}
		return Similarity;
	}
	
	public static GCNMatrix compareNetworksViaTOM (GCNMatrix Net1, GCNMatrix Net2){
		int D = Net1.getNumRows();
		GCNMatrix ReturnFrame = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			int i_k = Net1.findK(i, i);
			for(int j=0;j<D;j++){
				double T=0;
				if(i==j){
					double product=0;
					int j_k = Net2.findK(j,j);
					for(int u=0;u<D;u++){
						if((u != i) && (u != j) && (Net1.testValue(i, u)) && (Net2.testValue(j, u))){
							product += Net1.getValueByEntry(i,u) * Net2.getValueByEntry(j,u);
						}
					}
					int k_min = Math.min(i_k, j_k);
					double DFIJ=1;
					T=(product+DFIJ)/(k_min + 1 - DFIJ);
					//T=(product+DFIJ)/(k_min + 1);
				}else{
				/*	double product=0;
					int j_k = Net1.findK(j,j);
					for(int u=0;u<D;u++){
						if((u != i) && (u != j) && (Net1.testValue(i, u)) && (Net2.testValue(j, u))){
							product += InputFrame.getValueByEntry(i,u) * InputFrame.getValueByEntry(j,u);
						}
					}
					int k_min = Math.min(i_k, j_k);
					double DFIJ=InputFrame.getValueByEntry(i,j);
					T=(product+DFIJ)/(k_min + 1 - DFIJ);*/
				}
				ReturnFrame.setValueByEntry(T, i, j);
			}
		}
		return ReturnFrame;
	}
	
	public static GCNMatrix calculateTOM (GCNMatrix InputFrame){
		int D = InputFrame.getNumRows();
		GCNMatrix ReturnFrame = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			int i_k = InputFrame.findK(i, i);
			for(int j=0;j<D;j++){
				double T=0;
				if(i==j){
					T=1;
				}else{
					double product=0;
					int j_k = InputFrame.findK(j,j);
					for(int u=0;u<D;u++){
						if((u != i) && (u != j) && (InputFrame.testValue(i, u)) && (InputFrame.testValue(j, u))){
							product += InputFrame.getValueByEntry(i,u) * InputFrame.getValueByEntry(j,u);
						}
					}
					int k_min = Math.min(i_k, j_k);
					double DFIJ=InputFrame.getValueByEntry(i,j);
					T=(product+DFIJ)/(k_min + 1 - DFIJ);
				}
				ReturnFrame.setValueByEntry(T, i, j);
			}
		}
		return ReturnFrame;
	}
	
	public static GCNMatrix calculateSigmoidAdjacency (GCNMatrix Similarity, double mu, double alpha){
		int D = Similarity.getNumRows();
		GCNMatrix Adjacency = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				double adjacency=0.0;
				if(i==j){
					adjacency = 1.0;
				}else{
					adjacency = 1/(1+Math.exp(alpha*-1*(Math.abs(Similarity.getValueByEntry(i,j))-mu)));
				}
				Adjacency.setValueByEntry(adjacency, i, j);
				Adjacency.setValueByEntry(adjacency, j, i);
			}
		}
		return Adjacency;
	}
	
	public static GCNMatrix calculateDifference (GCNMatrix mat1, GCNMatrix mat2){
		int D = mat1.getNumRows();
		GCNMatrix Difference = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				double v1=mat1.getValueByEntry(i, j);
				double v2=mat2.getValueByEntry(i, j);
				//if((v1 != 0) & (v2 != 0)){
					double d1 =v1-v2;
					double d2 =v1-v2;
				//System.out.println("Val1: " + v1 +" Val2: " + v2 + " diff " + d);
					Difference.setValueByEntry(d1,i,j);
					Difference.setValueByEntry(d2,j,i);
				//}
			}
		}
		return Difference;
	}
	
	public static GCNMatrix calculateSimilarity (GCNMatrix Expression){
		int D = Expression.getNumRows();
		GCNMatrix Similarity = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				double correlation=0.0;
				if(i==j){
					correlation = 1.0;
				}else{
					double[] I_data = Expression.getRowByIndex(i);
					double[] J_data = Expression.getRowByIndex(j);
					PearsonsCorrelation corr = new PearsonsCorrelation();
					correlation = corr.correlation(I_data,J_data);
				}
				Similarity.setValueByEntry(correlation,i,j);
				Similarity.setValueByEntry(correlation,j,i);
			}
		}
		System.exit(0);
		return Similarity;
	}
	
	public static void generateHistogram (GCNMatrix DataFrame, String pathOut, String Title,String Xlab, String Ylab) {
		int H = DataFrame.getNumRows();
		int W = DataFrame.getNumColumns();
		double[] Histogram=new double[201];
		DecimalFormat df = new DecimalFormat("#.####");
		df.setRoundingMode(RoundingMode.HALF_UP);
		for(int i=0;i<H;i++){
			for(int j=0;j<W;j++){
				//System.out.println("Val: "+DataFrame[i][j]+"\n");
				if(DataFrame.getValueByEntry(i,j) != 0){
					try{
						Double v= ((Double.valueOf(df.format(DataFrame.getValueByEntry(i,j))))+1)*100;
						int value = v.intValue();
						Histogram[value]++;
					}catch(NumberFormatException ex){
						System.out.println("Obtain " + DataFrame.getValueByEntry(i,j) +" from matrix.");
						System.exit(1);
					}
				}
			}
		}
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries series = new XYSeries("Values");
		for(int x=0;x<201;x++){
			double a = ((double)x/100)-1.0;
			Double A = Double.valueOf(df.format(a));
			System.out.println("Index: " + x + " Value: "+ A +"\tObs: "+Histogram[x]);
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
	}
}
