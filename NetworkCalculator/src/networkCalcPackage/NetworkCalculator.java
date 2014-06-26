package networkCalcPackage;

import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.lang.Math;
import org.apache.commons.cli.*;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class NetworkCalculator {

	public static void main(String[] args) {
		CommandLineParser parser = new BasicParser();
		Options options = buildOptions();
		String pathIn = null;
		String SimOut = null;
		String AdjOut = null;
		String TomOut = null;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("d")){
				pathIn=cmd.getOptionValue("d");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("a")){
				AdjOut=cmd.getOptionValue("a");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("s")){
				SimOut=cmd.getOptionValue("s");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("t")){
				TomOut=cmd.getOptionValue("t");
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
		File file = null;
		try{
			file = new File(pathIn);
		}
		catch (NullPointerException e){
			System.err.println("No file found to read.\n");
			System.exit(0);
		}
		
		int[] FileDimensions = new int [2]; 
		FileDimensions = getFileDimensions(file);
		
		String[] Loci = new String[FileDimensions[0]];
		Loci = loadLoci(file,FileDimensions[0]);
		
		System.err.println("Loading Data File\n");
		
		double[][] DataFrame = new double[FileDimensions[0]][FileDimensions[1]];
		DataFrame = loadData(file,FileDimensions);
		
		/*System.err.println("Calculating Similarity\n");
		double[][] Similarity = new double[FileDimensions[0]][FileDimensions[0]]; 
		Similarity = calculateSimilarity(DataFrame,FileDimensions);
		
		System.err.println("Printing similarity to file...\n");
		printMatrixToFile(Similarity,Loci,SimOut,FileDimensions);
		
		System.err.println("Calculating Adjacency...\n");
		double[][] Adjacency = new double[FileDimensions[0]][FileDimensions[0]]; 
		Adjacency = calculateSigmoidAdjacency(Similarity,0.8,15);
		
		System.err.println("Printing Adjacency to file...\n");
		printMatrixToFile(Adjacency,Loci,AdjOut,FileDimensions);
				
		System.err.println("Masking Adjacency...\n");
		Adjacency = maskMatrix(Adjacency,0.01);
		*/
		
		System.err.println("Calculating Similarity\n");
		double[][] CurrentMatrix  = new double[FileDimensions[0]][FileDimensions[0]]; 
		CurrentMatrix = calculateSimilarity(DataFrame,FileDimensions);
		
		System.err.println("Calculating Adjacency...\n");
		CurrentMatrix = calculateSigmoidAdjacency(CurrentMatrix,0.8,15);
		
		System.err.println("Masking Adjacency...\n");
		CurrentMatrix = maskMatrix(CurrentMatrix,0.01);
		
		System.err.println("Calculating TOM...\n");
		CurrentMatrix = calculateTOM(CurrentMatrix);
		
		System.err.println("Printing TOM to file...\n");
		printMatrixToFile(CurrentMatrix,Loci,TomOut,FileDimensions);
	}
	
	private static Options buildOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option datafile = OptionBuilder.withArgName("datafile")
				.hasArg()
				.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values")
				.create("d");
		Option similarity = OptionBuilder.withArgName("similarity")
				.hasArg()
				.withDescription("File for output of similarity matrix")
				.create("s");
		Option adjacency = OptionBuilder.withArgName("adjacency")
				.hasArg()
				.withDescription("File for output of adjacency matrix")
				.create("a");
		Option tom = OptionBuilder.withArgName("tom")
				.hasArg()
				.withDescription("File for output of TOM matrix")
				.create("t");
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(similarity);
		options.addOption(adjacency);
		options.addOption(tom);
		return options;
	}
	
	private static double[][] maskMatrix (double[][] DataFrame,double maskLevel){
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		//// By keeping the dimensions separate (not assuming square matrix), enables method to mask dataframes as well as network frames
		for(int i=0;i<H;i++){
			for(int j=0;j<W;j++){
				if(DataFrame[i][j]<maskLevel){
					DataFrame[i][j]=0;
				}
			}
		}
		return DataFrame;
	}
	private static int findK (double[] Row,int j){
		int K=0;
		for(int i=0;i<Row.length;i++){
			if(i==j){
			}else{
				if(Row[i] != 0){
					K++;
				}
			}
		}
		return K;
	}
	
	private static double[][] calculateTOM (double[][] DataFrame){
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		if(H != W){
			System.err.println("TOM cannot be calculated on non-square matricies\n"); // will not... 
			System.exit(0);
		}
		double[][] ReturnFrame = new double[H][W];
		for(int i=0;i<H;i++){
			int i_k = findK(DataFrame[i],i);
			for(int j=0;j<H;j++){
				double T=0;
				if(i==j){
					T=1;
				}else{
					double product=0;
					int j_k = findK(DataFrame[j],j);
					for(int u=0;u<H;u++){
						if((u != i) && (u != j) && (DataFrame[i][u] != 0) && (DataFrame[j][u] != 0)){
							product += DataFrame[i][u] * DataFrame [j][u];
						}
					}
					int k_min = Math.min(i_k, j_k);
					double DFIJ=DataFrame[i][j];
					T=(product+DFIJ)/(k_min + 1 - DFIJ);
				}
				ReturnFrame[i][j]=T;
			}
		}
		return ReturnFrame;
	}
	
	private static double[][] calculateSigmoidAdjacency (double[][] DataFrame,double mu, double alpha){
		/* Lifted shamelessly from WGCNA:
			function (ss, mu = 0.8, alpha = 20){
				1/(1 + exp(-alpha * (ss - mu)))
			}
		 */
		int H = DataFrame.length;
		int W = DataFrame[0].length;
		double[][] Adjacency = new double[H][W];
		for(int i=0;i<H;i++){
			for(int j=i;j<H;j++){
				double adjacency=0.0;
				if(i==j){
					adjacency = 1.0;
				}else{
					adjacency = 1/(1+Math.exp(alpha*-1*(DataFrame[i][j]-mu)));
				}
				Adjacency[i][j] = adjacency;
				Adjacency[j][i] = adjacency;
			}
		}
		return Adjacency;
		// is DataFrame here automagically garbage collected by the JVM at the close of this method? Or does the block hang around?
	}
	
	private static void printMatrixToFile (double[][] DataFrame,String[] Loci, String path,int[] Dims){
		try {
			PrintWriter writer = new PrintWriter(path,"UTF-8");
			for(int i=0;i<Dims[0];i++){
				double[] Row = new double[Dims[1]];
				Row = DataFrame[i];
				String row = StringUtils.join(Row,",");
				writer.println(row+"\n");
			}
			writer.close();
		} catch (Exception e){
			// 
		}
		
	}
	/*
	private static void printMatrix (ArrayList<ArrayList<Double>> DataFrame){
		int S=DataFrame.size();
		for(int i=0;i<S;i++){
			ArrayList<Double> Row = new ArrayList<Double>();
			Row=DataFrame.get(i);
			System.out.println(ArrayUtils.toString(Row)+"\n");
		}
	}
	*/
	private static double[][] calculateSimilarity (double[][] DataFrame,int[] Dims){
		double[][] Similarity = new double[Dims[0]][Dims[0]];
		for(int i=0;i<Dims[0];i++){
			for(int j=i;j<Dims[0];j++){
				double correlation=0.0;
				if(i==j){
					correlation = 1.0;
				}else{
					double[] I_data = new double[Dims[1]];
					double[] J_data = new double[Dims[1]];
					I_data = DataFrame[i];
					J_data = DataFrame[j];
					PearsonsCorrelation corr = new PearsonsCorrelation(); 
					correlation = corr.correlation(I_data,J_data);
				}
				Similarity[i][j]=correlation;
				Similarity[j][i]=correlation;
			}
		}
		return Similarity;
	}
	
private static int[] getFileDimensions (File file) {
	int[] dimensions = new int[2];
	// pre-declaring sizes allows use of non-dynamic double[][] instead of nested ArrayLists. 
	// performance gain over ArrayList per-entry is very small, but with 7k gene#, we have 49 million entries - or at least (49 million * .5)ish
    try {
    	Scanner scanner = new Scanner(file);
		String header[] = scanner.nextLine().split("\t");
		dimensions[0]=0; // Frame height (minus header)
		dimensions[1]=header.length - 1;  // Frame width (minus rowID)
		while(scanner.hasNextLine()){
			String line=scanner.nextLine();
			String[] Line = line.split("\t");
			int this_width = Line.length - 1;
			if(this_width != dimensions[1]){
				fileDimErr();
			}
			dimensions[0]+=1;
		}
		scanner.close();
    } catch (FileNotFoundException e){
		e.printStackTrace();	    	
    } finally {
    }
    return dimensions;
}

private static void fileDimErr () {
	System.err.println("data file is not a square data file");
	System.exit(0);
}

private static String[] loadLoci (File file,int Dim) {
	String[] Loci = new String[Dim];
	try {
		Scanner scanner = new Scanner(file);
		String header[] = scanner.nextLine().split("\t");
		int it=0;
		while(scanner.hasNextLine()){
			String line=scanner.nextLine();
			String[] Line = line.split("\t");
			Loci[it] = Line[0];
			it++;
		}
		scanner.close();		
	} catch (FileNotFoundException e){
		e.printStackTrace();
	}
	return Loci;
}

private static double[][] loadData (File file, int[] Dims) {
	double[][] DataFrame = new double[Dims[0]][Dims[1]];
	try {
		Scanner scanner = new Scanner(file);
		String header[] = scanner.nextLine().split("\t");
		int it=0;
		while(scanner.hasNextLine()){
			String line=scanner.nextLine();
			String[] Line = line.split("\t");
			double[] data = new double[Dims[1]];
			for(int i=1;i<Line.length;i++){
				try {
					int I=i-1;
					double value = Double.parseDouble(Line[i]);
					data[I]=value;
				}catch(NumberFormatException e){
					e.printStackTrace();
				}
			}
			DataFrame[it]=data;
		}
		scanner.close();
	} catch (FileNotFoundException e){
		e.printStackTrace();
	}
	return DataFrame;
}

}
