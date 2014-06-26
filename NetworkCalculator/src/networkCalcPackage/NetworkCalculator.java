package networkCalcPackage;

import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.ArrayList;

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
		
		//ArrayList<ArrayList<Double>> DataFrame = loadData(file);
		
		System.err.println("Calculating Similarity\n");
		double[][] Similarity = new double[FileDimensions[0]][FileDimensions[1]]; 
		Similarity = calculateSimilarity(DataFrame,FileDimensions);
		
		System.err.println("Printing similarity to file...\n");
		printMatrixToFile(Similarity,Loci,SimOut,FileDimensions);
		
		System.err.println("Calculating Adjacency...\n");
		double[][] Adjacency = new double[FileDimensions[0]][FileDimensions[1]]; 
		Adjacency = calculateSigmoidAdjacency(Similarity,0.8,15);
		
		/*
		System.err.println("Printing Adjacency to file...\n");
		printMatrixToFile(Adjacency,Loci,AdjOut);
		*/		
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
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(similarity);
		options.addOption(adjacency);
		return options;
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
	}
	
	private static void printMatrixToFile (double[][] DataFrame,String[] Loci, String path,int[] Dims){
		try {
			PrintWriter writer = new PrintWriter(path,"UTF-8");
			for(int i=0;i<Dims[0];i++){
				double[] Row = new double[Dims[1]];
				Row = DataFrame[i];
				String row = StringUtils.join(Row,",");
				writer.println(row+"\n");
				//System.out.println(ArrayUtils.toString(Row)+"\n");
			}
			writer.close();
		} catch (Exception e){
			
		}
		
	}
	
	private static void printMatrix (ArrayList<ArrayList<Double>> DataFrame){
		int S=DataFrame.size();
		for(int i=0;i<S;i++){
			ArrayList<Double> Row = new ArrayList<Double>();
			Row=DataFrame.get(i);
			System.out.println(ArrayUtils.toString(Row)+"\n");
		}
	}
	
	private static double[][] calculateSimilarity (double[][] DataFrame,int[] Dims){
		double[][] Similarity = new double[Dims[0]][Dims[0]];
		for(int i=0;i<Dims[0];i++){
			//ArrayList<Double> Row = new ArrayList<Double>();
			for(int j=i;j<Dims[0];j++){
				double correlation=0.0;
				if(i==j){
					correlation = 1.0;
				}else{/*
					ArrayList<Double> I = DataFrame.get(i);
					ArrayList<Double> J = DataFrame.get(j);
					//Double[] I_data = new Double[I.size()];
					//Double[] J_data = new Double[J.size()];
					Double[] Id = new Double[I.size()];
					Double[] Jd = new Double[J.size()];
					Id = I.toArray(Id);
					Jd = J.toArray(Jd);
					double[] I_data = ArrayUtils.toPrimitive(Id);
					double[] J_data = ArrayUtils.toPrimitive(Jd);*/
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
	/*
	private static ArrayList<String> loadLoci (File file) {
		ArrayList<String> Loci = new ArrayList<String>();
		try {
			Scanner scanner = new Scanner(file);
			String header[] = scanner.nextLine().split("\t");
			while(scanner.hasNextLine()){
				String line=scanner.nextLine();
				String[] Line = line.split("\t");
				Loci.add(Line[0]);
			}
			scanner.close();		
		} catch (FileNotFoundException e){
			e.printStackTrace();
		}
		return Loci;
	}
	
	private static ArrayList<ArrayList<Double>> loadData (File file) {
		ArrayList<ArrayList<Double>> DataFrame = new ArrayList<ArrayList<Double>>();
		try {
			Scanner scanner = new Scanner(file);
			String header[] = scanner.nextLine().split("\t");
			while(scanner.hasNextLine()){
				String line=scanner.nextLine();
				String[] Line = line.split("\t");
				ArrayList<Double> Data = new ArrayList<Double>();
				for(int i=1;i<Line.length;i++){
					try {
						double value = Double.parseDouble(Line[i]);
						Data.add(value);
					}catch(NumberFormatException e){
						e.printStackTrace();
					}
				}
				DataFrame.add(Data);
			}
			scanner.close();
		} catch (FileNotFoundException e){
			e.printStackTrace();
		}
		return DataFrame;
	}

}
*/
private static int[] getFileDimensions (File file) {
	int[] dimensions = new int[2];
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

//private static ArrayList<ArrayList<Double>> loadData (File file) {
private static double[][] loadData (File file, int[] Dims) {
	double[][] DataFrame = new double[Dims[0]][Dims[1]];
//	ArrayList<ArrayList<Double>> DataFrame = new ArrayList<ArrayList<Double>>();
	try {
		Scanner scanner = new Scanner(file);
		String header[] = scanner.nextLine().split("\t");
		int it=0;
		while(scanner.hasNextLine()){
			String line=scanner.nextLine();
			String[] Line = line.split("\t");
			double[] data = new double[Dims[1]];
			//ArrayList<Double> Data = new ArrayList<Double>();
			for(int i=1;i<Line.length;i++){
				try {
					int I=i-1;
					double value = Double.parseDouble(Line[i]);
					data[I]=value;
					//Data.add(value);
				}catch(NumberFormatException e){
					e.printStackTrace();
				}
			}
			DataFrame[it]=data;
			//DataFrame.add(Data);
		}
		scanner.close();
	} catch (FileNotFoundException e){
		e.printStackTrace();
	}
	return DataFrame;
}

}
