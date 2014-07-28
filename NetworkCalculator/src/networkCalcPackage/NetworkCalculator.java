package networkCalcPackage;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.lang.Math;

import org.apache.commons.cli.*;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class NetworkCalculator {

	private static void makeSimilarity(String[] args){
		CommandLineParser parser = new BasicParser();
		Options options = buildSimilarityOptions();
		String pathIn=null;
		String Out=null;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("d")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("o")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			pathIn=cmd.getOptionValue("d");
			Out=cmd.getOptionValue("o");
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
		
		// *** TODO : MAKE SIMILARITY FUNC
		File file = null;
		try{
			file = new File(pathIn);
		}
		catch (NullPointerException e){
			System.err.println("No file found to read.\n");
			System.exit(0);
		}
		String sep ="\t";
		int[] FileDimensions = new int [2]; 
		FileDimensions = getFileDimensions(file, sep);
		
		System.err.println("Loading Data File\n");
		
		GCNMatrix DataFrame = new GCNMatrix(FileDimensions[0],FileDimensions[1]);

		DataFrame = loadData(file,FileDimensions,sep);
		
		System.err.println("Calculating Similarity\n");
		GCNMatrix CurrentMatrix  = new GCNMatrix(FileDimensions[0],FileDimensions[0]); 
		CurrentMatrix = calculateSimilarity(DataFrame);
		CurrentMatrix.printMatrixToFile(Out,",");
		System.exit(0);
		
	}
	
	private static void makeNetwork(String[] args){
		CommandLineParser parser = new BasicParser();
		Options options = buildConstructOptions();
		String pathIn=null;
		String Out=null;
		double alpha=0.0;
		double mu=0.0;
		double Mask=0.0;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("d")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("o")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("a")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("m")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("M")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			pathIn=cmd.getOptionValue("d");
			Out=cmd.getOptionValue("o");
			alpha=Double.parseDouble(cmd.getOptionValue("a"));
			mu=Double.parseDouble(cmd.getOptionValue("m"));
			Mask = Double.parseDouble(cmd.getOptionValue("M"));
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
		
		// *** TODO: make network construction method
		File file = null;
		String sep = "\t";
		try{
			file = new File(pathIn);
		}
		catch (NullPointerException e){
			System.err.println("No file found to read.\n");
			System.exit(0);
		}
		
		File dir = new File(Out);
		try{
			dir.mkdir();
		} catch(SecurityException se){
			System.out.println("Cannot create directory!");
			System.exit(0);
		}
		
		int[] FileDimensions = new int [2]; 
		FileDimensions = getFileDimensions(file,sep);
		
		System.err.println("Loading Data File\n");
		
		GCNMatrix DataFrame = new GCNMatrix(FileDimensions[0],FileDimensions[1]);

		DataFrame = loadData(file,FileDimensions,sep);
		
		System.err.println("Calculating Similarity\n");
		GCNMatrix CurrentMatrix  = new GCNMatrix(FileDimensions[0],FileDimensions[0]); 
		CurrentMatrix = calculateSimilarity(DataFrame);
		String ThisOut = Out + "/Similarity.dist.jpeg";
		CurrentMatrix.generateHistogram(ThisOut, "Similarity Distribution", "Pearsons Correlation", "# Relationships");
		
		System.err.println("Calculating Adjacency...\n");
		CurrentMatrix = calculateSigmoidAdjacency(CurrentMatrix,mu,alpha);
		ThisOut = Out + "/Adjacency.dist.jpeg";
		CurrentMatrix.generateHistogram(ThisOut, "Adjacency Distribution", "Sigmoid Adjacency Value", "# Relationships");
		
//		System.err.println("Masking Adjacency...\n");
		CurrentMatrix.maskMatrix(Mask);
		
		System.err.println("Calculating TOM...\n");
		CurrentMatrix = calculateTOM(CurrentMatrix);
		
	//	System.err.println("Printing TOM to file...\n");
	//	System.err.println("Printing Matrix...\n");
	//	CurrentMatrix.printMatrixToFile(Out,",");
		CurrentMatrix.maskMatrix(Mask);
		ThisOut = Out + "/TOM.dist.jpeg";
		CurrentMatrix.generateHistogram(ThisOut,"Masked Distribution of Topological Overlaps","Topological Overlap","# of relationships");
		System.exit(0);
	}
	
	private static void compareNetworks(String[] args){
		CommandLineParser parser = new BasicParser();
		Options options = buildCompareOptions();
		/// **** TODO : need to even conceive of whats going on here
		String pathIn=null;
		String Out=null;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("d")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("o")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			pathIn=cmd.getOptionValue("d");
			Out=cmd.getOptionValue("o");
			/*
			DefaultCategoryDataset dataset = new DefaultCategoryDataset();
			dataset.setValue(6, "Profit", "Jane");
			dataset.setValue(7, "Profit", "Tom");
			dataset.setValue(8, "Profit", "Jill");
			dataset.setValue(5, "Profit", "John");
			dataset.setValue(12, "Profit", "Fred");
			JFreeChart chart = ChartFactory.createBarChart("Comparison between Salesman", 
					 "Salesman", "Profit", dataset, PlotOrientation.VERTICAL, 
					 false, true, false);
					 
			try {
				ChartUtilities.saveChartAsJPEG(new File(Out), chart, 500, 300);
			} catch (IOException e) {
				System.err.println("Problem occurred creating chart.");
			}*/
			File file = null;
			String sep = ",";
			try{
				file = new File(pathIn);
			}
			catch (NullPointerException e){
				System.err.println("No file found to read.\n");
				System.exit(0);
			}
			
			int[] FileDimensions = new int [2]; 
			FileDimensions = getFileDimensions(file, sep);
			
			System.err.println("Loading Data File\n");
			
			GCNMatrix DataFrame = new GCNMatrix(FileDimensions[0],FileDimensions[1]);

			DataFrame = loadData(file,FileDimensions,sep);
		//	System.out.println("Columns loaded... " + DataFrame.getNumColumns() +"\n");
		//	System.out.println("Rows loaded... " + DataFrame.getNumRows() +"\n");
			DataFrame.generateHistogram(Out,"Title","X label","Y label");
			System.exit(0);
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
		
		// *** TODO: make network construction method
	}
	
	private static void viewNetwork(String[] args){
		CommandLineParser parser = new BasicParser();
		Options options = buildViewOptions();
		// **** TODO -- again, what is this for??
		String pathIn=null;
		String Out=null;
		double alpha=0.0;
		double mu=0.0;
		try{
			CommandLine cmd = parser.parse(options, args);
			HelpFormatter formatter = new HelpFormatter();
			if(cmd.hasOption("h")){
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("n")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("o")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("a")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			if(cmd.hasOption("m")){
			}else{
				formatter.printHelp( "java -jar jarfile.jar", options );
				System.exit(0);
			}
			pathIn=cmd.getOptionValue("n");
			Out=cmd.getOptionValue("o");
			alpha=Double.parseDouble(cmd.getOptionValue("a"));
			mu=Double.parseDouble(cmd.getOptionValue("m"));
		}
		catch(ParseException exp){
			System.err.println("Problem parsing arguments:\n" + exp.getMessage());
			System.err.println("Exiting...\n");
			System.exit(0);
		}
	}
	
	private static void baseOptions(String[] args){
		Options options = buildOptions();
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp( "java -jar jarfile.jar", options );
		System.exit(0);
	}
	
	public static void main(String[] args) {
		if(args.length==0){
			baseOptions(args);
		}else{
			switch (args[0]){
				case "similarity": makeSimilarity(args);
				break;
				case "construct": makeNetwork(args);
				break;
				case "compare": compareNetworks(args);
				break;
				case "view": viewNetwork(args);
				break;
				default: baseOptions(args);
				break;
			}
		}
		/*
		Options options = buildOptions();
		
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
		
		System.err.println("Loading Data File\n");
		
		GCNMatrix DataFrame = new GCNMatrix(FileDimensions[0],FileDimensions[1]);

		DataFrame = loadData(file,FileDimensions);
		
		System.err.println("Calculating Similarity\n");
		GCNMatrix CurrentMatrix  = new GCNMatrix(FileDimensions[0],FileDimensions[0]); 
		CurrentMatrix = calculateSimilarity(DataFrame);
		
		System.err.println("Calculating Adjacency...\n");
		CurrentMatrix = calculateSigmoidAdjacency(CurrentMatrix,0.8,15);
		
		System.err.println("Masking Adjacency...\n");
		CurrentMatrix.maskMatrix(0.01);
		
		System.err.println("Calculating TOM...\n");
		CurrentMatrix = calculateTOM(CurrentMatrix);
		
		System.err.println("Printing TOM to file...\n");
		CurrentMatrix.printMatrixToFile(TomOut,",");
		*/
	}
	
	private static Options buildOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option similarity = OptionBuilder.withArgName("similarity")
				.hasArg()
				.withDescription("construct a similarity matrix and print it to a file")
				.create("similarity");
		Option construct = OptionBuilder.withArgName("construct")
				.hasArg()
				.withDescription("construct a network via application of a topological overlap calculation")
				.create("construct");
		Option compare = OptionBuilder.withArgName("compare")
				.hasArg()
				.withDescription("compare several network matricies")
				.create("compare");
		Option view = OptionBuilder.withArgName("view")
				.hasArg()
				.withDescription("query a single network for node properties")
				.create("view");
		options.addOption(help);
		options.addOption(construct);
		options.addOption(similarity);
		options.addOption(compare);
		options.addOption(view);
		return options;
	}
	private static Options buildCompareOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option datafile = OptionBuilder.withArgName("datafile")
				.hasArg()
				.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values")
				.create("d");
		Option output = OptionBuilder.withArgName("similarity")
				.hasArg()
				.withDescription("File for output of output")
				.create("o");
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(output);
		//System.out.println("This method is not yet implemented\n");
		//System.exit(0);
		return options;
	}
	private static Options buildViewOptions (){
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
		System.out.println("This method is not yet implemented\n");
		System.exit(0);
		return options;
	}
	private static Options buildSimilarityOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option datafile = OptionBuilder.withArgName("datafile")
				.hasArg()
				.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values")
				.create("d");
		Option output = OptionBuilder.withArgName("output")
				.hasArg()
				.withDescription("File for output of similarity matrix")
				.create("o");
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(output);
		return options;
	}
	private static Options buildConstructOptions (){
		Options options = new Options();
		Option help = new Option( "h", "print this message" );
		Option datafile = OptionBuilder.withArgName("datafile")
				.hasArg()
				.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values")
				.create("d");
		Option alpha = OptionBuilder.withArgName("alpha")
				.hasArg()
				.withDescription("alpha parameter for sigmoid adjacency calculation (i.e. 20)")
				.create("a");
		Option mu = OptionBuilder.withArgName("mu")
				.hasArg()
				.withDescription("mu parameter for sigmoid adjacency calculation (i.e. 0.8)")
				.create("m");
		Option Mask = OptionBuilder.withArgName("Mask")
				.hasArg()
				.withDescription("parameter for masking matricies - values below this will be set to zero")
				.create("M");
		Option output = OptionBuilder.withArgName("output")
				.hasArg()
				.withDescription("Temporary Directory")
				.create("o");
		options.addOption(help);
		options.addOption(datafile);
		options.addOption(alpha);
		options.addOption(mu);
		options.addOption(Mask);
		options.addOption(output);
		return options;
	}

	private static GCNMatrix calculateTOM (GCNMatrix InputFrame){
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

	private static GCNMatrix calculateSigmoidAdjacency (GCNMatrix Similarity, double mu, double alpha){
		int D = Similarity.getNumRows();
		GCNMatrix Adjacency = new GCNMatrix(D,D);
		for(int i=0;i<D;i++){
			for(int j=i;j<D;j++){
				double adjacency=0.0;
				if(i==j){
					adjacency = 1.0;
				}else{
					adjacency = 1/(1+Math.exp(alpha*-1*(Similarity.getValueByEntry(i,j)-mu)));
				}
				Adjacency.setValueByEntry(adjacency, i, j);
				Adjacency.setValueByEntry(adjacency, j, i);
			}
		}
		return Adjacency;
	}
	
	private static GCNMatrix calculateSimilarity (GCNMatrix Expression){
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
		return Similarity;
	}
	
	private static int[] getFileDimensions (File file, String sep) {
	int[] dimensions = new int[2];
	// pre-declaring sizes allows use of non-dynamic double[][] instead of nested ArrayLists. 
	// performance gain over ArrayList per-entry is very small, but with 7k gene#, we have 49 million entries - or at least (49 million * .5)ish
    try {
    	Scanner scanner = new Scanner(file);
		String header[] = scanner.nextLine().split(sep);
		dimensions[0]=0; // Frame height (minus header)
		dimensions[1]=header.length - 1;  // Frame width (minus rowID)
		while(scanner.hasNextLine()){
			String line=scanner.nextLine();
			String[] Line = line.split(sep);
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
	System.err.println("data file is not a rectilinear data file");
	System.exit(0);
}

	private static GCNMatrix loadData (File file, int[] Dims,String sep) {
		GCNMatrix Expression = new GCNMatrix(Dims[0],Dims[1]);
		//System.err.println("Loading file of " + Dims[0] + " by " + Dims[1] + "\n");
		try {
			Scanner scanner = new Scanner(file);
			String[] header = scanner.nextLine().split(sep);
			String[] loci = new String[Dims[0]];
			Expression.setColumnNames(header);
			int it=0;
			while(scanner.hasNextLine()){
				String line=scanner.nextLine();
				String[] Line = line.split(sep);
				loci[it]=Line[0];
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
				it++;
				Expression.addRow(data);
			}
			Expression.setRowNames(loci);
			scanner.close();
		} catch (FileNotFoundException e){
			e.printStackTrace();
		}
		return Expression;
	}

}
