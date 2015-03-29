package networkCalcPackage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;

import org.apache.commons.cli.*;

public class NetworkCalculator {

    private static void speedtest(String[] args) {
        CommandLineParser parser = new BasicParser();
        Options options = buildConstructOptions();
        String pathIn = null;
        String Out = null;
        float alpha = 0.0f;
        float mu = 0.0f;
        double Mask = 0.0d;
        int threads = 0;
        String corr = null;
        try {
            CommandLine cmd = parser.parse(options, args);
            HelpFormatter formatter = new HelpFormatter();
            if (cmd.hasOption("h")) {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("d")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("c")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("o")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("a")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("m")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("M")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            if (cmd.hasOption("t")) {
            } else {
                formatter.printHelp("java -jar jarfile.jar", options);
                System.exit(0);
            }
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            pathIn = cmd.getOptionValue("d");
            Out = cmd.getOptionValue("o");
            corr = cmd.getOptionValue("c");
            alpha = Float.parseFloat(cmd.getOptionValue("a"));
            mu = Float.parseFloat(cmd.getOptionValue("m"));
            Mask = Double.parseDouble(cmd.getOptionValue("M"));
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }
        String sep = "\t";

        int[] FileDimensions = new int[2];
        FileDimensions = _getFileDimensions(pathIn, sep);

        System.err.println("Loading Data File\n");

        ExpressionFrame DataFrame = _loadData(pathIn, FileDimensions, sep);

        GCNMatrix CurrentMatrix;

        System.err.println("Calculating Similarity\n");
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", mu, alpha, threads,true);
        File theDir = new File(Out);
        if (!theDir.exists()) {
            System.out.println("creating directory: " + Out);
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                //TODO handle it
            }
        }
        System.err.println("Calculating Adjacency...\n");
        CurrentMatrix.calculateKs();
        System.err.println("Calculating TOM...\n");
        CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);
        System.err.println("done");
        System.exit(0);
    }

    private static void _checkCompareOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d1")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d2")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("dnA")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("dpA")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("o")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }

    private static void _checkPermuteOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d1")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d2")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("o")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("p")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }

    private static void _checkScalefreeOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d1")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d2")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("o")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }    
    
    private static void _checkConstructOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);

            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("o")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("a")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("m")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("M")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }

    private static void _checkDetermineOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d")) {
        } else {
            System.err.println("Please specify datafile via -d option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            System.err.println("Please specify similarity method via -c option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("aL")) {
        } else {
            System.err.println("Please specify lower bound for alpha parameter via -aL option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("aH")) {
        } else {
            System.err.println("Please specify upper bound for alpha parameter via -aH option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("mL")) {
        } else {
            System.err.println("Please specify lower bound for mu parameter via -mL option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("mH")) {
        } else {
            System.err.println("Please specify upper bound for mu parameter via -mH option.");
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }

    private static void _checkTestOptions(CommandLine cmd, Options options) {
        HelpFormatter formatter = new HelpFormatter();
        if (cmd.hasOption("h")) {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("d")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);

            System.exit(0);
        }
        if (cmd.hasOption("c")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("o")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("aL")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("aH")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("mL")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("mH")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("M")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
        if (cmd.hasOption("t")) {
        } else {
            formatter.printHelp("java -jar jarfile.jar", options);
            System.exit(0);
        }
    }

    private static void constructNetwork(String[] args) {
        String IAM = "Construction";
        CommandLineParser parser = new BasicParser();
        Options options = buildConstructOptions();
        String pathIn = null;
        String Out = null;
        float alpha = 0.0f;
        float mu = 0.0f;
        float Mask = 0.0f;
        int threads = 0;
        String corr = null;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkConstructOptions(cmd, options);
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            pathIn = cmd.getOptionValue("d");
            Out = cmd.getOptionValue("o");
            corr = cmd.getOptionValue("c");
            alpha = Float.parseFloat(cmd.getOptionValue("a"));
            mu = Float.parseFloat(cmd.getOptionValue("m"));
            Mask = Float.parseFloat(cmd.getOptionValue("M"));

        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }

        String sep = "\t";

        int[] FileDimensions = new int[2];
        FileDimensions = _getFileDimensions(pathIn, sep);
        System.err.println("Beginning network construction...\n");
        System.err.println("Loading Data File\n");
        Operations.createDirectory(Out);
        ExpressionFrame DataFrame = _loadData(pathIn, FileDimensions, sep);
        String FrameOut = Out + "/InputExpression.matrix.tab";
        DataFrame.printMatrixToFile(FrameOut, sep);

        GCNMatrix CurrentMatrix = new GCNMatrix(FileDimensions[0], FileDimensions[0]);
        System.err.println("Calculating Similarity & Adjacency...\nMu : " + mu + "\nAlpha: " + alpha + "\n");
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", mu, alpha, threads,true);
        String ThisOut = Out + "/Adjacency.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        String MatrixOut = Out + "/Adjacency.cytoscape.raw.tab";
        CurrentMatrix.printMatrixToCytoscape(MatrixOut, "\t", Mask);
    
/*        
        String MatrixOut = Out + "/Adj.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);
*/
        
/*
        ThisOut = Out + "/TOM.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        MatrixOut = Out + "/TOM.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);
*/
        
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.HALF_UP);
        CurrentMatrix.calculateKs();
        double[] Return = CurrentMatrix.determineScaleFreeCritereon();
        double RSquared = Return[0];
        double Slope = Return[1];
        float mean = CurrentMatrix.getMeanK();
        if(Double.isNaN(RSquared)) RSquared=0.0d;
        RSquared=(Double.valueOf(df.format(RSquared)));
        mean = (Float.valueOf(df.format(mean)));
        System.out.println("Scale Free Criteron of resultant network: " + RSquared);
        System.out.println("Slope of scale free connectivity: " + Slope);
        System.out.println("Average Connectivity of resultant network: " + mean);
            
        System.err.println("Calculating TOM...\n");
        CurrentMatrix.calculateKs();
        CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);
        MatrixOut = Out + "/TOM.cytoscape.raw.tab";
        CurrentMatrix.printMatrixToCytoscape(MatrixOut, "\t", Mask);        
    //    String MatrixOut = Out + "/Network.Cytoscape.Raw.tab";
    //    CurrentMatrix.printMatrixToCytoscape(MatrixOut, "\t", Mask);
        System.out.println("Calculating clusters...");
        
        int MinSize = 50;
        Cluster Clustering = new Cluster(CurrentMatrix, 4);
        String ClustOut = Out+"/Clusters/";
        ArrayList<int[]> Clusters = Clustering.dynamicTreeCut(MinSize);
        _clustersToFile(CurrentMatrix, Clusters, MinSize, ClustOut);
        System.out.println("Done.");
        System.exit(0);
    }

    private static void determine(String[] args) {
        String IAM = "determine";
        CommandLineParser parser = new BasicParser();
        Options options = buildDetermineOptions();
        String pathIn = null;
        float alpha_Low = 0.0f;
        float alpha_High = 0.0f;
        float mu_Low = 0.0f;
        float mu_High= 0.0f;
        int threads = 0;
        String corr = null;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkDetermineOptions(cmd, options);
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            pathIn = cmd.getOptionValue("d");
            corr = cmd.getOptionValue("c");
            alpha_Low = Float.parseFloat(cmd.getOptionValue("aL"));
            alpha_High= Float.parseFloat(cmd.getOptionValue("aH"));
            mu_Low = Float.parseFloat(cmd.getOptionValue("mL"));
            mu_High= Float.parseFloat(cmd.getOptionValue("mH"));
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }

        String sep = "\t";

        int[] FileDimensions = new int[2];
        FileDimensions = _getFileDimensions(pathIn, sep);
        System.err.println("Beginning parameter testing\n\n");
        System.err.println("Loading Data File\n");
        ExpressionFrame DataFrame = _loadData(pathIn, FileDimensions, sep);

        GCNMatrix CurrentMatrix = new GCNMatrix(FileDimensions[0], FileDimensions[0]);
        System.err.println("Calculating Initial Similarity...\n");
        // passthrough, just calculate similarity
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", 0.0f, 0.0f, threads,true);
        //CurrentMatrix.maskMatrix(0.05f);
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.HALF_UP);
        String header="";
        String rsq_all = "";
        String mean_all = "";
        String slope_all= "";
        for(float MU=mu_Low;MU<=mu_High;MU+=0.05f){
            Float V = (Float.valueOf(df.format(MU)));
            String rsq_Line=String.valueOf(V);
            String mean_Line= rsq_Line;
            String slope_Line= rsq_Line;
            for(float A=alpha_Low;A<=alpha_High;A+=2.0f){
                if(MU==mu_Low){
                    header = header + "," + A;
                }
                System.err.println("Mu: "+MU+" Alpha: "+A);
                System.err.println("Calculating Adjacency...");
                GCNMatrix this_Adjacency = Operations.applySigmoid(CurrentMatrix, V, A, threads);
                this_Adjacency.calculateKs();
                double[] Returns = this_Adjacency.determineScaleFreeCritereon();
                Double RSquared = Returns[0];
                Double Slope = Returns[1];
                float mean = this_Adjacency.getMeanK();
                if(Double.isNaN(RSquared)) RSquared=0.0d;
                RSquared=(Double.valueOf(df.format(RSquared)));
                mean = (Float.valueOf(df.format(mean)));
                slope_Line = slope_Line + "," + Slope;
                rsq_Line = rsq_Line + "," + RSquared;
                mean_Line = mean_Line + "," + mean;
            }
            rsq_all = rsq_all + rsq_Line + "\n";
            mean_all = mean_all + mean_Line + "\n";
            slope_all = slope_all + slope_Line + "\n";
        }
        rsq_all = header + "\n" + rsq_all;
        mean_all = header + "\n" + mean_all;
        slope_all = header + "\n" + slope_all;
        System.out.println("Scale Free Criterion:\nCorrelation to Log-Log model");
        System.out.println(rsq_all);
        System.out.println();
        System.out.println("Regression Slope:");
        System.out.println(slope_all);
        System.out.println();
        System.out.println("Mean node connectivity");
        System.out.println(mean_all);
        System.out.println("Done.");
        System.exit(0);
    }

    private static void test(String[] args) {
        CommandLineParser parser = new BasicParser();
        Options options = buildConstructOptions();
        String pathIn = null;
        String Out = null;
        float alpha = 0.0f;
        float mu = 0.0f;
        float Mask = 0.0f;
        int threads = 0;
        String corr = null;
        threads = 16;

        Out = "test";
        corr = "pcc";
        alpha = 20;
        mu = 0.8f;
        Mask = 0.0f;

        GCNMatrix CurrentMatrix = new GCNMatrix(10, 10);
        System.err.println("Calculating Similarity & Adjacency...\n");
        float[] Row = {0.91f, 0.92f, 0.93f, 0.1f, 0.1f, 0.1f, 0.1f, 0.15f, 0.15f, 0.15f};
        CurrentMatrix.addRow(Row); // 0
        CurrentMatrix.addRow(Row); // 1
        CurrentMatrix.addRow(Row); // 2
        float[] Row1 = {0.1f, 0.1f, 0.1f, 0.91f, 0.92f, 0.93f, 0.94f, 0.15f, 0.15f, 0.15f};
        CurrentMatrix.addRow(Row1); // 3
        CurrentMatrix.addRow(Row1); // 4
        CurrentMatrix.addRow(Row1); // 5
        CurrentMatrix.addRow(Row1); // 6
        float[] Row2 = {0.15f, 0.15f, 0.15f, 0.1f, 0.1f, 0.1f, 0.1f, 0.91f, 0.92f, 0.93f};
        CurrentMatrix.addRow(Row2); // 7 
        CurrentMatrix.addRow(Row2); // 8
        CurrentMatrix.addRow(Row2); // 9

        //CurrentMatrix.calculateKs();
        //CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);
        //CurrentMatrix.maskMatrix(Mask);
        Cluster Clusters = new Cluster(CurrentMatrix, 4);
        System.out.println("Calculating clusters...");
        Clusters.dynamicTreeCut(25);
        //Dendrogram.getDendrogram(4); // Method for clustering 4 == average -- Parameterize
        //Dendrogram.getClusters(50); // 50 = min cluster size -- Parameterize
        System.out.println("Done.");
        System.exit(0);
    }

    private static void permute(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        String IAM = "Comparison";
        CommandLineParser parser = new BasicParser();
        Options options = buildPermuteOptions();
        String dir1;
        String dir2;
        String corr = "pcc";
        float mu1;
        float mu2;
        float alpha1;
        float alpha2;
        String out;
        int threads;
        int permutations;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkPermuteOptions(cmd, options);
            dir1 = cmd.getOptionValue("d1");
            dir2 = cmd.getOptionValue("d2");
            mu1  = Float.parseFloat(cmd.getOptionValue("m1"));
            mu2  = Float.parseFloat(cmd.getOptionValue("m2"));
            alpha1=Float.parseFloat(cmd.getOptionValue("a1"));
            alpha2=Float.parseFloat(cmd.getOptionValue("a2"));
            corr = cmd.getOptionValue("c");
            permutations = Integer.parseInt(cmd.getOptionValue("p"));
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            out = cmd.getOptionValue("o");
            Operations.createDirectory(out);

            String Exp1 = dir1 + "/InputExpression.matrix.tab";
            String Exp2 = dir2 + "/InputExpression.matrix.tab";
            String sep = "\t";
            System.err.println("Identifying run parameters...");
            int[] FD_1 = new int[2];
            FD_1 = _getFileDimensions(Exp1, sep);
            int[] FD_2 = new int[2];
            FD_2 = _getFileDimensions(Exp2, sep);

            if (FD_1[0] == FD_2[0]) {
            } else {
                System.err.println("Matrix files are not the same size");
                System.exit(0);
            }
            System.err.println("Comparing Networks...\n\n");
            System.err.println("Loading original Expression data... (1)");
            ExpressionFrame ExpF1 = _loadData(Exp1, FD_1, sep);
            System.err.println("Loading original Expression data... (2)");
            ExpressionFrame ExpF2 = _loadData(Exp2, FD_2, sep);
            
            float CUTOFF = 0.0f;
            float[] RESULT =  new float[FD_1[0]+1];
            if(permutations >0){
                System.err.println("Beginning permuation analysis...");
                CUTOFF = Operations.permuteData(ExpF1, ExpF2, permutations, out, corr, mu1, mu2, alpha1, alpha2, threads,true);
            }
            System.err.println("Permutations done. Obtained Cutoff of dTOM = " + CUTOFF +"\n");
            System.err.println("See "+out+"/PermutationDetails.tab for more details on calculated False Discovery Rates\n");
            
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }
    }
    
    private static void scaleFree(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        String IAM = "Comparison";
        CommandLineParser parser = new BasicParser();
        Options options = buildScalefreeOptions();
        String dir1;
        String dir2;
        String corr = "pcc";
        float mu1;
        float mu2;
        float alpha1;
        float alpha2;
        String out;
        int threads;
        int permutations;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkScalefreeOptions(cmd, options);
            dir1 = cmd.getOptionValue("d1");
            dir2 = cmd.getOptionValue("d2");
            mu1  = Float.parseFloat(cmd.getOptionValue("m1"));
            mu2  = Float.parseFloat(cmd.getOptionValue("m2"));
            alpha1=Float.parseFloat(cmd.getOptionValue("a1"));
            alpha2=Float.parseFloat(cmd.getOptionValue("a2"));
            corr = cmd.getOptionValue("c");
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            out = cmd.getOptionValue("o");
            Operations.createDirectory(out);

            String Exp1 = dir1 + "/InputExpression.matrix.tab";
            String Exp2 = dir2 + "/InputExpression.matrix.tab";
            String sep = "\t";
            System.err.println("Identifying run parameters...");
            int[] FD_1 = new int[2];
            FD_1 = _getFileDimensions(Exp1, sep);
            int[] FD_2 = new int[2];
            FD_2 = _getFileDimensions(Exp2, sep);

            if (FD_1[0] == FD_2[0]) {
            } else {
                System.err.println("Matrix files are not the same size");
                System.exit(0);
            }
            System.err.println("Comparing Networks...\n\n");
            System.err.println("Loading original Expression data... (1)");
            ExpressionFrame ExpF1 = _loadData(Exp1, FD_1, sep);
            System.err.println("Loading original Expression data... (2)");
            ExpressionFrame ExpF2 = _loadData(Exp2, FD_2, sep);
            
            float[] RESULT =  new float[2];
            System.err.println("Beginning scale-free analysis...");
            RESULT = Operations.determineCutoffSF(ExpF1, ExpF2, out, corr, mu1, mu2, alpha1, alpha2, threads);
            System.err.println("See files in "+out+"/ for more details on calculated scale free criteria\n");
            System.err.println("Cutoff for negative plasticity: " + RESULT[0]);
            System.err.println("Cutoff for positive plasticity: " + RESULT[1]);
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }
    }
    
    private static void compareNetworks(String[] args) {
        String IAM = "Comparison";
        CommandLineParser parser = new BasicParser();
        Options options = buildCompareOptions();
        String dir1;
        String dir2;
        String corr = "pcc";
        float mu1;
        float mu2;
        float alpha1;
        float alpha2;
        float CUTOFF_pos;
        float CUTOFF_neg;
        String out;
        int threads;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkCompareOptions(cmd, options);
            dir1 = cmd.getOptionValue("d1");
            dir2 = cmd.getOptionValue("d2");
            mu1  = Float.parseFloat(cmd.getOptionValue("m1"));
            mu2  = Float.parseFloat(cmd.getOptionValue("m2"));
            alpha1=Float.parseFloat(cmd.getOptionValue("a1"));
            alpha2=Float.parseFloat(cmd.getOptionValue("a2"));
            CUTOFF_pos=Float.parseFloat(cmd.getOptionValue("dpA"));
            CUTOFF_neg=Float.parseFloat(cmd.getOptionValue("dnA"));
            corr = cmd.getOptionValue("c");
            threads = Integer.parseInt(cmd.getOptionValue("t"));
            out = cmd.getOptionValue("o");
            Operations.createDirectory(out);

            String Exp1 = dir1 + "/InputExpression.matrix.tab";
            String Exp2 = dir2 + "/InputExpression.matrix.tab";
            String sep = "\t";
            System.err.println("Identifying run parameters...");
            int[] FD_1 = new int[2];
            FD_1 = _getFileDimensions(Exp1, sep);
            int[] FD_2 = new int[2];
            FD_2 = _getFileDimensions(Exp2, sep);

            if (FD_1[0] == FD_2[0]) {
            } else {
                System.err.println("Matrix files are not the same size");
                System.exit(0);
            }
            System.err.println("Comparing Networks...\n\n");
            System.err.println("Loading original Expression data... (1)");
            ExpressionFrame ExpF1 = _loadData(Exp1, FD_1, sep);
            System.err.println("Loading original Expression data... (2)");
            ExpressionFrame ExpF2 = _loadData(Exp2, FD_2, sep);

            System.err.println("Calculating Network Plasticity...");
            System.err.println("Calculation adjacencies...");
            GCNMatrix SNetworkA = Operations.calculateAdjacency(ExpF1, corr, "sigmoid", mu1, alpha1, threads,true);
            GCNMatrix SNetworkB = Operations.calculateAdjacency(ExpF2, corr, "sigmoid", mu2, alpha2, threads,true);
            GCNMatrix USNetworkA = Operations.calculateAdjacency(ExpF1, corr, "sigmoid", mu1, alpha1, threads,false);
            GCNMatrix USNetworkB = Operations.calculateAdjacency(ExpF2, corr, "sigmoid", mu2, alpha2, threads,false);
            SNetworkA.calculateKs();
            SNetworkB.calculateKs();
            USNetworkA.calculateKs();
            USNetworkB.calculateKs();
            String[] names = SNetworkA.getRowNames();
            
            System.err.println("Calculating final plasticity network...");
            
            GCNMatrix SDifference = Operations.calculateDifferenceThreaded(SNetworkA, SNetworkB, threads);
            /*
            NetworkA = Operations.calculateTOM(NetworkA, threads);
            NetworkB = Operations.calculateTOM(NetworkB, threads)
            float[] MWW = Operations.compareNetworksViaMWW(NetworkA,NetworkB);
            float[] averagePlasticity = Difference.getCurrentAverageEdgeStrength();
            _MWWToFile(MWW,averagePlasticity,names,out,"NodePlasticity.MWW.tab");
            */
            String O3 = out + "/Adjacency.unmasked.plasticity.cytoscape.tab";
            SDifference.printMatrixToCytoscape(O3, "\t", 0.01f);
            
            //System.err.println("Masking plasticity network based on significance cutoff: " + CUTOFF);
            //Difference.maskOutside(CUTOFF);
            DecimalFormat df = new DecimalFormat("#.###");
            df.setRoundingMode(RoundingMode.HALF_UP);
            //Difference = Operations.calculateTOM(Difference, threads);
            // Clustering
            System.err.println("Finding negative plasticity network... cutoff: " + CUTOFF_neg);
            //SDifference.maskAbove(CUTOFF_neg);
            SDifference = Operations.findPlasticity(SDifference,CUTOFF_neg,SNetworkA,SNetworkB,"negative");
            SDifference.calculateKs();

            O3 = out + "/Adjacency.negPlasticity.cytoscape.tab";
            SDifference.printMatrixToCytoscape(O3, "\t", 0.01f);
            double[] Return = SDifference.determineScaleFreeCritereon();
            double RSquared = Return[0];
            double Slope = Return[1];
            float mean = SDifference.getMeanK();
            if(Double.isNaN(RSquared)) RSquared=0.0d;
            RSquared=(Double.valueOf(df.format(RSquared)));
            mean = (Float.valueOf(df.format(mean)));
            System.out.println("Metrics for Negative Plasticity matrix:");
            System.out.println("Scale Free Criteron of resultant plasticity network: " + RSquared);
            System.out.println("Slope of scale free connectivity: " + Slope);
            System.out.println("Average Connectivity of resultant plasticity network: " + mean);
            System.err.println("Clustering negative plasticity...");
            //GCNMatrix USDifference = Operations.calculateDifferenceThreaded(USNetworkA, USNetworkB, threads); // yeah this is a hack
            //USDifference.maskAbove(CUTOFF_neg);
            //USDifference.calculateKs();
            SDifference.prepareForTOM();
            SDifference = Operations.calculateTOM(SDifference, threads);
            O3 = out + "/TOM.negPlasticity.cytoscape.tab";
            SDifference.printMatrixToCytoscape(O3, "\t", 0.01f);
            int MinSize = 50;
            String ClustOut = out+"/Clusters_Neg/";
            Cluster Clustering = new Cluster(SDifference, 4);
            ArrayList<int[]> Clusters = Clustering.dynamicTreeCut(MinSize);
            _clustersToFile(SDifference, Clusters, MinSize, ClustOut);            
            // For some reason, using the same variable and overwriting below did not work
            // Not sure why... even copy constructors don't seem to work. hacking.
            SDifference = Operations.calculateDifferenceThreaded(SNetworkA, SNetworkB,threads);
            System.err.println("Finding positive plasticity network... cutoff: " + CUTOFF_pos);
            //SDifference.maskBelow(CUTOFF_pos); // MaskedDif is now the pos matrix
            O3 = out + "/Adjacency.posPlasticity.cytoscape.tab";
            SDifference = Operations.findPlasticity(SDifference,CUTOFF_pos,SNetworkA,SNetworkB,"positive");
            SDifference.printMatrixToCytoscape(O3, "\t", 0.01f);
            SDifference.calculateKs();
            Return = SDifference.determineScaleFreeCritereon();
            RSquared = Return[0];
            Slope = Return[1];
            mean = SDifference.getMeanK();
            if(Double.isNaN(RSquared)) RSquared=0.0d;
            RSquared=(Double.valueOf(df.format(RSquared)));
            mean = (Float.valueOf(df.format(mean)));
            System.out.println("Metrics for Positive Plasticity matrix:");
            System.out.println("Scale Free Criteron of resultant plasticity network: " + RSquared);
            System.out.println("Slope of scale free connectivity: " + Slope);
            System.out.println("Average Connectivity of resultant plasticity network: " + mean);
            System.err.println("Clustering positive plasticity...");
            //USDifference = Operations.calculateDifferenceThreaded(USNetworkA, USNetworkB, threads);// yeah this is a hack
            //USDifference.calculateKs();
            //USDifference.maskBelow(CUTOFF_pos);
            SDifference.prepareForTOM();
            SDifference = Operations.calculateTOM(SDifference, threads);
            O3 = out + "/TOM.posPlasticity.cytoscape.tab";
            SDifference.printMatrixToCytoscape(O3, "\t", 0.01f);
            ClustOut = out+"/Clusters_Pos/";
            Clustering = new Cluster(SDifference, 4);
            Clusters = Clustering.dynamicTreeCut(MinSize);
            _clustersToFile(SDifference, Clusters, MinSize, ClustOut);
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }
    }

    private static void baseOptions(String[] args) {
        Options options = buildOptions();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("java -jar jarfile.jar", options);
        System.exit(0);
    }

    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        if (args.length == 0) {
            baseOptions(args);
        } else {
            switch (args[0]) {
                case "construct":
                    constructNetwork(args);
                    break;
                case "speedtest":
                    speedtest(args);
                    break;
                case "compare":
                    compareNetworks(args);
                    break;
                case "test":
                    test(args);
                    break;
                case "determine":
                    determine(args);
                    break;
                case "permute":
                    permute(args);
                    break;
                case "scalefree":
                    scaleFree(args);
                    break;
                default:
                    baseOptions(args);
                    break;
            }
        }
    }

    private static Options testOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("datafile");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values");
        Option datafile = OptionBuilder.create("d");

        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric to use ('gini' or 'pcc')");
        Option corr = OptionBuilder.create("c");

        OptionBuilder.withArgName("alpha");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha parameter for sigmoid adjacency calculation (i.e. 20)");
        Option alpha = OptionBuilder.create("a");

        OptionBuilder.withArgName("mu");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu parameter for sigmoid adjacency calculation (i.e. 0.8)");
        Option mu = OptionBuilder.create("m");

        OptionBuilder.withArgName("Mask");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("parameter for masking matricies - values below this will be set to zero");
        Option Mask = OptionBuilder.create("M");

        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Temporary Directory");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of Compute Threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(datafile);
        options.addOption(alpha);
        options.addOption(corr);
        options.addOption(mu);
        options.addOption(Mask);
        options.addOption(output);
        options.addOption(threads);

        return options;
    }

    private static Options buildOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("similarity");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("construct a similarity matrix and print it to a file");
        Option similarity = OptionBuilder.create("similarity");

        OptionBuilder.withArgName("construct");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("construct a network via application of a topological overlap calculation");
        Option construct = OptionBuilder.create("construct");

        OptionBuilder.withArgName("permute");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("perform background permutations to determine significance cutoff");
        Option permute = OptionBuilder.create("permute");
        
        OptionBuilder.withArgName("scalefree");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Determine plasticity cutoff via conformation to scale free criteria");
        Option scalefree = OptionBuilder.create("scalefree");
        
        OptionBuilder.withArgName("compare");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("compare several network matricies");
        Option compare = OptionBuilder.create("compare");

        OptionBuilder.withArgName("determine");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Determine best parameters for dataset");
        Option determine = OptionBuilder.create("determine");        
        
        OptionBuilder.withArgName("test");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("test routines");
        Option test = OptionBuilder.create("test");

        options.addOption(help);
        options.addOption(determine);
        options.addOption(construct);
        options.addOption(permute);
        options.addOption(scalefree);
        options.addOption(compare);
        options.addOption(similarity);
        options.addOption(test);
        

        return options;
    }

    private static Options buildCompareOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("matrix1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 1");
        Option dir1 = OptionBuilder.create("d1");

        OptionBuilder.withArgName("matrix2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 2");
        Option dir2 = OptionBuilder.create("d2");

        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("output file");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("mu1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 1");
        Option mu1 = OptionBuilder.create("m1");

        OptionBuilder.withArgName("mu2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 2");
        Option mu2 = OptionBuilder.create("m2");
        
        OptionBuilder.withArgName("alpha1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 1");
        Option alpha1 = OptionBuilder.create("a1");
        
        OptionBuilder.withArgName("alpha2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 2");
        Option alpha2 = OptionBuilder.create("a2");
        
        OptionBuilder.withArgName("positive delta-Adjacency");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("delta-Adjacency Overlap Cutoff for positive plasticity determination");
        Option pos_delta_Adj = OptionBuilder.create("dpA");
        
        OptionBuilder.withArgName("negative delta-Adjacency");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("delta-Adjacency Overlap Cutoff for negative plasticity determination");
        Option neg_delta_Adj = OptionBuilder.create("dnA");

        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric used for both networks");
        Option c = OptionBuilder.create("c");
        
        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("number of processing threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(dir1);
        options.addOption(dir2);
        options.addOption(mu1);
        options.addOption(mu2);
        options.addOption(alpha1);
        options.addOption(alpha2);
        options.addOption(neg_delta_Adj);
        options.addOption(pos_delta_Adj);
        options.addOption(c);
        options.addOption(output);
        options.addOption(threads);
        return options;
    }
    
        private static Options buildPermuteOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("matrix1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 1");
        Option dir1 = OptionBuilder.create("d1");

        OptionBuilder.withArgName("matrix2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 2");
        Option dir2 = OptionBuilder.create("d2");

        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("output file");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("mu1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 1");
        Option mu1 = OptionBuilder.create("m1");

        OptionBuilder.withArgName("mu2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 2");
        Option mu2 = OptionBuilder.create("m2");
        
        OptionBuilder.withArgName("alpha1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 1");
        Option alpha1 = OptionBuilder.create("a1");
        
        OptionBuilder.withArgName("alpha2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 2");
        Option alpha2 = OptionBuilder.create("a2");
        
        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric used for both networks");
        Option c = OptionBuilder.create("c");
        
        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("number of processing threads");
        Option threads = OptionBuilder.create("t");

        OptionBuilder.withArgName("permutations");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("number of permutations to determine significance");
        Option permutations = OptionBuilder.create("p");

        options.addOption(help);
        options.addOption(dir1);
        options.addOption(dir2);
        options.addOption(mu1);
        options.addOption(mu2);
        options.addOption(alpha1);
        options.addOption(alpha2);
        options.addOption(c);
        options.addOption(output);
        options.addOption(threads);
        options.addOption(permutations);
        return options;
    }
    
    private static Options buildScalefreeOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("matrix1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 1");
        Option dir1 = OptionBuilder.create("d1");

        OptionBuilder.withArgName("matrix2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing run files for network 2");
        Option dir2 = OptionBuilder.create("d2");

        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("output file");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("mu1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 1");
        Option mu1 = OptionBuilder.create("m1");

        OptionBuilder.withArgName("mu2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu for network 2");
        Option mu2 = OptionBuilder.create("m2");
        
        OptionBuilder.withArgName("alpha1");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 1");
        Option alpha1 = OptionBuilder.create("a1");
        
        OptionBuilder.withArgName("alpha2");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha for network 2");
        Option alpha2 = OptionBuilder.create("a2");
        
        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric used for both networks");
        Option c = OptionBuilder.create("c");
        
        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("number of processing threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(dir1);
        options.addOption(dir2);
        options.addOption(mu1);
        options.addOption(mu2);
        options.addOption(alpha1);
        options.addOption(alpha2);
        options.addOption(c);
        options.addOption(output);
        options.addOption(threads);
        return options;
    }
 
    private static Options buildConstructOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("datafile");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values");
        Option datafile = OptionBuilder.create("d");

        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric to use ('gini','pcc','spearman')");
        Option corr = OptionBuilder.create("c");

        OptionBuilder.withArgName("alpha");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha parameter for sigmoid adjacency calculation (i.e. 20)");
        Option alpha = OptionBuilder.create("a");

        OptionBuilder.withArgName("mu");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu parameter for sigmoid adjacency calculation (i.e. 0.8)");
        Option mu = OptionBuilder.create("m");

        OptionBuilder.withArgName("Mask");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("parameter for masking matricies - values below this will be set to zero");
        Option Mask = OptionBuilder.create("M");

        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Temporary Directory");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of Compute Threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(datafile);
        options.addOption(alpha);
        options.addOption(corr);
        options.addOption(mu);
        options.addOption(Mask);
        options.addOption(output);
        options.addOption(threads);

        return options;
    }

    private static Options buildTestOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("datafile");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values");
        Option datafile = OptionBuilder.create("d");

        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric to use ('gini' or 'pcc')");
        Option corr = OptionBuilder.create("c");

        OptionBuilder.withArgName("alpha low");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Lower bound for alpha parameter - must be a mutliple of 2");
        Option alpha_low = OptionBuilder.create("aL");

        OptionBuilder.withArgName("alpha high");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Upper bound for alpha parameter - must be a multiple of 2");
        Option alpha_high = OptionBuilder.create("aH");
        
        OptionBuilder.withArgName("mu low");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Lower bound for mu parameter - must be a multiple of .1, on the interval [0,1]");
        Option mu_low = OptionBuilder.create("mL");

        OptionBuilder.withArgName("mu high");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Upper bound for mu parameter - must be a multiple of .1, on the interval [0,1]");
        Option mu_high = OptionBuilder.create("mH");             
        
        OptionBuilder.withArgName("output");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Temporary Directory");
        Option output = OptionBuilder.create("o");

        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of Compute Threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(datafile);
        options.addOption(alpha_low);
        options.addOption(alpha_high);
        options.addOption(corr);
        options.addOption(mu_low);
        options.addOption(mu_high);
        options.addOption(output);
        options.addOption(threads);

        return options;
    }
    
    private static Options buildDetermineOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("datafile");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Data frame, tab delimited, with header, of per-gene, per-condition expression values");
        Option datafile = OptionBuilder.create("d");

        OptionBuilder.withArgName("correlation");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("correlation metric to use ('gini','pcc' or 'spearman')");
        Option corr = OptionBuilder.create("c");

        OptionBuilder.withArgName("alpha low");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Lower bound for alpha parameter - must be a mutliple of 2");
        Option alpha_low = OptionBuilder.create("aL");

        OptionBuilder.withArgName("alpha high");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Upper bound for alpha parameter - must be a multiple of 2");
        Option alpha_high = OptionBuilder.create("aH");
        
        OptionBuilder.withArgName("mu low");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Lower bound for mu parameter - must be a multiple of .1, on the interval [0,1]");
        Option mu_low = OptionBuilder.create("mL");

        OptionBuilder.withArgName("mu high");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Upper bound for mu parameter - must be a multiple of .1, on the interval [0,1]");
        Option mu_high = OptionBuilder.create("mH");             

        OptionBuilder.withArgName("threads");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Number of Compute Threads");
        Option threads = OptionBuilder.create("t");

        options.addOption(help);
        options.addOption(datafile);
        options.addOption(alpha_low);
        options.addOption(alpha_high);
        options.addOption(corr);
        options.addOption(mu_low);
        options.addOption(mu_high);
        options.addOption(threads);

        return options;
    }

    private static int[] _getFileDimensions(String pathIn, String sep) {
        int[] dimensions = new int[2];
        // pre-declaring sizes allows use of non-dynamic double[][] instead of nested ArrayLists. 
        // performance gain over ArrayList per-entry is very small, but we have N(N-1)/2 entries...
        try {
            File file = null;

            try {
                file = new File(pathIn);
            } catch (NullPointerException e) {
                System.err.println("No file found to read.\n");
                System.exit(0);
            }

            Scanner scanner = new Scanner(file);
            String header[] = scanner.nextLine().split(sep);
            dimensions[0] = 0; // Frame height (minus header)
            dimensions[1] = header.length - 1;  // Frame width (minus rowID)
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                String[] Line = line.split(sep);
                int this_width = Line.length - 1;
                if (this_width != dimensions[1]) {
                    _fileDimErr();
                }
                dimensions[0] += 1;
            }
            scanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            System.exit(0);
        } finally {
        }
        return dimensions;
    }

    private static void _fileDimErr() {
        System.err.println("data file is not a rectilinear data file");
        System.exit(0);
    }

    private static ExpressionFrame _loadData(String pathIn, int[] Dims, String sep) {
        ExpressionFrame Expression = new ExpressionFrame(Dims[0], Dims[1]);
        try {
            File file = null;
            try {
                file = new File(pathIn);
            } catch (NullPointerException e) {
                System.err.println("No file found to read.\n");
                System.exit(0);
            }

            Scanner scanner = new Scanner(file);
            String[] header = scanner.nextLine().split(sep);
            String[] loci = new String[Dims[0]];
            Expression.setColumnNames(header);
            int it = 0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                String[] Line = line.split(sep);
                loci[it] = Line[0];
                float[] data = new float[Dims[1]];
                for (int i = 1; i < Line.length; i++) {
                    try {
                        int I = i - 1;
                        float value = Float.parseFloat(Line[i]);
                        data[I] = value;
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                    }
                }
                it++;
                Expression.addRow(data);
            }
            Expression.setRowNames(loci);
            scanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return Expression;
    }

    private static GCNMatrix _loadNetwork(String pathIn, int[] Dims, String sep) {
        GCNMatrix Matrix = new GCNMatrix(Dims[0], Dims[1]);
        // TODO : Completely change this datatype:
        // GCNMatricies and Expression Frames are not at all the same
        //System.err.println("Loading file of " + Dims[0] + " by " + Dims[1] + "\n");
        try {
            File file = null;
            try {
                file = new File(pathIn);
            } catch (NullPointerException e) {
                System.err.println("No file found to read.\n");
                System.exit(0);
            }

            Scanner scanner = new Scanner(file);
            String[] header = scanner.nextLine().split(sep);
            String[] loci = new String[Dims[0]];
            Matrix.setColumnNames(header);
            int it = 0;
            while (scanner.hasNextLine()) {
                String line = scanner.nextLine();
                String[] Line = line.split(sep);
                loci[it] = Line[0];
                float[] data = new float[Dims[1]];
                for (int i = 1; i < Line.length; i++) {
                    try {
                        int I = i - 1;
                        float value = Float.parseFloat(Line[i]);
                        data[I] = value;
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                    }
                }
                it++;
                Matrix.addRow(data);
            }
            Matrix.setRowNames(loci);
            scanner.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return Matrix;
    }

    private static void _clustersToFile(GCNMatrix Similarity, ArrayList<int[]> Clusters, int MinSize, String ClustDir) {
        Iterator<int[]> it = Clusters.iterator();
        int iter = 1;
        int count = 0;
        while (it.hasNext()) {
            int[] cluster = it.next();
            if (cluster.length < MinSize) {
                continue;
            }
            count++;
            //System.out.println("Final cluster size: " + cluster.length);
            //String ClustDir = OutDir + "/Clusters/";
            Operations.createDirectory(ClustDir);
            String nPath = ClustDir + "/" + "Cluster." + iter + ".txt";
            iter++;
            try {
                PrintWriter writer = new PrintWriter(nPath, "UTF-8");
                for (int i = 0; i < cluster.length; i++) {
                    int node = cluster[i];
                    String name = Similarity.getRowName(node);
                    writer.println(name);
                }
                writer.close();
            } catch (Exception e) {
                //
            }
        }
        System.err.println("Found " + count + " clusters.");
    }

    private static void _cTOMsToFile(float[] rTOMs, String[] names, String OutDir, String fileName) {
        String Path = OutDir + "/" + fileName ;
        try {
            PrintWriter writer = new PrintWriter(Path, "UTF-8");
            for (int i = 1; i < names.length; i++) {
                int index = i - 1;
                String line = names[index] + "\t" + rTOMs[index];
                writer.println(line);
            }
            writer.close();
        } catch (Exception e) {
            //
        }

    }
    private static void _MWWToFile(float[] rTOMs, float[] average, String[] names, String OutDir, String fileName) {
        String Path = OutDir + "/" + fileName ;
        try {
            PrintWriter writer = new PrintWriter(Path, "UTF-8");
            for (int i = 1; i < names.length; i++) {
                int index = i - 1;
                String line = names[index] + "\t" + rTOMs[index] + "\t" + average[index];
                writer.println(line);
            }
            writer.close();
        } catch (Exception e) {
            //
        }

    }
}
