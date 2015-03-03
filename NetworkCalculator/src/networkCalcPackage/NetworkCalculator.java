package networkCalcPackage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
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
        double Mask = 0.0;
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
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", mu, alpha, threads);
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
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", mu, alpha, threads);
/*
        String ThisOut = Out + "/Adjacency.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        String MatrixOut = Out + "/Adj.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);
*/
        System.err.println("Calculating TOM...\n");
        CurrentMatrix.calculateKs();
        CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);
/*
        ThisOut = Out + "/TOM.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        MatrixOut = Out + "/TOM.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);
*/
        
        String MatrixOut = Out + "/Network.Cytoscape.Raw.tab";
        CurrentMatrix.printMatrixToCytoscape(MatrixOut, "\t", Mask);
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
        float Mask = 0.0f;
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
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", 0.0f, 0.0f, threads);
        //CurrentMatrix.maskMatrix(0.05f);
        DecimalFormat df = new DecimalFormat("#.###");
        df.setRoundingMode(RoundingMode.HALF_UP);
        String header="";
        String rsq_all = "";
        String mean_all = "";
        for(float MU=mu_Low;MU<=mu_High;MU+=0.05f){
            Float V = (Float.valueOf(df.format(MU)));
            String rsq_Line=String.valueOf(V);
            String mean_Line= rsq_Line;
            for(float A=alpha_Low;A<=alpha_High;A+=2.0f){
                if(MU==mu_Low){
                    header = header + "," + A;
                }
                GCNMatrix this_Adjacency = Operations.applySigmoid(CurrentMatrix, V, A, threads);
                this_Adjacency.calculateKs();
                
                double RSquared = this_Adjacency.determineScaleFreeCritereon();
                float mean = this_Adjacency.getMeanK();
                if(Double.isNaN(RSquared)) RSquared=0.0d;
                RSquared=(Double.valueOf(df.format(RSquared)));
                mean = (Float.valueOf(df.format(mean)));
                
                //System.out.println("Mu: " + MU + " Alpha: " + A + " RSQ " + RSquared);
                rsq_Line = rsq_Line + "," + RSquared;
                mean_Line = mean_Line + "," + mean;
            }
            rsq_all = rsq_all + rsq_Line + "\n";
            mean_all = mean_all + mean_Line + "\n";
        }
        rsq_all = header + "\n" + rsq_all;
        mean_all = header + "\n" + mean_all;
        System.out.println("Scale Free Criterion:\nCorrelation to Log-Log model");
        System.out.println(rsq_all);
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
        String out;
        int threads;
        int permutations;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkCompareOptions(cmd, options);
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
            /*
             * TODO : need to implement a few things 
             * 1: print and load expression 
             */

            String Exp1 = dir1 + "/InputExpression.matrix.tab";
            String Exp2 = dir2 + "/InputExpression.matrix.tab";
            float mu = 0.8f;
            float alpha = 20f;
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
            
            float pmu = (mu1+mu2)/2;
            float pa  = (alpha1+alpha2)/2;
            float CUTOFF = 0.0f;
            float[] RESULT =  new float[FD_1[0]+1];
            if(permutations >0){
                System.err.println("Beginning permuation analysis...");
                RESULT = Operations.permuteData(ExpF1, ExpF2, permutations, out, corr, pmu, pa, threads);
                CUTOFF = RESULT[0];
            }
            System.err.println("Permutations done. Obtained Cutoff of dTOM = " + CUTOFF);

            System.err.println("Calculating actual values...");
            GCNMatrix NetworkA = Operations.calculateAdjacency(ExpF1, corr, "sigmoid", mu1, alpha1, threads);
            GCNMatrix NetworkB = Operations.calculateAdjacency(ExpF2, corr, "sigmoid", mu2, alpha2, threads);
            NetworkA.calculateKs();
            NetworkB.calculateKs();
            //float[] rcTOMs = Operations.compareNetworksViaTOM(NetworkA, NetworkB);
            float[] rcTOMs = Operations.compareNetworksViaAverage(NetworkA, NetworkB);
            String[] names = NetworkA.getRowNames();
            if(permutations >0) _cTOMsToFile(rcTOMs, RESULT, names, out);
            String ThisOut = out + "/dTOM.dist.tab";
            NetworkA = Operations.calculateTOM(NetworkA, threads);
            NetworkB = Operations.calculateTOM(NetworkB, threads);
            GCNMatrix Difference = Operations.calculateDifference(NetworkA, NetworkB);
            
            DecimalFormat df = new DecimalFormat("#.###");
            df.setRoundingMode(RoundingMode.HALF_UP);
            Difference.calculateKs();
            double RSquared = Difference.determineScaleFreeCritereon();
            float mean = Difference.getMeanK();
            if(Double.isNaN(RSquared)) RSquared=0.0d;
            RSquared=(Double.valueOf(df.format(RSquared)));
            mean = (Float.valueOf(df.format(mean)));
            System.out.println("Scale Free Criteron of resultant plasticity matrix: " + RSquared);
            System.out.println("Average Connectivity of resultant plasticity matrix: " + mean);
            
            Difference.generateDistributionToFile(ThisOut);
            String O3 = out + "/Cytoscape.sigEdge.tab";
            Difference.printMatrixToCytoscape(O3, "\t", CUTOFF);
            O3 = out + "/Cytoscape.raw.tab";
            Difference.printMatrixToCytoscape(O3, "\t", 0.01f);
            
            // Clustering
            System.err.println("Clustering negative plasticity...");
            Difference.maskAbove(0.0f);
            int MinSize = 50;
            String ClustOut = out+"/Clusters_Neg/";
            Cluster Clustering = new Cluster(Difference, 4);
            ArrayList<int[]> Clusters = Clustering.dynamicTreeCut(MinSize);
            _clustersToFile(Difference, Clusters, MinSize, ClustOut);            
            // For some reason, using the same variable and overwriting below did not work
            // Not sure why... even copy constructors don't seem to work. hacking.
            Difference = Operations.calculateDifference(NetworkA, NetworkB);
            System.err.println("Clustering positive plasticity...");
            Difference.maskBelow(0.0f); // MaskedDif is now the pos matrix
            ClustOut = out+"/Clusters_Pos/";
            Clustering = new Cluster(Difference, 4);
            Clusters = Clustering.dynamicTreeCut(MinSize);
            _clustersToFile(Difference, Clusters, MinSize, ClustOut);

        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");

            System.exit(0);
        }

        // *** TODO: make network construction method
    }

    private static void clusterNetwork(String[] args) {
        CommandLineParser parser = new BasicParser();
        Options options = buildClusterOptions();
        // **** TODO -- again, what is this for??
        String directory = null;
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
            directory = cmd.getOptionValue("d");
        } catch (ParseException exp) {
            System.err.println("Problem parsing arguments:\n" + exp.getMessage());
            System.err.println("Exiting...\n");
            System.exit(0);
        }
        String pathIn = directory + "/tom.matrix";
        String sep = ",";
        int[] FileDimensions = new int[2];
        FileDimensions = _getFileDimensions(pathIn, sep);

        System.err.println("Loading Data File\n");

        GCNMatrix DataFrame = new GCNMatrix(FileDimensions[0], FileDimensions[1]);
        DataFrame = _loadNetwork(pathIn, FileDimensions, sep);

        System.exit(0);
    }

    private static void baseOptions(String[] args) {
        Options options = buildOptions();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("java -jar jarfile.jar", options);
        System.exit(0);
    }

    public static void main(String[] args) {
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
                case "view": // TODO make into query
                    clusterNetwork(args);
                    break;
                case "test":
                    test(args);
                    break;
                case "determine":
                    determine(args);
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

        OptionBuilder.withArgName("compare");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("compare several network matricies");
        Option compare = OptionBuilder.create("compare");

        OptionBuilder.withArgName("view");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("query a single network for node properties");
        Option view = OptionBuilder.create("view");

        OptionBuilder.withArgName("determine");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Determine best parameters for dataset");
        Option determine = OptionBuilder.create("determine");        
        
        OptionBuilder.withArgName("test");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("test routines");
        Option test = OptionBuilder.create("test");

        options.addOption(help);
        options.addOption(construct);
        options.addOption(similarity);
        options.addOption(compare);
        options.addOption(view);
        options.addOption(test);
        options.addOption(determine);

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

    private static Options buildClusterOptions() {
        Options options = new Options();
        Option help = new Option("h", "print this message");

        OptionBuilder.withArgName("Analysis Directory");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("Directory containing the outputs of the Construct command");
        Option directory = OptionBuilder.create("d");

        options.addOption(help);
        options.addOption(directory);

        System.out.println("This method is not yet implemented\n");
        System.exit(0);
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
            System.out.println("Final cluster size: " + cluster.length);
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

    private static void _cTOMsToFile(float[] rTOMs, float[] ratios, String[] names, String OutDir) {
        String Path = OutDir + "/CrossNetwork.TOM.tab";
        try {
            PrintWriter writer = new PrintWriter(Path, "UTF-8");
            for (int i = 1; i < ratios.length; i++) {
                int index = i - 1;
                String line = names[index] + "\t" + rTOMs[index] + "\t" + ratios[i];
                writer.println(line);
            }
            writer.close();
        } catch (Exception e) {
            //
        }

    }
}
