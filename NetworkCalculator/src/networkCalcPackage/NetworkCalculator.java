package networkCalcPackage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

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

        String ThisOut = Out + "/Adjacency.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        String MatrixOut = Out + "/Adj.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);

        System.err.println("Calculating TOM...\n");
        CurrentMatrix.calculateKs();
        CurrentMatrix = Operations.calculateTOM(CurrentMatrix, threads);

        ThisOut = Out + "/TOM.dist.tab";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        MatrixOut = Out + "/TOM.matrix.tab";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);

        System.out.println("Calculating clusters...");
        int MinSize = 50;
        Cluster Clustering = new Cluster(CurrentMatrix, 4);
        ArrayList<int[]> Clusters = Clustering.dynamicTreeCut(MinSize);
        _clustersToFile(CurrentMatrix, Clusters, MinSize, Out);
        // Cluster.getClusters(50); // 50 = min cluster size -- Parameterize
        System.out.println("Done.");
        System.exit(0);
    }

    private static void testMetrics(String[] args) {
        String IAM = "test";
        CommandLineParser parser = new BasicParser();
        Options options = buildTestOptions();
        String pathIn = null;
        String Out = null;
        float alpha = 0.0f;
        float mu = 0.0f;
        float Mask = 0.0f;
        int threads = 0;
        String corr = null;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkTestOptions(cmd, options);
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
        System.err.println("Beginning testing run...\n\n");
        System.err.println("Loading Data File\n");
        Operations.createDirectory(Out);
        ExpressionFrame DataFrame = _loadData(pathIn, FileDimensions, sep);
        String FrameOut = Out + "/InputExpression.matrix.tab";
        DataFrame.printMatrixToFile(FrameOut, sep);

        GCNMatrix CurrentMatrix = new GCNMatrix(FileDimensions[0], FileDimensions[0]);
        System.err.println("Calculating Similarity & Adjacency...\nMu : " + mu + "\nAlpha: " + alpha + "\n");
        CurrentMatrix = Operations.calculateAdjacency(DataFrame, corr, "sigmoid", mu, alpha, threads);
        // TODO : implement K distribution printing
        String ThisOut = Out + "/Test.dist.txt";
        CurrentMatrix.generateDistributionToFile(ThisOut);
        String MatrixOut = Out + "/Adj.matrix.txt";
        CurrentMatrix.printMatrixToFile(MatrixOut, sep);
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
        String out;
        int threads;
        int permutations;
        try {
            CommandLine cmd = parser.parse(options, args);
            _checkCompareOptions(cmd, options);
            dir1 = cmd.getOptionValue("d1");
            dir2 = cmd.getOptionValue("d2");
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
            String matrix1 = dir1 + "/Adj.matrix.tab";
            String matrix2 = dir2 + "/Adj.matrix.tab";
            float mu = 0.8f;
            float alpha = 20f;
            String sep = "\t";

            int[] FD_1 = new int[2];
            FD_1 = _getFileDimensions(matrix1, sep);
            int[] FD_2 = new int[2];
            FD_2 = _getFileDimensions(matrix2, sep);

            if ((FD_1[0] == FD_2[0]) & (FD_1[1] == FD_2[1])) {
            } else {
                System.err.println("Matrix files are not the same size");
                System.exit(0);
            }
            System.err.println("Comparing Networks...\n\n");
            System.err.println("Loading original Expression data... (1)");
            int[] ExpDim = _getFileDimensions(Exp1, sep);
            ExpressionFrame ExpF1 = _loadData(Exp1, ExpDim, sep);
            System.err.println("Loading original Expression data... (2)");
            ExpDim = _getFileDimensions(Exp2, sep);
            ExpressionFrame ExpF2 = _loadData(Exp2, ExpDim, sep);
            String corr = "pcc";
            System.err.println("Beginning permuation analysis...");
            float[] RESULT = Operations.permuteData(ExpF1, ExpF2, permutations, out, corr, mu, alpha, threads);
            float CUTOFF = RESULT[0];
            System.err.println("Permutations done. Obtained Cutoff of dTOM = " + CUTOFF);

            System.err.println("Calculating actual values...");
            GCNMatrix NetworkA = Operations.calculateAdjacency(ExpF1, corr, "sigmoid", mu, alpha, threads);
            GCNMatrix NetworkB = Operations.calculateAdjacency(ExpF2, corr, "sigmoid", mu, alpha, threads);
            NetworkA.calculateKs();
            NetworkB.calculateKs();
            //float[] rcTOMs = Operations.compareNetworksViaTOM(NetworkA, NetworkB);
            float[] rcTOMs = Operations.compareNetworksViaAverage(NetworkA, NetworkB);
            String[] names = NetworkA.getRowNames();
            _cTOMsToFile(rcTOMs, RESULT, names, out);
            // TODO : add back in self-wise TOM
            String O2 = out + "/Selfwise.actual.jpeg";
            String ThisOut = out + "/dTOM.dist.tab";
            GCNMatrix Difference = Operations.calculateDifference(NetworkA, NetworkB);
            Difference.generateDistributionToFile(ThisOut);
            Operations.generateHistogramHM(Difference, O2, "Cross-network Selfwise Topological Overlap Zm vs Sv", "selfwise TOM", "Count", false);
            String O3 = out + "/Cytoscape.sigEdge.tab";
            Difference.printMatrixToCytoscape(O3, "\t", CUTOFF);
            O3 = out + "/Cytoscape.raw.tab";
            Difference.printMatrixToCytoscape(O3, "\t", 0.0f);

            int MinSize = 50;
            Cluster Clustering = new Cluster(Difference, 4);
            ArrayList<int[]> Clusters = Clustering.dynamicTreeCut(MinSize);
            _clustersToFile(Difference, Clusters, MinSize, out);

            /*
             NetworkA = Operations.calculateTOM(NetworkA, threads);
             NetworkB = Operations.calculateTOM(NetworkB, threads);
             O2 = out + "/TOMA.actual.jpeg";
             Operations.generateHistogramHM(NetworkA, O2, "Topological Overlaps Network A", "TOM", "Count", false);
             O2 = out + "/TOMB.actual.jpeg";
             Operations.generateHistogramHM(NetworkB, O2, "Topological Overlaps Network B", "TOM", "Count", false);
             Difference = Operations.calculateDifference(NetworkA, NetworkB);
             String O1 = out + "/Pairwise.actual.jpeg";
             Operations.generateHistogramHM(Difference, O1, "Pairwise Edge-based TOM Differences Zm vs Sv", "cross-pair Delta-TOM", "Count", false);
             //            Operations.generateHistogramHM(CurrentMatrix, ThisOut, "Masked Distribution of Topological Overlaps", "Topological Overlap", "# Edges", false);
             System.exit(0);
             System.exit(0);
             */
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

        OptionBuilder.withArgName("alpha");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("alpha parameter for sigmoid adjacency calculation (i.e. 20)");
        Option alpha = OptionBuilder.create("a");

        OptionBuilder.withArgName("mu");
        OptionBuilder.hasArg();
        OptionBuilder.withDescription("mu parameter for sigmoid adjacency calculation (i.e. 0.8)");
        Option mu = OptionBuilder.create("m");

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
        options.addOption(output);
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

    private static void _clustersToFile(GCNMatrix Similarity, ArrayList<int[]> Clusters, int MinSize, String OutDir) {
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
            String ClustDir = OutDir + "/Clusters/";
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
