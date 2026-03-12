package main;

import java.io.*;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;

import main.utils.Constants;
import main.utils.Pair;
import main.utils.post_processing_utils;
import static main.utils.parse_utils.*;
import static main.utils.output_utils.*;
import static main.utils.general_utils.*;
import static main.utils.Constants.*;

public class PredictionGORall {

    public static void main(String[] args) throws IOException {
        try {
            boolean outputProbabilities = false;
            boolean coloredOutput = false;
            boolean usePSa = false;
            boolean postProcessing = false;
            String model = null;
            String format = null;
            String seqFile = null;
            String mafFolder = null;

            for (int i = 0; i < args.length; i++) {
                switch (args[i]) {
                    case "--probabilities":
                        outputProbabilities = true;
                        break;
                    case "--colored":
                        coloredOutput = true;
                        break;
                    case "--usepsa":
                        usePSa = true;
                        break;
                    case "--postprocessing":
                        postProcessing = true;
                        break;
                    case "--model":
                        if (i + 1 >= args.length) stop("--model requires a value");
                        model = args[++i];
                        break;
                    case "--format":
                        if (i + 1 >= args.length) stop("--format requires a value");
                        format = args[++i];
                        if (!format.equals("txt") && !format.equals("html")) {
                            stop("Invalid format: " + format + ". Must be 'txt' or 'html'.");
                        }
                        break;
                    case "--seq":
                        if (i + 1 >= args.length) stop("--seq requires a value");
                        seqFile = args[++i];
                        break;
                    case "--maf":
                        if (i + 1 >= args.length) stop("--maf requires a value");
                        mafFolder = args[++i];
                        break;
                    default:
                        stop("Unknown option: " + args[i]);
                }
            }

            if (model == null) stop("Missing required option: --model");
            if (format == null) stop("Missing required option: --format");
            if (seqFile == null && mafFolder == null) stop("Either --seq or --maf must be provided");
            if (seqFile != null && mafFolder != null) stop("Cannot provide both --seq and --maf");

            // output regardless of prediction method
            HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities = new HashMap<>(); // hash map of seqId: (sequence, probabilities); probabilities is [length][3]

            if (seqFile != null) {
                int method = identifyModelFile(model);

                // parse input fasta
                HashMap<String, String> idsSequences = main.utils.parse_utils.parse_fasta(seqFile); // map of seq_id: sequence

                if (method == 1) {
                    idsSequencesProbabilities = predictGOR1(model, idsSequences, usePSa);
                } else if (method == 3) {
                    idsSequencesProbabilities = predictGOR3(model, idsSequences, usePSa);
                } else if (method == 4) {
                    idsSequencesProbabilities = predictGOR4(model, idsSequences);
                }
            } else {
                idsSequencesProbabilities = predictGOR5(model, mafFolder);
            }

            String output = buildOutput(idsSequencesProbabilities, outputProbabilities, coloredOutput, format, postProcessing, (rawPrediction, windowSize) -> post_processing_utils.applySmoothing(rawPrediction, windowSize));
            System.out.println(output);
        } catch (Exception e) {
            String exceptionType = e.getClass().getName();
            String exceptionMessage = e.getMessage();

            StringWriter sw = new StringWriter();
            e.printStackTrace(new PrintWriter(sw));
            String stackTrace = sw.toString();

            String cause = e.getCause() != null ? e.getCause().toString() : "none";

            String message = "Exception: " + exceptionType + "\nMessage: " + exceptionMessage + "\nCause: " + cause + "\nStacktrace: " + stackTrace;

            try {
                System.out.println("Trying to find script...");
                Path script = findScript();
                System.out.println("Found script: " + script);
                ProcessBuilder pb = new ProcessBuilder("python3", "-u", script.toString(), "--message", message);
                pb.redirectErrorStream(true);
                Process process = pb.start();

                // output ticket creation process to stdout
                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println("[ticket script] " + line);
                }
                process.waitFor();
            } catch (Exception scriptEx) {
                scriptEx.printStackTrace();
            }
        }
    }

    public static HashMap<String, Pair<String, double[][]>> predictGOR1(String model, HashMap<String, String> idsSequences, boolean usePSa) {
        // load model from file
        int[][][] rawCountMatrices = main.gor1_utils.prediction_utils.read_matrices_from_file(model);
        double[][][] scoreMatrices = main.gor1_utils.prediction_utils.counts_to_scores(rawCountMatrices, usePSa);

        HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities = new HashMap<>();

        // calculate probabilities for each sequence
        for (String seqId : idsSequences.keySet()) {
            double[][] probabilities;
            if (idsSequences.get(seqId).length() >= 17) { // could be ternary, but then the line is super long
                probabilities = main.gor1_utils.prediction_utils.predictSequence(idsSequences.get(seqId), scoreMatrices);
            } else {
                probabilities = null;
            }
            idsSequencesProbabilities.put(seqId, new Pair<>(idsSequences.get(seqId), probabilities));
        }

        return idsSequencesProbabilities;
    }

    public static HashMap<String, Pair<String, double[][]>> predictGOR3(String model, HashMap<String, String> idsSequences, boolean usePSa) {
        // load model from file
        int[][][][] rawCountMatrices = main.gor3_utils.prediction_utils.read_matrices_from_file(model);
        double[][][][] scoreMatrices = main.gor3_utils.prediction_utils.counts_to_scores(rawCountMatrices, usePSa);

        HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities = new HashMap<>();

        // calculate probabilities for each sequence
        for (String seqId : idsSequences.keySet()) {
            double[][] probabilities;
            if (idsSequences.get(seqId).length() >= 17) { // could be ternary, but then the line is super long
                probabilities = main.gor3_utils.prediction_utils.predictSequence(idsSequences.get(seqId), scoreMatrices);
            } else {
                probabilities = null;
            }
            idsSequencesProbabilities.put(seqId, new Pair<>(idsSequences.get(seqId), probabilities));
        }

        return idsSequencesProbabilities;
    }

    public static HashMap<String, Pair<String, double[][]>> predictGOR4(String model, HashMap<String, String> idsSequences) {
        // load raw counts from file
        ArrayList<Object> rawCountMatrices = main.gor4_utils.prediction_utils.read_matrices_from_file_6d_4d(model);
        if (rawCountMatrices == null || rawCountMatrices.size() != 2) {
            stop("Failed to read model file: " + model);
        }
        int [][][][][][] rawCountMatrix6D = (int[][][][][][]) rawCountMatrices.get(0);
        int [][][][] rawCountMatrix4D = (int[][][][]) rawCountMatrices.get(1);

        double[][][][][][] scoreMatrix6D = main.gor4_utils.prediction_utils.counts_to_scores_6d(rawCountMatrix6D);
        double[][][][] scoreMatrix4D = main.gor4_utils.prediction_utils.counts_to_scores_4d(rawCountMatrix4D);
        Pair<double[][][][][][], double[][][][]> scoreMatrices = new Pair<>(scoreMatrix6D, scoreMatrix4D);

        HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities = new HashMap<>();

        // calculate probabilities for each sequence
        for (String seqId : idsSequences.keySet()) {
            double[][] probabilities;
            if (idsSequences.get(seqId).length() >= 17) { // could be ternary, but then the line is super long
                probabilities = main.gor4_utils.prediction_utils.predictSequence(idsSequences.get(seqId), scoreMatrices);
            } else {
                probabilities = null;
            }
            idsSequencesProbabilities.put(seqId, new Pair<>(idsSequences.get(seqId), probabilities));
        }

        return idsSequencesProbabilities;
    }

    public static HashMap<String, Pair<String, double[][]>> predictGOR5(String model, String mafDirectory) {
        // identify model file and parse accordingly
        int modelType = identifyModelFile(model);
        Object scoringMatrices = null;

        // get correct model with generics
        if (modelType == 1) {
            int[][][] rawCountMatrices = main.gor1_utils.prediction_utils.read_matrices_from_file(model);
            scoringMatrices = main.gor1_utils.prediction_utils.counts_to_scores(rawCountMatrices);
        } else if (modelType == 3) {
            int[][][][] rawCountMatrices = main.gor3_utils.prediction_utils.read_matrices_from_file(model);
            scoringMatrices = main.gor3_utils.prediction_utils.counts_to_scores(rawCountMatrices);
        } else if (modelType == 4) {
            ArrayList<Object> raw = main.gor4_utils.prediction_utils.read_matrices_from_file_6d_4d(model);
            int [][][][][][] rawCountMatrix6D = (int[][][][][][]) raw.get(0);
            int [][][][] rawCountMatrix4D = (int[][][][]) raw.get(1);
            double[][][][][][] scoreMatrix6D = main.gor4_utils.prediction_utils.counts_to_scores_6d(rawCountMatrix6D);
            double[][][][] scoreMatrix4D = main.gor4_utils.prediction_utils.counts_to_scores_4d(rawCountMatrix4D);
            scoringMatrices = new main.utils.Pair<>(scoreMatrix6D, scoreMatrix4D);

        } else {
            stop("This is not a valid model file!");
        }

        // final output hashmap
        HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities = new HashMap<>();

        // get all alignment files in folder
        ArrayList<String> filePaths = getAlignmentFiles(mafDirectory);
        for (String filePath : filePaths) {
            // parse alignment file
            Pair<ArrayList<String>, Pair<String, String>> alignmentFileContent = pastaAln(filePath);
            ArrayList<String> alignmentSequences = alignmentFileContent.getV1();
            String seqId = alignmentFileContent.getV2().getV1();
            String originalSequence = alignmentFileContent.getV2().getV2();

            // predict sequence probabilities
            double[][] averageProbabilities;
            if (originalSequence.length() >= 17) {
                Double[][][] allProbabilities = main.gor5_utils.prediction_utils.generateProbabilitiesAllSequences(
                        alignmentSequences,
                        (window, matrices) -> {
                            if (matrices instanceof double[][][]) {
                                return main.gor1_utils.prediction_utils.predictSingleResidue(window, (double[][][]) matrices);
                            } else if (matrices instanceof double[][][][]) {
                                return main.gor3_utils.prediction_utils.predictSingleResidue(window, (double[][][][]) matrices);
                            } else if (matrices instanceof Pair) {
                                @SuppressWarnings("unchecked")
                                main.utils.Pair<double[][][][][][], double[][][][]> pairMatrices = (Pair<double[][][][][][], double[][][][]>) matrices;
                                return main.gor4_utils.prediction_utils.predictSingleResidue(window, pairMatrices);
                            } else {
                                return null;
                            }
                        },
                        scoringMatrices);
                // average probabilities across all alignment sequences
                averageProbabilities = main.gor5_utils.prediction_utils.averageProbabilities(alignmentSequences, allProbabilities);
            } else {
                averageProbabilities = null;
            }

            idsSequencesProbabilities.put(seqId, new Pair<>(originalSequence, averageProbabilities));
        }

        return idsSequencesProbabilities;
    }
}


