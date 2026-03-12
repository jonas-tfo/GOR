
package main;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Objects;
import java.util.HashMap;

import static main.utils.general_utils.findScript;
import static main.utils.parse_utils.parse_seclib;
import main.utils.Pair;

public class TrainingGORall {

    public static void main(String[] args) {
        try {
            String db = null; // training file
            String method = null; // gor method
            String model = null; // output file for model (matrices with the raw counts)

            for (int i = 0; i < args.length; i++) {
                switch (args[i]) {
                    case "--db":
                        if (i + 1 >= args.length) stop("--db requires a value");
                        db = args[++i];
                        break;
                    case "--method":
                        if (i + 1 >= args.length) stop("--method requires a value");
                        method = args[++i];
                        break;
                    case "--model":
                        if (i + 1 >= args.length) stop("--model requires a value");
                        model = args[++i];
                        break;
                    default:
                        stop("Unknown option: " + args[i]);
                }
            }

            if (db == null) stop("Missing required option: --db");
            if (method == null) stop("Missing required option: --method");
            if (model == null) stop("Missing required option: --model");

            // parse data
            //ArrayList<ArrayList<String>> parsedData = parse_seclib(db);
            HashMap<String, Pair<String, String>> idsSequencesStructures = parse_seclib(db);

            if (idsSequencesStructures == null) stop("Failed to parse seclib file: " + db);
            // get the data from parser output

            ArrayList<String> aminoAcidSequences = new ArrayList<>();
            ArrayList<String> labelSequences = new ArrayList<>();

            for (String seqId : idsSequencesStructures.keySet()) {
                aminoAcidSequences.add(idsSequencesStructures.get(seqId).getV1());
                labelSequences.add(idsSequencesStructures.get(seqId).getV2());
            }

            // get raw counts in 3d matrix for gor1
            if (Objects.equals(method, "gor1")) {

                // compute raw counts for 3 states
                int[][][] rawCounts = main.gor1_utils.training_utils.compute_raw_counts(aminoAcidSequences, labelSequences);

                int exit = main.gor1_utils.training_utils.write_matrices_to_file(rawCounts, model);
                if (exit == 1) stop("Failed to write model to file: " + model + " check if the path is correct");

                // get raw counts in 4d matrix for gor3
            } else if (Objects.equals(method, "gor3")) {

                // compute raw counts for 3 states
                int[][][][] rawCounts = main.gor3_utils.training_utils.compute_raw_counts(aminoAcidSequences, labelSequences);

                int exit = main.gor3_utils.training_utils.write_matrices_to_file(rawCounts, model);
                if (exit == 1) stop("Failed to write model to file: " + model + " check if the path is correct.");

                // get raw counts in 6d and 4d matrices for gor4
            } else if (Objects.equals(method, "gor4")) {

                // compute raw counts for 3 states
                int[][][][][][] rawCounts6D = main.gor4_utils.training_utils.compute_6d_counts(aminoAcidSequences, labelSequences);
                int[][][][] rawCounts4D = main.gor4_utils.training_utils.compute_4d_counts(aminoAcidSequences, labelSequences); // TODO: merge this and gor3 functions

                int exit = main.gor4_utils.training_utils.write_matrices_to_file(rawCounts6D, rawCounts4D, model);
                if (exit == 1) stop("Failed to write model to file: " + model + " check if the path is correct");

                // TODO: prediction

            }
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

    // helper to stop the program with error and usage stuff
    private static void stop(String message) {
        System.err.println("Error: " + message);
        System.err.println("Usage: java -jar train.jar --db <seclib-file> --method <gor1|gor3|gor4> --model <model-file>");
        System.exit(1);
    }
}


