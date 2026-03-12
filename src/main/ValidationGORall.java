
package main;

import static main.utils.general_utils.findScript;
import static main.utils.parse_utils.parse_prd;
import static main.utils.parse_utils.parse_seclib;
import static main.validation_utils.val_utils.write_detailed_file;
import static main.validation_utils.val_utils.write_summary_file;
import static main.validation_utils.val_utils.write_summary_file_from_detailed_file;
import main.utils.Pair;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;

public class ValidationGORall{

    public static void main(String[] args) {
        try {
            String predictions = null;
            String seclib_file = null;
            String summary_file = null;
            String detailed_file = null;
            String format = null;

            for (int i = 0; i < args.length; i++) {
                switch (args[i]) {
                    case "-p":
                        if (i + 1 >= args.length) stop("-p requires a value");
                        predictions = args[++i];
                        break;
                    case "-r":
                        if (i + 1 >= args.length) stop("-r requires a value");
                        seclib_file = args[++i];
                        break;
                    case "-s":
                        if (i + 1 >= args.length) stop("-s requires a value");
                        summary_file = args[++i];
                        break;
                    case "-d":
                        if (i + 1 >= args.length) stop("-d requires a value");
                        detailed_file = args[++i];
                        break;
                    case "-f":
                        if (i + 1 >= args.length) stop("-f requires a value");
                        format = args[++i];
                        if (!format.equals("txt") && !format.equals("html")) {
                            stop("Invalid format: " + format + ". Must be 'txt' or 'html'.");
                        }
                        break;
                    default:
                        stop("Unknown option: " + args[i]);
                }
            }

            try {

                if (predictions == null) stop("Missing required option: -p (predictions file)");
                if (seclib_file == null) stop("Missing required option: -r (seclib file)");
                if (summary_file == null) stop("Missing required option: -s (summary file)");
                if (detailed_file == null) stop("Missing required option: -d (detailed file)");
                if (format == null) stop("Missing required option: -f (format)");

                // parse the files
                HashMap<String, Pair<String, String>> preds = parse_seclib(predictions);
                HashMap<String, Pair<String, String>> idsSequencesStructures = parse_seclib(seclib_file);
                HashMap<String, Pair<String, Pair<String, String>>> idsSequencesObservedStructuresPredictedStructures = new HashMap<>();

                for (String seqId : preds.keySet()) {
                    String aaSequence = preds.get(seqId).getV1();
                    String predictedSequence = preds.get(seqId).getV2();
                    String observedSequence = idsSequencesStructures.get(seqId).getV2();

                    idsSequencesObservedStructuresPredictedStructures.put(seqId, new Pair<String, Pair<String, String>>(aaSequence, new Pair<String, String>(observedSequence, predictedSequence)));
                }

                ArrayList<String> ids = new ArrayList<>();
                ArrayList<String> predictedLabels = new ArrayList<>();
                ArrayList<String> trueLabels = new ArrayList<>();
                ArrayList<String> trueSequences = new ArrayList<>();

                for (String seqId : idsSequencesObservedStructuresPredictedStructures.keySet()) {
                    ids.add(seqId);
                    predictedLabels.add(idsSequencesObservedStructuresPredictedStructures.get(seqId).getV2().getV2());
                    trueLabels.add(idsSequencesObservedStructuresPredictedStructures.get(seqId).getV2().getV1());
                    trueSequences.add(idsSequencesObservedStructuresPredictedStructures.get(seqId).getV1());
                }

                // compute metrics
                // write_summary_file(summary_file, ids, predictedLabels, trueLabels);
                write_detailed_file(detailed_file, ids, trueSequences, predictedLabels, trueLabels);
                write_summary_file_from_detailed_file(summary_file, ids, predictedLabels, trueLabels, detailed_file);


            } catch (Exception e) {
                stop("Error: " + e.getMessage());
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
        System.err.println("Usage: java ValidationGORall -p <predictions file path> -r <seclib file path> -s <summary file path> -d <detailed file path> -f <txt|html>");
        System.exit(1);
    }

}


