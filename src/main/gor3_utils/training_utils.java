package main.gor3_utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class training_utils {

    // alphabetical to write to the count matrices later
    public static final String AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";

    /*
     * Write the GOR matrices with the raw counts to a file
     */
    public static int write_matrices_to_file(int[][][][] gorMatrices, String filename) {
        try {
            Path outputPath = Paths.get(filename);
            StringBuilder output = new StringBuilder();
            char[] states = {'C', 'E', 'H'};

            // build output
            output.append("// Matrix4D\n\n");

            for (int centerAA = 0; centerAA < 20; centerAA++) { // each aa has 3 matrices in model file
                for (int stateIndex = 0; stateIndex < 3; stateIndex++) { // each state for any given aa
                    output.append("=").append(AA_ALPHABET.charAt(centerAA)).append(",").append(states[stateIndex]).append("=\n\t\n"); // matrix header
                    for (int offsetAA = 0; offsetAA < 20; offsetAA++) { // the rows of the matrix
                        output.append(AA_ALPHABET.charAt(offsetAA)).append("\t");
                        for (int column = 0; column < 17; column++) {
                            output.append(gorMatrices[stateIndex][centerAA][offsetAA][column]).append("\t");
                        }

                        output.append("\n");
                    }

                    output.append("\n\n");
                }
            }

            Files.writeString(outputPath, output.toString());

            return 0;
        } catch (IOException e) {
            return 1;
        }
    }

    /*
     * get raw counts of amino acids at each position in each window for each
     * state, gives a 3d array of all of the raw counts that can be used later for
     * information content
     */
    public static int[][][][] compute_raw_counts(List<String> sequences, List<String> labels) {
        // 3: states, 20: center AA, 20: offset AA, 17: window positions
        int[][][][] counts = new int[3][20][20][17];

        char[] states = {'C', 'E', 'H'};

        for (int listIndex = 0; listIndex < sequences.size(); listIndex++) {
            String currentSeq = sequences.get(listIndex);
            String currentLabels = labels.get(listIndex);

            int seqLength = currentSeq.length();

            // because of the windows we need to look at windows that contain 8 res in the seq before and after the center res
            for (int i = 8; i < seqLength - 8; i++) {

                char centerAA = currentSeq.charAt(i);
                int centerAAIndex = AA_ALPHABET.indexOf(centerAA);

                if (centerAAIndex == -1) {
                    continue;
                }

                char currentLabel = currentLabels.charAt(i);

                // find the index of the current center label in the states list
                int stateIndex = -1;
                for (int s = 0; s < states.length; s++) {
                    if (currentLabel == states[s]) {
                        stateIndex = s;
                        break;
                    }
                }

                if (stateIndex == -1) { // not valid state
                    continue;
                }

                // window
                for (int m = -8; m <= 8; m++) {
                    int residueIndex = i + m;

                    // check amino acid at position i + m
                    char offsetAA = currentSeq.charAt(residueIndex);
                    int offsetAAIndex = AA_ALPHABET.indexOf(offsetAA);
                    if (offsetAAIndex == -1) {
                        continue;
                    }

                    // get correct index for count matrix and increment
                    int columnIndex = m + 8;
                    counts[stateIndex][centerAAIndex][offsetAAIndex][columnIndex]++;
                }
            }
        }

        return counts;
    }

}