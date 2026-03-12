package main.gor1_utils;

import java.util.List;

public class training_utils {

    // alphabetical to write to the count matrices later
    public static final String AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";

    /*
     * Write the GOR matrices with the raw counts to a file
     */
    public static int write_matrices_to_file(int[][][] gorMatrices, String filename) {
        try {
            java.io.BufferedWriter writer = new java.io.BufferedWriter(new java.io.FileWriter(filename));

            writer.write("// Matrix3D\n\n");

            char[] states = { 'H', 'E', 'C' };

            // get states TODO: this should be made better
            for (int s = 2; s >= 0; s--) {

                String header = "=" + states[s] + "=\n\t\n";
                writer.write(header);

                // get amino acids
                for (int a = 0; a < AA_ALPHABET.length(); a++) {

                    StringBuilder line = new StringBuilder();

                    // need the amino acid and the values after
                    line.append(AA_ALPHABET.charAt(a) + "\t");
                    for (int m = 0; m < 17; m++) {
                        line.append(gorMatrices[s][a][m] + "\t");
                    }
                    writer.write(line.toString());
                    writer.write("\n");
                }
                writer.write("\n\n");
            }
            writer.close();
            return 0;
        } catch (java.io.IOException e) {
            e.printStackTrace();
            return 1;
        }
    }

    /*
     * get raw counts of amino acids at each position in each window for each
     * state, gives a 3d array of all of the raw counts that can be used later for
     * information content
     */
    public static int[][][] compute_raw_counts(List<String> sequences, List<String> labels) {
        int[][][] counts = new int[3][20][17];

        char[] states = { 'H', 'E', 'C' };

        // go thotough all seqs
        for (int listIndex = 0; listIndex < sequences.size(); listIndex++) {

            String currentSeq = sequences.get(listIndex);
            String currentLabels = labels.get(listIndex);

            int seqLength = currentSeq.length();

            // because of the windows we need to look at windows that contain 8 res in the seq before and after the center res
            for (int i = 8; i < seqLength - 8; i++) {

                char currentLabel = currentLabels.charAt(i);

                // find the index of the current center label in the states list (0, 1 or 2)
                int stateIndex = -1;
                for (int s = 0; s < states.length; s++) {
                    if (currentLabel == states[s]) {
                        stateIndex = s;
                        break;
                    }
                }
                if (stateIndex == -1) continue; // not valid state

                // window
                for (int m = -8; m <= 8; m++) {
                    int residueIndex = i + m;

                    // check which amino acid is at the position i + m
                    char aa = currentSeq.charAt(residueIndex);
                    int aaIndex = AA_ALPHABET.indexOf(aa);
                    if (aaIndex == -1) continue;

                    // get the correct index for the count matrix and increment the position we found
                    int columnIndex = m + 8;
                    counts[stateIndex][aaIndex][columnIndex]++;
                }
            }
        }
        return counts;
    }

}
