package main.gor4_utils;

import java.util.List;

  // int[][][][][][] matrix6D = new int[3][20][20][17][20][17]; // [state][AA_p1][AA_center][p1_pos][AA_p2][p2_pos]
  // int[][][][] matrix4D    = new int[20][3][20][17];           // [center_AA][state][window_AA][window_pos]


public class training_utils {

    // alphabetical to write to the count matrices later
    public static final String AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY";

    /*
     * Write the 4d and 6d matrices with the raw counts to a file
     */
    public static int write_matrices_to_file(int[][][][][][] matrix6D, int[][][][] matrix4D, String filename) {
        try {
            java.io.BufferedWriter writer = new java.io.BufferedWriter(new java.io.FileWriter(filename));
            char[] states = { 'H', 'E', 'C' };

            writer.write("// Matrix6D\n");
            // for all states
            for (int s = 2; s >= 0; s--) {
                // the center amino acid
                for (int centerAA = 0; centerAA < 20; centerAA++) {
                    // the amino acids at position 1 of the pair
                    for (int aa_p1 = 0; aa_p1 < 20; aa_p1++) {
                        // the position of the first amino acid in the window
                        for (int p1_pos = 0; p1_pos < 17; p1_pos++) {
                            // write the header for state, amino acid at p1, center aa, window position of p1 (-8 beacuse of the 0 to 17 indexing)
                            // for loops in this order because of the order of the header
                            writer.write("\n=" + states[s] + "," + AA_ALPHABET.charAt(centerAA) + "," + AA_ALPHABET.charAt(aa_p1) + "," + (p1_pos - 8) + "=" + "\n\n");
                            for (int aa_p2 = 0; aa_p2 < 20; aa_p2++) {
                                writer.write(AA_ALPHABET.charAt(aa_p2) + "\t");
                                for (int p2_pos = 0; p2_pos < 17; p2_pos++) {
                                    int count = matrix6D[s][aa_p1][centerAA][p1_pos][aa_p2][p2_pos];
                                    writer.write(count + "\t");
                                }
                                writer.write("\n");
                            }
                        }
                    }
                }
            }

            writer.write("\n// Matrix4D\n\n");
            for (int centerAA = 0; centerAA < 20; centerAA++) {
                for (int state = 2; state >= 0; state--) {
                    writer.write("=" + AA_ALPHABET.charAt(centerAA) + "," + states[state] + "=" + "\n\t");
                    for (int windowAA = 0; windowAA < 20; windowAA++) {
                        writer.write("\n" + AA_ALPHABET.charAt(windowAA) + "\t");
                        for (int window_pos = 0; window_pos < 17; window_pos++) {
                            int count = matrix4D[centerAA][state][windowAA][window_pos];
                            writer.write(count + "\t");
                        }
                    }
                    writer.write("\n\n\n");
                }
            }
            writer.close();
            return 0;
        } catch (java.io.IOException e) {
            e.printStackTrace();
            return 1;
        }
    }

    public static int[][][][] compute_4d_counts(List<String> sequences, List<String> labels) {

        int[][][][] counts = new int[20][3][20][17]; //[centerAA_idx][state_idx][windowAA_idx][window_pos]

        char[] states = { 'H', 'E', 'C' };

        // go thotough all seqs
        for (int listIndex = 0; listIndex < sequences.size(); listIndex++) {

            String currentSeq = sequences.get(listIndex);
            String currentLabels = labels.get(listIndex);

            int seqLength = currentSeq.length();

            // because of the windows we need to look at windows that contain 8 res in the seq before and after the center res
            for (int i = 8; i < seqLength - 8; i++) {

                char currentLabel = currentLabels.charAt(i);

                // find the index of the current center label in the states list (0, 1, 2)
                int stateIndex = -1;
                for (int s = 0; s < states.length; s++) {
                    if (currentLabel == states[s]) {
                        stateIndex = s;
                        break;
                    }
                }
                if (stateIndex == -1) continue; // not valid state

                char centerAA = currentSeq.charAt(i);
                int centerAAIndex = AA_ALPHABET.indexOf(centerAA);
                if (centerAAIndex == -1) continue;

                // window
                for (int m = -8; m <= 8; m++) {
                    int residueIndex = i + m;
                    // needs to be in the sequence
                    if (residueIndex < 0 || residueIndex >= seqLength) {
                        continue;
                    }

                    // check which amino acid is at the position i + m
                    char aa = currentSeq.charAt(residueIndex);
                    int aaIndex = AA_ALPHABET.indexOf(aa);
                    if (aaIndex == -1) continue;

                    // get the correct index for the count matrix and increment the position we found
                    int window_pos_idx = m + 8;
                    counts[centerAAIndex][stateIndex][aaIndex][window_pos_idx]++;
                }
            }
        }
        return counts;
    }

    public static int[][][][][][] compute_6d_counts(List<String> sequences, List<String> labels) {

        int[][][][][][] counts = new int[3][20][20][17][20][17]; // [state][AA_p1][AA_center][p1_window_idx][AA_p2][p2_window_idx]

        char[] states = { 'H', 'E', 'C' };

        // go thotough all seqs
        for (int listIndex = 0; listIndex < sequences.size(); listIndex++) {

            String currentSeq = sequences.get(listIndex);
            String currentLabels = labels.get(listIndex);

            int seqLength = currentSeq.length();

            // because of the windows we need to look at windows that contain 8 res in the seq before and after the center res
            for (int i = 8; i < seqLength - 8; i++) {

                char currentLabel = currentLabels.charAt(i);

                // find the index of the current center label in the states list (0, 1, 2)
                int stateIndex = -1;
                for (int s = 0; s < states.length; s++) {
                    if (currentLabel == states[s]) {
                        stateIndex = s;
                        break;
                    }
                }
                if (stateIndex == -1) continue; // not valid state

                char centerAA = currentSeq.charAt(i);
                int centerAAIndex = AA_ALPHABET.indexOf(centerAA);
                if (centerAAIndex == -1) continue;

                // want to look at all pais of amino acids in the window around the center amino acid
                // for the first amino acid loop through the window using m1
                for (int m1 = -8; m1 <= 8; m1++) {
                    int residueIndex1 = i + m1;
                    // needs to be in the sequence
                    if (residueIndex1 < 0 || residueIndex1 >= seqLength) {
                        continue;
                    }

                    char aa_p1 = currentSeq.charAt(residueIndex1);
                    int aa_p1_Index = AA_ALPHABET.indexOf(aa_p1);
                    if (aa_p1_Index == -1) continue;

                    // for the second amino acid we need to loop through the window using m2
                    // but m2 needs to start after m1 to avoid counting same pair twice and also avoid same position pairs
                    for (int m2 = m1 + 1; m2 <= 8; m2++) {
                        int residueIndex2 = i + m2;
                        // needs to be in the sequence
                        if (residueIndex2 < 0 || residueIndex2 >= seqLength) {
                            continue;
                        }

                        char aa_p2 = currentSeq.charAt(residueIndex2);
                        int aa_p2_Index = AA_ALPHABET.indexOf(aa_p2);
                        if (aa_p2_Index == -1) continue;

                        // get the correct index for the count matrix and increment the position we found
                        int windowIndex1 = m1 + 8;
                        int windowIndex2 = m2 + 8;

                        counts[stateIndex][aa_p1_Index][centerAAIndex][windowIndex1][aa_p2_Index][windowIndex2]++;
                    }
                }
            }
        }
        return counts;
    }

}
