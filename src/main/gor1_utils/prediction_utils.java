package main.gor1_utils;

import static main.utils.Constants.*;

public class prediction_utils {

    /*
     * Read GOR raw count matrices from a file written by write_gor_matrices_to_file.
     */
    public static int[][][] read_matrices_from_file(String filename) {
        int[][][] gorMatrices = new int[3][20][17];

        try {
            java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(filename));

            int currentState = -1;

            String line;
            while ((line = reader.readLine()) != null) {

                line = line.trim();
                if (line.isEmpty()) continue;


                if (line.startsWith("=") && line.endsWith("=")) {
                    currentState = STATES.indexOf("" + line.charAt(1));
                }
                if (currentState == -1) continue; // no valid state header found yet

                String[] parts = line.split("\t");
                if (parts.length < 18) continue; // need 1 for amino acid and 17 for counts + 1 empty at the end

                char aa = parts[0].charAt(0);
                int aaIndex = AA_ALPHABET.indexOf(aa);
                if (aaIndex == -1) continue; // invalid amino acid

                // fill the window for the according amino acid in the matrix for the according state by going through m vals
                for (int m = 0; m < 17; m++) {
                    gorMatrices[currentState][aaIndex][m] = Integer.parseInt(parts[m + 1]);
                }
            }
            reader.close();
        } catch (java.io.IOException e) {
            e.printStackTrace();
            return null;
        }
        return gorMatrices;
    }


    public static double calculate_information_content_score(int fSAm, int fNotSAm, int fS, int fNotS) {
        // information content score formula
        // I(S, A, m) = log( f(S, A, m) / f(not S, A, m) ) − log( f(S) / f(not S) )
        // or I(S, A, m) = log( f(S, A, m) / f(not S, A, m) ) + log( f(not S) / f(S) ) -> using this

        // pseudocounts, avoid division by zero and log of zero
        fSAm += GOR1_PSEUDOCOUNT;
        fNotSAm += GOR1_PSEUDOCOUNT;
        fS += GOR1_PSEUDOCOUNT;
        fNotS += GOR1_PSEUDOCOUNT;

        return Math.log((double) fSAm / fNotSAm) + Math.log((double) fNotS / fS);
    }

    public static double calculatePSa(int fSAm, int fNotSAm) {
        return ((double)(fSAm + GOR1_PSEUDOCOUNT)) / ((double)(fSAm + fNotSAm + 2 * GOR1_PSEUDOCOUNT));
    }

    // overload for counts_to_scores to not use PSa by default
    public static double[][][] counts_to_scores(int[][][] counts) {
        return counts_to_scores(counts, false);
    }

    /* gets the information content scores from the raw counts, gives a 3d array of all of the scores for each state, amino acid and position */
    public static double[][][] counts_to_scores(int[][][] counts, boolean usePSa) {

        double[][][] scores = new double[3][20][17];

        // precompute fS and fNotS for each state, can do this because only need to do it for center position
        int [] fS = new int[3]; // how often state S is found at center
        int [] fNotS = new int[3]; // how often states other than S are found at center

        for (int S = 0; S < 3; S++) {
            for (int A = 0; A < 20; A++) {
                fS[S] += counts[S][A][8]; // count how often S is found at center residue for any A
                for (int otherS = 0; otherS < 3; otherS++) {
                    if (otherS != S) {
                        fNotS[S] += counts[otherS][A][8]; // count how often states other than S are found at center res for any A
                    }
                }
            }
        }

        // do the calc for all 3 states
        for (int S = 0; S < 3; S++) {
            // look at the amino acids (y axis on slide)
            for (int A = 0; A < 20; A++) {
                // look at position in the window (not -8 to +8 for simplicity)
                for (int m = 0; m < 17; m++) {

                    int fSAm = counts[S][A][m]; // how often A is found at offset m in state S
                    int fNotSAm = 0; // how often A is found at offset m in states other than S
                    for (int otherS = 0; otherS < 3; otherS++) {
                        if (otherS != S) {
                            fNotSAm += counts[otherS][A][m];
                        }
                    }
                    // sum up over the whole window so we get sum of the scores to save
                    // use precomputed fS and fNotS for the center position
                    scores[S][A][m] += usePSa ? calculatePSa(fSAm, fS[S]) : calculate_information_content_score(fSAm, fNotSAm, fS[S], fNotS[S]);
                }
            }
        }
        return scores;
    }

    /*
     * Predict the secondary structure for a single residue based on a given
     * window, gives a list of all of the scores for the 3 states
     */
    public static double[] predictSingleResidue(String window, double[][][] scoringMatrices) {

        double[] scores = new double[3];

        // go through the states matrices
        for (int state_idx = 0; state_idx < 3; state_idx++) {

            // go through the positions in the window, third dim of the matrix
            for (int pos_idx = 0; pos_idx < 17; pos_idx++) {

                // get center amino acid
                char aminoacid = window.charAt(pos_idx);
                int aminoacid_idx = AA_ALPHABET.indexOf(aminoacid);

                // ignore invalid
                if (aminoacid_idx == -1) continue;

                scores[state_idx] += scoringMatrices[state_idx][aminoacid_idx][pos_idx];

            }
        }
        return scores;

    }

    /*
    Calculates probabilities for a single sequence. The output double[][] array only has predictions for valid positions and has length len(fullSequence) - 16
     */
    public static double[][] predictSequence(String fullSequence, double[][][] scoringMatrices) {
        double[][] scores = new double[fullSequence.length() - 16][3];

        for (int i = 8; i < fullSequence.length() - 8; i++) {
            String window = fullSequence.substring(i - 8, i + 9);
            scores[i - 8] = predictSingleResidue(window, scoringMatrices);
        }

        return scores;
    }
}
