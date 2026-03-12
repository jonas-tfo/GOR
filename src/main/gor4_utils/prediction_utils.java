package main.gor4_utils;

import java.util.ArrayList;
import main.utils.Pair;
import static main.utils.Constants.*;

public class prediction_utils {
    // save init scores for the states to add to scores before prediction
    public static double[] GLOBAL_PRIORS = new double[3];

    /*
     * Reads the 4d and 6d matrices with the raw counts from training from a file and returns both
     */
    public static ArrayList<Object> read_matrices_from_file_6d_4d(String filename) {

        int[][][][][][] matrix6D = new int[3][20][20][17][20][17]; // [state][AA_p1][AA_center][p1_pos][AA_p2][p2_pos]
        int[][][][] matrix4D = new int[20][3][20][17];           // [center_AA][state][window_AA][window_pos]

        try {
            java.io.BufferedReader reader = new java.io.BufferedReader(new java.io.FileReader(filename));

            String line;
            boolean reading6D = false;
            boolean reading4D = false;

            // init for chekcing if were set later
            int stateIndex = -1;
            int centerAAIndex = -1;
            int aa_p1_Index = -1;
            int p1_window_pos = -1;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("// Matrix6D")) {
                    reading6D = true;
                    reading4D = false;
                    continue;
                } else if (line.startsWith("// Matrix4D")) {
                    reading6D = false;
                    reading4D = true;
                    continue;
                }

                if (reading6D) {

                    // parse 6D matrix lines
                    if (line.startsWith("=") && line.endsWith("=")) {
                        String[] headerParts = line.split(",");
                        stateIndex = STATES.indexOf(headerParts[0].charAt(1)); // gets the state, after the =
                        centerAAIndex = AA_ALPHABET.indexOf(headerParts[1].charAt(0));
                        aa_p1_Index = AA_ALPHABET.indexOf(headerParts[2].charAt(0));
                        p1_window_pos = Integer.parseInt(headerParts[3].substring(0, headerParts[3].length() - 1)) + 8; // convert to normal indexing
                        if (stateIndex == -1 || centerAAIndex == -1 || aa_p1_Index == -1 || p1_window_pos < 0 || p1_window_pos >= 17) {
                            continue; // invalid header
                        }
                    } else {
                        String[] parts = line.split("\t");
                        if (parts.length < 18) continue; // need 1 for amino acid and 17 for counts + 1 empty at the end

                        // get aa from beginning of line
                        int aa_p2_Index = AA_ALPHABET.indexOf(parts[0].charAt(0));
                        if (aa_p2_Index == -1) continue; // invalid amino acid

                        // go through the counts in the window for the current amino acid (p2) and save with window position
                        for (int p2_window_pos = 0; p2_window_pos < 17; p2_window_pos++) {
                            // get each count, skip the aa letter at pos 0
                            int count = Integer.parseInt(parts[p2_window_pos + 1]);
                            matrix6D[stateIndex][aa_p1_Index][centerAAIndex][p1_window_pos][aa_p2_Index][p2_window_pos] = count;
                        }
                    }
                } else if (reading4D) {

                    // parse 4D matrix lines
                    if (line.startsWith("=") && line.endsWith("=")) {
                        String[] headerParts = line.split(",");

                        centerAAIndex = AA_ALPHABET.indexOf(headerParts[0].charAt(1));
                        stateIndex = STATES.indexOf(headerParts[1].charAt(0));

                        if (centerAAIndex == -1 || stateIndex == -1) {
                            continue; // invalid header
                        }
                    } else {
                        String[] parts = line.split("\t");
                        if (parts.length < 18) continue; // need 1 for amino acid and 17 for counts + 1 empty at the end

                        int windowAAIndex = AA_ALPHABET.indexOf(parts[0].charAt(0));
                        if (windowAAIndex == -1) continue; // invalid amino acid

                        // go through the counts for the current amino acid and save with window position
                        for (int window_pos = 0; window_pos < 17; window_pos++) {
                            // get each count, skip the aa letter at pos 0
                            int count = Integer.parseInt(parts[window_pos + 1]);
                            // this needs to be the same shape of whats used for the count to score in the pred
                            matrix4D[centerAAIndex][stateIndex][windowAAIndex][window_pos] = count;
                        }
                    }
                }
            }

            reader.close();

        } catch (java.io.IOException e) {
            e.printStackTrace();
            return null;
        }

        ArrayList<Object> result = new ArrayList<>();
        result.add(matrix6D);
        result.add(matrix4D);
        return result;
    }

    public static double calculate_information_content_score_gor3(double fSAmaa, double fNotSAmaa, double fSaa, double fNotSaa) {
        // information content score formula
        // I(S, cA, oA, m) = log( f(S, cA, oA, m) / f(not S, cA, oA, m) ) + log( f(not S, cA) / f(S, cA) )

        // pseudocounts, avoid division by zero and log of zero
        fSAmaa += GOR4_PSEUDOCOUNT;
        fNotSAmaa += GOR4_PSEUDOCOUNT;
        fSaa += GOR4_PSEUDOCOUNT;
        fNotSaa += GOR4_PSEUDOCOUNT;

       return Math.log(fSAmaa / fNotSAmaa) + Math.log(fNotSaa / fSaa);
    }

    public static double calculate_information_content_score_gor4(int fSPair, int fNotSPair, long fS, long fNotS) {

        double fSPairPseudo = fSPair + GOR4_PSEUDOCOUNT;
        double fNotSPairPseudo = fNotSPair + GOR4_PSEUDOCOUNT;
        double fNotSPseudo = fNotS + GOR4_PSEUDOCOUNT;
        double fSPseudo = fS + GOR4_PSEUDOCOUNT;

        return Math.log(fSPairPseudo / fNotSPairPseudo) + Math.log(fNotSPseudo / fSPseudo);
    }

    public static double[][][][] counts_to_scores_4d(int[][][][] counts) {

        double[][][][] scores = new double[20][3][20][17];

        long[] globalFS = new long[3];
        long[] globalFNotS = new long[3];

        // calculate global state frequencies
        // look at 8 only, prevent overcounting
        for (int centerAA = 0; centerAA < 20; centerAA++) {
            for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
                for (int offsetAA = 0; offsetAA < 20; offsetAA++) {
                    globalFS[stateIndex] += counts[centerAA][stateIndex][offsetAA][8];
                }
            }
        }

        //  set the global prior for prediction later for all of the states
        for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
            globalFNotS[stateIndex] = (globalFS[0] + globalFS[1] + globalFS[2]) - globalFS[stateIndex]; // total minus globalFS to get globalFNotS
            // get global prior for each state, second part of the info content formula on slide
            GLOBAL_PRIORS[stateIndex] = Math.log((globalFS[stateIndex] + GOR4_PSEUDOCOUNT) / (globalFNotS[stateIndex] + GOR4_PSEUDOCOUNT));
        }

        for (int centerAA = 0; centerAA < 20; centerAA++) {
            for (int stateIndex = 0; stateIndex <  3; stateIndex++) {
                for (int offsetAA = 0; offsetAA < 20; offsetAA++) {
                    for (int m = 0; m < 17; m++) {
                        int fSAmaa = counts[centerAA][stateIndex][offsetAA][m];
                        int fNotSAmaa = 0;
                        for (int otherStateIndex = 0; otherStateIndex < 3; otherStateIndex++) {
                            if (otherStateIndex != stateIndex) {
                                fNotSAmaa += counts[centerAA][otherStateIndex][offsetAA][m];
                            }
                        }
                        // use global frequencies here, not the center-AA specific frequencies, that was making performance horrible
                        scores[centerAA][stateIndex][offsetAA][m] += calculate_information_content_score_gor3(
                                fSAmaa, fNotSAmaa, globalFS[stateIndex], globalFNotS[stateIndex]);
                    }
                }
            }
        }
        return scores;
    }

    public static double[][][][][][] counts_to_scores_6d(int[][][][][][] counts) {

        double[][][][][][] scores = new double[3][20][20][17][20][17];

        long[] globalFS = new long[3];
        long[] globalFNotS = new long[3];

        // calculate globsal state frequencies
        // only look at 8 and 9 in the window here to prevent overcounting like in 4d
        // doesnt matter which positions, should always be same count
        for (int S = 0; S < 3; S++) {
            for (int AA_p1 = 0; AA_p1 < 20; AA_p1++) {
                for (int AA_center = 0; AA_center < 20; AA_center++) {
                    for (int AA_p2 = 0; AA_p2 < 20; AA_p2++) {
                        globalFS[S] += counts[S][AA_p1][AA_center][8][AA_p2][9];
                    }
                }
            }
        }

        for (int S = 0; S < 3; S++) {
            globalFNotS[S] = (globalFS[0] + globalFS[1] + globalFS[2]) - globalFS[S];
            GLOBAL_PRIORS[S] = Math.log((globalFS[S] + GOR4_PSEUDOCOUNT) / (globalFNotS[S] + GOR4_PSEUDOCOUNT));
        }

        for (int S=0; S<3; S++) {
            for (int AA_p1 = 0; AA_p1 < 20; AA_p1++) {
                for (int AA_center = 0; AA_center < 20; AA_center++) {
                    for (int p1_pos = 0; p1_pos < 17; p1_pos++) {
                        for (int AA_p2 = 0; AA_p2 < 20; AA_p2++) {
                            for (int p2_pos = 0; p2_pos < 17; p2_pos++) {

                                int fSPair = counts[S][AA_p1][AA_center][p1_pos][AA_p2][p2_pos];
                                int fNotSPair = 0;
                                for (int otherS=0; otherS<3; otherS++) {
                                    if (otherS != S) {
                                        fNotSPair += counts[otherS][AA_p1][AA_center][p1_pos][AA_p2][p2_pos];
                                    }
                                }
                                scores[S][AA_p1][AA_center][p1_pos][AA_p2][p2_pos] = calculate_information_content_score_gor4(
                                        fSPair, fNotSPair, globalFS[S], globalFNotS[S]);
                            }
                        }
                    }
                }
            }
        }
        return scores;
    }


    /* gets the information content scores from the raw counts, gives a 3d array of all of the scores for each state, amino acid and position */
    public static double[][][] counts_to_scores(int[][][] counts) {

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
                    scores[S][A][m] += calculate_information_content_score_gor3(fSAm, fNotSAm, fS[S], fNotS[S]);
                }
            }
        }
        return scores;
    }

    public static double[] predictSingleResidue(String window, main.utils.Pair<double[][][][][][], double[][][][]> matrices) {

        // gives init values to the scores for different states
        // these are negative, so are subtracted from the actual scores
        double[] scores = new double[3];
        for (int i = 0; i < 3; i++) {
            scores[i] = GLOBAL_PRIORS[i];
        }

        double[][][][][][] scoringMatrix6d = (double[][][][][][]) matrices.getV1();
        double[][][][] scoringMatrix4d = (double[][][][]) matrices.getV2();

        int AA_center_idx = AA_ALPHABET.indexOf(window.charAt(8));
        if (AA_center_idx == -1) {
            return scores;
            //return new double[] { helixScore, sheetScore, coilScore }; // Fallback
        }

        for (int S = 0; S < 3; S++) {
            double summed_score = 0.0;

            for (int m1 = 0; m1 < 17; m1++) {
                char aa1 = window.charAt(m1);
                int aa1_idx = AA_ALPHABET.indexOf(aa1);
                if (aa1_idx == -1) continue;

                for (int m2 = m1 + 1; m2 < 17; m2++) {
                    char aa2 = window.charAt(m2);
                    int aa2_idx = AA_ALPHABET.indexOf(aa2);
                    if (aa2_idx == -1) continue;

                    summed_score += (2.0 / 17.0) * scoringMatrix6d[S][aa1_idx][AA_center_idx][m1][aa2_idx][m2];
                }
                summed_score -= (15.0 / 17.0) * scoringMatrix4d[AA_center_idx][S][aa1_idx][m1];
            }

            scores[S] += summed_score; // add the global prior for the state to the score
        }
        return scores;
    }

    public static double[][] predictSequence(String fullSequence, Pair<double[][][][][][], double[][][][]> scoringMatrices) {
        double[][] scores = new double[fullSequence.length() - 16][3];

        for (int i = 8; i < fullSequence.length() - 8; i++) {
            String window = fullSequence.substring(i - 8, i + 9);
            scores[i - 8] = predictSingleResidue(window, scoringMatrices);
        }

        return scores;
    }
}
