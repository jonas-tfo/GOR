package main.gor3_utils;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import main.utils.Pair;
import static main.utils.Constants.*;

public class prediction_utils {

    /*
     * Read GOR raw count matrices from a file written by write_gor_matrices_to_file.
     */
    public static int[][][][] read_matrices_from_file(String filename) {
        int[][][][] gorMatrices = new int[3][20][20][17];

        try {
            Path modelPath = Paths.get(filename);
            String[] lines = Files.readString(modelPath).split("\n");

            int currentState = -1;
            int currentCenterAA = -1;
            int currentOffsetAA;

            for (String line : lines) {

                line = line.trim();
                if (line.isEmpty()) {
                    continue;
                }

                if (line.startsWith("=") && line.endsWith("=")) {
                    currentCenterAA = AA_ALPHABET.indexOf(line.charAt(1));
                    currentState = STATES.indexOf(line.charAt(3));
                }
                if (currentState == -1 || currentCenterAA == -1) {
                    continue;
                }

                String[] values = line.split("\t");
                if (values.length != 18) {
                    continue;
                }

                currentOffsetAA = AA_ALPHABET.indexOf(values[0]);
                if (currentOffsetAA == -1) {
                    continue;
                }

                for (int m = 0; m < 17; m++) {
                    gorMatrices[currentState][currentCenterAA][currentOffsetAA][m] = Integer.parseInt(values[m + 1]);
                }
            }

            return gorMatrices;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }


    public static double calculate_information_content_score(double fSAmaa, double fNotSAmaa, double fSaa, double fNotSaa) {
        // information content score formula
        // I(S, cA, oA, m) = log( f(S, cA, oA, m) / f(not S, cA, oA, m) ) + log( f(not S, cA) / f(S, cA) )

        // pseudocounts, avoid division by zero and log of zero
        fSAmaa += GOR3_PSEUDOCOUNT;
        fNotSAmaa += GOR3_PSEUDOCOUNT;
        fSaa += GOR3_PSEUDOCOUNT;
        fNotSaa += GOR3_PSEUDOCOUNT;

        return Math.log(fSAmaa / fNotSAmaa) + Math.log(fNotSaa / fSaa);
    }

    public static double calculatePSa(double fSAmaa, double fNotSAmaa) {
        return (fSAmaa + GOR3_PSEUDOCOUNT) / (fSAmaa + fNotSAmaa + 2 * GOR3_PSEUDOCOUNT);
    }

    public static double[][][][] counts_to_scores(int[][][][] counts) {
        return counts_to_scores(counts, false);
    }

    /* gets the information content scores from the raw counts, gives a 3d array of all of the scores for each state, amino acid and position */
    /*
    public static double[][][][] counts_to_scores(int[][][][] counts, boolean usePSa) {

        double[][][][] scores = new double[3][20][20][17];

        // precompute fS and fNotS for each state and amino acid, can do this because only need to do it for center position
        int[][] fSaa = new int[3][20]; // how often state S is found at center for each amino acid
        int[][] fNotSaa = new int[3][20]; // how often states other than S are found at center for each amino acid

        for (int stateIndex = 0; stateIndex < 3; stateIndex++) { // for each state
            for (int centerAA = 0; centerAA < 20; centerAA++) { // for each amino acid, count how many times that amino acid is at the center position
                for (int offsetAA = 0; offsetAA < 20; offsetAA++) { // sum all pairs
                    fSaa[stateIndex][centerAA] += counts[stateIndex][centerAA][offsetAA][8];
                    for (int otherStateIndex = 0; otherStateIndex < 3; otherStateIndex++) { // sum other two states in the same way
                        if (otherStateIndex != stateIndex) {
                            fNotSaa[stateIndex][centerAA] += counts[otherStateIndex][centerAA][offsetAA][8];
                        }
                    }
                }
            }
        }

        for (int stateIndex = 0; stateIndex < 3; stateIndex++) { // for each state
            for (int centerAA = 0; centerAA < 20; centerAA++) { // for every state, center aa combination calculate
                for (int offsetAA = 0; offsetAA < 20; offsetAA++) { // the information content that every pair of center to other aa has
                    for (int m = 0; m < 17; m++) {
                        /*if (m == 8) {
                            // single-residue term: only stored for offsetAA == centerAA; log(fSAAj / fNotSAAj) [per slide should add global background log(fNotS / fS), however since this is added to all amino acids and doesn't make a difference as it is added equally for each AA]
                            if (offsetAA == centerAA) { // note: during inference the center column should only ever be accessed for the current amino acid, but this extra check doesn't hurt much
                                double fS = fSaa[stateIndex][centerAA] + GOR3_PSEUDOCOUNT;
                                double fNotS = fNotSaa[stateIndex][centerAA] + GOR3_PSEUDOCOUNT;
                                scores[stateIndex][centerAA][offsetAA][8] = Math.log(fS / fNotS); // center term in formula is not
                            }
                        } else {
                            // regular calculation inside sum
                            int fSAmaa = counts[stateIndex][centerAA][offsetAA][m];
                            int fNotSAmaa = 0;
                            for (int otherStateIndex = 0; otherStateIndex < 3; otherStateIndex++) {
                                if (otherStateIndex != stateIndex) {
                                    fNotSAmaa += counts[otherStateIndex][centerAA][offsetAA][m];
                                }
                            }
                            scores[stateIndex][centerAA][offsetAA][m] = usePSa ? calculatePSa(fSAmaa, fNotSAmaa) : calculate_information_content_score(fSAmaa, fNotSAmaa, fSaa[stateIndex][centerAA], fNotSaa[stateIndex][centerAA]);
                        }*/
                        /*
                        if (m == 8) {
                            // Standard GOR1-style score for the center
                            scores[stateIndex][centerAA][offsetAA][8] = Math.log(fS / fNotS) + Math.log(globalNotS / globalS);
                        } else {
                            // The "Extra" information provided by the neighbor
                            // I(S; Rj, Rj+m) - I(S; Rj)
                            double pairInfo = Math.log(fSAmaa / fNotSAmaa);
                            double centerInfo = Math.log(fS / fNotS);
                            scores[stateIndex][centerAA][offsetAA][m] = pairInfo - centerInfo;
                        }
                    }
                }
            }
        }

        return scores;
    }
    */

    public static double[][][][] counts_to_scores(int[][][][] counts, boolean usePSa) {
        double[][][][] scores = new double[3][20][20][17];

        // 1. Precompute frequencies for normalization
        int[] globalS = new int[3];
        int globalTotal = 0;
        int[][] fSaa = new int[3][20];
        int[][] fNotSaa = new int[3][20];

        for (int s = 0; s < 3; s++) {
            for (int cAA = 0; cAA < 20; cAA++) {
                // At m=8, the count is effectively f(S, centerAA)
                // We sum over offsetAA because center row is sparse (only offset==center is non-zero)
                for (int oAA = 0; oAA < 20; oAA++) {
                    int c = counts[s][cAA][oAA][8];
                    fSaa[s][cAA] += c;
                    globalS[s] += c;
                    globalTotal += c;
                }
            }
        }

        for (int s = 0; s < 3; s++) {
            for (int cAA = 0; cAA < 20; cAA++) {
                for (int otherS = 0; otherS < 3; otherS++) {
                    if (s != otherS) fNotSaa[s][cAA] += fSaa[otherS][cAA];
                }
            }
        }

        // 2. Calculate scores using Information Theory
        for (int s = 0; s < 3; s++) {
            double probS = (double) globalS[s] / globalTotal;
            double probNotS = 1.0 - probS;

            for (int cAA = 0; cAA < 20; cAA++) {
                // Background information for the center residue: I(S; Rj)
                double fS = fSaa[s][cAA] + GOR3_PSEUDOCOUNT;
                double fNotS = fNotSaa[s][cAA] + GOR3_PSEUDOCOUNT;
                double infoCenter = Math.log(fS / fNotS) + Math.log(probNotS / probS);

                for (int oAA = 0; oAA < 20; oAA++) {
                    for (int m = 0; m < 17; m++) {
                        if (m == 8) {
                            // Only assign center info to the diagonal to avoid overcounting
                            scores[s][cAA][oAA][8] = (cAA == oAA) ? infoCenter : 0;
                        } else {
                            // Pair information: I(S; Rj, Rj+m)
                            int fSAmaa = counts[s][cAA][oAA][m];
                            int fNotSAmaa = 0;
                            for (int otherS = 0; otherS < 3; otherS++) {
                                if (s != otherS) fNotSAmaa += counts[otherS][cAA][oAA][m];
                            }

                            if (usePSa) {
                                scores[s][cAA][oAA][m] = (fSAmaa + GOR3_PSEUDOCOUNT) /
                                                         (fSAmaa + fNotSAmaa + 2 * GOR3_PSEUDOCOUNT);
                            } else {
                                double infoPair = Math.log((fSAmaa + GOR3_PSEUDOCOUNT) /
                                                           (fNotSAmaa + GOR3_PSEUDOCOUNT)) +
                                                  Math.log(probNotS / probS);

                                // GOR III key: Neighbor info is infoPair MINUS infoCenter
                                scores[s][cAA][oAA][m] = infoPair - infoCenter;
                            }
                        }
                    }
                }
            }
        }
        return scores;
    }

    /*
     * Predict the secondary structure for a single residue based on a given
     * window, gives a list of all of the scores for the 3 states
     */
    public static double[] predictSingleResidue(String window, double[][][][] scoringMatrices) {

        double[] scores = new double[3];
        char centerAA = window.charAt(8);
        int centerAAIndex = AA_ALPHABET.indexOf(centerAA);

        if (centerAAIndex == -1) {
            return new double[]{0.0, 0.0, 0.0};
        }

        // go through the states matrices
        for (int stateIndex = 0; stateIndex < 3; stateIndex++) {

            // go through the positions in the window, third dim of the matrix
            for (int m = 0; m < 17; m++) {

                // get offset aa
                char offsetAA = window.charAt(m);
                int offsetAAIndex = AA_ALPHABET.indexOf(offsetAA);

                // ignore invalid
                if (offsetAAIndex == -1) {
                    continue;
                }

                double score = scoringMatrices[stateIndex][centerAAIndex][offsetAAIndex][m];

                scores[stateIndex] += score;

            }
        }
        return scores;

    }

    /*
    Calculates probabilities for a single sequence. The output double[][] array only has predictions for valid positions and has length len(fullSequence) - 16
     */
    public static double[][] predictSequence(String fullSequence, double[][][][] scoringMatrices) {
        double[][] scores = new double[fullSequence.length() - 16][3];

        for (int i = 8; i < fullSequence.length() - 8; i++) {
            String window = fullSequence.substring(i - 8, i + 9);
            scores[i - 8] = predictSingleResidue(window, scoringMatrices);
        }

        return scores;
    }

    /*
     * get preds for the secondary structure for full sequence, returns string of
     * predicted labels with probabilities
     */
    public static Pair<ArrayList<String>, ArrayList<String>> predictSequenceWithProbabilities(String fullSequence, double[][][][] scoringMatrices) {

        ArrayList<String> aaAndPSLines = new ArrayList<>();
        ArrayList<String> probLines = new ArrayList<>();

        String formattedSequence = "AS " + fullSequence;
        StringBuilder predictedLabels = new StringBuilder();
        predictedLabels.append("PS --------");
        StringBuilder probabilitiesH = new StringBuilder();
        //probabilitiesH.append("PH --------");
        StringBuilder probabilitiesE = new StringBuilder();
        //probabilitiesE.append("PE --------");
        StringBuilder probabilitiesC = new StringBuilder();
        //probabilitiesC.append("PC --------");

        for (int i = 8; i < fullSequence.length() - 8; i++) {

            String window = fullSequence.substring(i - 8, i + 9);
            double[] scores = predictSingleResidue(window, scoringMatrices);

            // probabilities for states normed to 0-9 via softmax
            double expC = Math.exp(scores[0]);
            double expE = Math.exp(scores[1]);
            double expH = Math.exp(scores[2]);
            double sumExp = expH + expE + expC;

            probabilitiesH.append((int) Math.round((expH / sumExp) * 9));
            probabilitiesE.append((int) Math.round((expE / sumExp) * 9));
            probabilitiesC.append((int) Math.round((expC / sumExp) * 9));

            // for the PS line
            // find the state with the max score
            int maxIndex = 0;
            for (int j = 1; j < scores.length; j++) {
                if (scores[j] > scores[maxIndex]) {
                    maxIndex = j;
                }
            }
            // get label from index
            predictedLabels.append(STATES.charAt(maxIndex));
        }

        predictedLabels.append("--------");
        //probabilitiesH.append("--------");
        //probabilitiesE.append("--------");
        //probabilitiesC.append("--------");

        aaAndPSLines.add(formattedSequence);
        aaAndPSLines.add(predictedLabels.toString());
        probLines.add(probabilitiesH.toString());
        probLines.add(probabilitiesE.toString());
        probLines.add(probabilitiesC.toString());

        return new Pair<>(aaAndPSLines, probLines);
    }

}
