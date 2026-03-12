package main.gor5_utils;

import java.util.ArrayList;

public class prediction_utils {

    public static <T> Double[][] generateProbabilitiesSingleSequence(String fullSequence, predict_residue<T> predictionFunction, T scoringMatrices) {
        Double[][] scores = new Double[fullSequence.length() - 16][3]; // complete score matrix for sequence across 3 states; -16 because no predictions for edges are made

        for (int sequencePosition = 8; sequencePosition < fullSequence.length() - 8; sequencePosition++) {
            // skip if center position is gap
            if (fullSequence.charAt(sequencePosition) == '-') {
                scores[sequencePosition - 8] = null; // -8 to compensate for reduced size compared to full sequence
                continue;
            }
          
            String window = fullSequence.substring(sequencePosition - 8, sequencePosition + 9);
            double[] tmp = predictionFunction.predictResidue(window, scoringMatrices);

            // copy values to reference type // TODO find better method than this
            for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
                scores[sequencePosition - 8][stateIndex] = tmp[stateIndex]; // -8 to compensate for reduced size compared to full sequene
            }

        }

        return scores;
    }

    public static <T> Double[][][] generateProbabilitiesAllSequences(ArrayList<String> alignedSequences, predict_residue<T> predictionFunction, T scoringMatrices) {
        Double[][][] scores = new Double[alignedSequences.size()][alignedSequences.get(0).length() - 16][3]; // -16 because no predictions are made for both edges
        for (int sequenceNumber = 0; sequenceNumber < alignedSequences.size(); sequenceNumber++) {
            scores[sequenceNumber] = generateProbabilitiesSingleSequence(alignedSequences.get(sequenceNumber), predictionFunction, scoringMatrices);
        }

        return scores;
    }

    public static double[] averagePosition(char[] lettersAtPosition, Double[][] predictionsAtPosition) {
        double[] averageProbabilities = new double[3];
        int total = 0;

        for (int sequenceNumber = 0; sequenceNumber < lettersAtPosition.length; sequenceNumber++) {
            // Skip if the residue is a gap OR if the specific sequence had no prediction (null)
            if (lettersAtPosition[sequenceNumber] == '-' || predictionsAtPosition[sequenceNumber] == null) {
                continue;
            }

            for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
                averageProbabilities[stateIndex] += predictionsAtPosition[sequenceNumber][stateIndex];
            }
            total++;
        }

        // Guard against division by zero if a column is all gaps
        if (total > 0) {
            averageProbabilities[0] /= total;
            averageProbabilities[1] /= total;
            averageProbabilities[2] /= total;
        }

        return averageProbabilities;
    }

    public static double[][] averageProbabilities(ArrayList<String> alignedSequences, Double[][][] allProbabilities) {
        double[][] averageProbabilities = new double[allProbabilities[0].length][3];

        // iterate through every position
        for (int predictionPosition = 0; predictionPosition < allProbabilities[0].length; predictionPosition++) {
            // collect residues at this position
            char[] residuesAtPosition = new char[alignedSequences.size()];
            for (int sequenceNumber = 0; sequenceNumber < alignedSequences.size(); sequenceNumber++) {
                residuesAtPosition[sequenceNumber] = alignedSequences.get(sequenceNumber).charAt(predictionPosition + 8); // account for window offset
            }

            // collect predictions for each position from all sequences
            Double[][] predictionsAtPosition = new Double[alignedSequences.size()][3];
            for (int sequenceNumber = 0; sequenceNumber < alignedSequences.size(); sequenceNumber++) {
                predictionsAtPosition[sequenceNumber] = allProbabilities[sequenceNumber][predictionPosition];
            }

            averageProbabilities[predictionPosition] = averagePosition(residuesAtPosition, predictionsAtPosition);
        }

        return averageProbabilities;
    }
}

