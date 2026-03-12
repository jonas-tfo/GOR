package main.utils;

public class post_processing_utils {

    /* applies smoothing to a sequence of labels to replace biologically unrealistic patterns */
    public static String smoothingRegex(String predictedStructure) {
        return predictedStructure.replaceAll("CHC|CEC", "CCC").replaceAll("CHHC", "CCCC").replaceAll("CHHHC", "CCCCC").replaceAll("EHE", "EEE").replaceAll("ECE", "EEE").replaceAll("HEH", "HHH");
    }

    public static String smoothingWindow(String predictedStructure, int windowSize) {

        return predictedStructure;

    }

    public static String smoothingWindow_old(String predictedStructure, int windowSize) {
        StringBuilder smoothedStructure = new StringBuilder(predictedStructure);

        for (int i = windowSize; i < predictedStructure.length() - windowSize; i++) {
            int helixCount = 0;
            int sheetCount = 0;
            int coilCount = 0;

            // count the number of each state in the window
            for (int j = Math.max(0, i - windowSize); j <= Math.min(predictedStructure.length() - 1, i + windowSize); j++) {
                char state = predictedStructure.charAt(j);
                if (state == 'H') {
                    helixCount++;
                } else if (state == 'E') {
                    sheetCount++;
                } else if (state == 'C') {
                    coilCount++;
                }
            }

            // set state of i to most freq state
            if (helixCount > sheetCount && helixCount > coilCount) {
                smoothedStructure.setCharAt(i, 'H');
            } else if (sheetCount > helixCount && sheetCount > coilCount) {
                smoothedStructure.setCharAt(i, 'E');
            } else {
                smoothedStructure.setCharAt(i, 'C');
            }
        }
        return smoothedStructure.toString();
    }

    public static String applySmoothing(String predictedStructure, int windowSize) {
        //return smoothingWindow(predictedStructure, windowSize);
        return smoothingWindow_old(smoothingRegex(predictedStructure), windowSize);
    }
}
