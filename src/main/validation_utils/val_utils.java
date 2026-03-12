package main.validation_utils;

import static main.utils.Math2.max;
import static main.utils.Math2.mean;
import static main.utils.Math2.median;
import static main.utils.Math2.min;
import static main.utils.Math2.quantile;
import static main.utils.Math2.standard_deviation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

import main.utils.Math2;
import main.utils.Pair;

public class val_utils {

    public static final String STATES = "CEH";

    public static double q3(String prediction, String observed) {
        int correct = 0;
        int total = prediction.length();

        for (int i = 0; i < total; i++) {
            if (prediction.charAt(i) == observed.charAt(i)) {
                correct++;
            }
        }
        double res = (double) correct / total * 100;
        return res;
    }

    public static double q3_i(String prediction, String observed, String state) {

        int currentStateidx = STATES.indexOf(state);

        int predictedCount = 0;
        int totalCount = 0;

        for (int i = 0; i < prediction.length(); i++) {

            if (STATES.indexOf(observed.charAt(i)) == currentStateidx) {
                totalCount++;
                if (STATES.indexOf(prediction.charAt(i)) == currentStateidx) {
                    predictedCount++;
                }
            }
        }
        double res = (double) predictedCount / totalCount * 100;
        return res;
    }


    public static double q3_i_avg(ArrayList<String> predictions, ArrayList<String> observed, String state) {
        int currentStateidx = STATES.indexOf(state);

        int predictedCount = 0;
        int totalCount = 0;

        for (int i = 0; i < predictions.size(); i++) {

            if (STATES.indexOf(observed.get(i).charAt(0)) == currentStateidx) {
                totalCount++;
                if (STATES.indexOf(predictions.get(i).charAt(0)) == currentStateidx) {
                    predictedCount++;
                }
            }
        }
        double res = (double) predictedCount / totalCount * 100;
        return res;
    }

    /* calculate the overall Q3 score by averaging the Q3_i scores for each state */
    public static double q3_avg(ArrayList<String> predictions, ArrayList<String> observed) {

        double q3sum = 0;

        for (int i = 0; i < STATES.length(); i++) {
            double q3Score = q3_i_avg(predictions, observed, String.valueOf(STATES.charAt(i)));
            q3sum += q3Score;
        }
        double res = q3sum / STATES.length();
        return res;
    }



    public static ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> getOverlappingPairs(String observedSecondaryStructure, String predictedSecondaryStructure, char targetState) {
        // arraylist of overlapping pairs: Pair of s1, s2; s1 = startIndex, stopIndex; s2 = startIndex, stopIndex (stopIndex is exclusive)
        ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> overlappingPairs = new ArrayList<>();

        // start by getting list of secondary structure in both strands
        ArrayList<Pair<Integer, Integer>> observedStructures = new ArrayList<>();
        ArrayList<Pair<Integer, Integer>> predictedStructures = new ArrayList<>();

        int oStartIndex = 0;
        int pStartIndex = 0;
        boolean oInSS = false;
        boolean pInSS = false;
        for (int sequenceIndex = 0; sequenceIndex < observedSecondaryStructure.length(); sequenceIndex++) {
            char oSS = observedSecondaryStructure.charAt(sequenceIndex);
            char pSS = predictedSecondaryStructure.charAt(sequenceIndex);

            // if strands are not in target secondary structure, open them
            if (!oInSS && oSS == targetState) {
                oStartIndex = sequenceIndex;
                oInSS = true;
            }

            if (!pInSS && pSS == targetState) {
                pStartIndex = sequenceIndex;
                pInSS = true;
            }

            // check for closing secondary structures
            if (oInSS && oSS != targetState) {
                observedStructures.add(new Pair<>(oStartIndex, sequenceIndex));
                oStartIndex = sequenceIndex;
                oInSS = false;
            }

            if (pInSS && pSS != targetState) {
                predictedStructures.add(new Pair<>(pStartIndex, sequenceIndex));
                pStartIndex = sequenceIndex;
                pInSS = false;
            }
        }

        // get final overlapping pair
        if (oInSS) observedStructures.add(new Pair<>(oStartIndex, observedSecondaryStructure.length()));
        if (pInSS) predictedStructures.add(new Pair<>(pStartIndex, observedSecondaryStructure.length()));
        
        int oStopIndex = 0; //oStartIndex already exists
        int pStopIndex = 0;
        // check for overlapping secondary structures
        for (int observedStructuresIndex = 0; observedStructuresIndex < observedStructures.size(); observedStructuresIndex++) {
            oStartIndex = observedStructures.get(observedStructuresIndex).getV1();
            oStopIndex = observedStructures.get(observedStructuresIndex).getV2();
            
            for (int predictedStructuresIndex = 0; predictedStructuresIndex < predictedStructures.size(); predictedStructuresIndex++) {
                pStartIndex = predictedStructures.get(predictedStructuresIndex).getV1();
                pStopIndex = predictedStructures.get(predictedStructuresIndex).getV2();

                // check for overlap
                if ((oStartIndex >= pStartIndex && oStartIndex < pStopIndex) || (pStartIndex >= oStartIndex && pStartIndex < oStopIndex)) {
                    overlappingPairs.add(new Pair<>(observedStructures.get(observedStructuresIndex), predictedStructures.get(predictedStructuresIndex)));
                }
            }
        }

        return overlappingPairs;
    }

    public static ArrayList<Pair<Integer, Integer>> getLonelyObserved(String observedSecondaryStructure, String predictedSecondaryStructure, char targetState) {
        ArrayList<Pair<Integer, Integer>> lonelyStructures = new ArrayList<>();
        
        // start by getting list of secondary structure in both strands
        ArrayList<Pair<Integer, Integer>> observedStructures = new ArrayList<>();
        ArrayList<Pair<Integer, Integer>> predictedStructures = new ArrayList<>();

        int oStartIndex = 0;
        int pStartIndex = 0;
        boolean oInSS = false;
        boolean pInSS = false;
        for (int sequenceIndex = 0; sequenceIndex < observedSecondaryStructure.length(); sequenceIndex++) {
            char oSS = observedSecondaryStructure.charAt(sequenceIndex);
            char pSS = predictedSecondaryStructure.charAt(sequenceIndex);

            // if strands are not in target secondary structure, open them
            if (!oInSS && oSS == targetState) {
                oStartIndex = sequenceIndex;
                oInSS = true;
            }

            if (!pInSS && pSS == targetState) {
                pStartIndex = sequenceIndex;
                pInSS = true;
            }

            // check for closing secondary structures
            if (oInSS && oSS != targetState) {
                observedStructures.add(new Pair<>(oStartIndex, sequenceIndex));
                oStartIndex = sequenceIndex;
                oInSS = false;
            }

            if (pInSS && pSS != targetState) {
                predictedStructures.add(new Pair<>(pStartIndex, sequenceIndex));
                pStartIndex = sequenceIndex;
                pInSS = false;
            }
        }

        // get final overlapping pair
        if (oInSS) observedStructures.add(new Pair<>(oStartIndex, observedSecondaryStructure.length()));
        if (pInSS) predictedStructures.add(new Pair<>(pStartIndex, observedSecondaryStructure.length()));
        
        int oStopIndex = 0; //oStartIndex already exists
        int pStopIndex = 0;
        // check for overlapping secondary structures
        for (int observedStructuresIndex = 0; observedStructuresIndex < observedStructures.size(); observedStructuresIndex++) {
            oStartIndex = observedStructures.get(observedStructuresIndex).getV1();
            oStopIndex = observedStructures.get(observedStructuresIndex).getV2();

            boolean foundOverlap = false;
            
            for (int predictedStructuresIndex = 0; predictedStructuresIndex < predictedStructures.size(); predictedStructuresIndex++) {
                pStartIndex = predictedStructures.get(predictedStructuresIndex).getV1();
                pStopIndex = predictedStructures.get(predictedStructuresIndex).getV2();

                // check for overlap
                if ((oStartIndex >= pStartIndex && oStartIndex < pStopIndex) || (pStartIndex >= oStartIndex && pStartIndex < oStopIndex)) {
                    foundOverlap = true;
                    break;
                }
            }

            if (!foundOverlap) lonelyStructures.add(new Pair<>(oStartIndex, oStopIndex));
        }

        return lonelyStructures;
    }

    public static int maxov(Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> overlappingPair) {
        int oStart = overlappingPair.getV1().getV1();
        int oStop = overlappingPair.getV1().getV2();
        int pStart = overlappingPair.getV2().getV1();
        int pStop = overlappingPair.getV2().getV2();

        int minStart = oStart < pStart ? oStart : pStart;
        int maxStop = oStop > pStop ? oStop : pStop;

        return maxStop - minStart;
    }

    public static int minov(Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> overlappingPair) {
        int oStart = overlappingPair.getV1().getV1();
        int oStop = overlappingPair.getV1().getV2();
        int pStart = overlappingPair.getV2().getV1();
        int pStop = overlappingPair.getV2().getV2();

        int maxStart = oStart > pStart ? oStart : pStart;
        int minStop = oStop < pStop ? oStop : pStop;

        return minStop - maxStart;
    }

    public static int sovDelta(Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> overlappingPair) {
        ArrayList<Double> values = new ArrayList<>();
        values.add((double) (maxov(overlappingPair) - minov(overlappingPair)));
        values.add((double) (minov(overlappingPair)));
        values.add((double) Math.round((overlappingPair.getV1().getV2() - overlappingPair.getV1().getV1()) / 2));
        values.add((double) Math.round((overlappingPair.getV2().getV2() - overlappingPair.getV2().getV1()) / 2));

        return (int) Math2.min(values);
    }

    public static int sovN_i(ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> overlappingPairs, ArrayList<Pair<Integer, Integer>> lonelyStructures) {
        int sum = 0;
        for (Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> overlappingPair : overlappingPairs) {
            sum += overlappingPair.getV1().getV2() - overlappingPair.getV1().getV1();
        }

        for (Pair<Integer, Integer> lonelyStructure : lonelyStructures) {
            sum += lonelyStructure.getV2() - lonelyStructure.getV1();
        }

        return sum;
    } 

    public static double sovLargeSum_i(ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> overlappingPairs) {
        double sum = 0.0;

        for (Pair<Pair<Integer, Integer>, Pair<Integer, Integer>> overlappingPair : overlappingPairs) {
            double v = (((double) (minov(overlappingPair) + sovDelta(overlappingPair))) / ((double) maxov(overlappingPair))) * ((double) (overlappingPair.getV1().getV2() - overlappingPair.getV1().getV1()));
            sum += v;
        }

        return sum;
    }

    public static double sov_i(String observedSecondaryStructure, String predictedSecondaryStructure, char targetState) {
        ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>> overlappingPairs = getOverlappingPairs(observedSecondaryStructure, predictedSecondaryStructure, targetState);
        ArrayList<Pair<Integer, Integer>> lonelyStructures = getLonelyObserved(observedSecondaryStructure, predictedSecondaryStructure, targetState);
        if (overlappingPairs.isEmpty()) return 0; // skip division by 0
        return 100.0 * (1.0/ (double) sovN_i(overlappingPairs, lonelyStructures)) * ((double) sovLargeSum_i(overlappingPairs));
    }

    public static double sov(String observedSecondaryStructure, String predictedSecondaryStructure) {

        // calculate overlaps and lonelies across all three states
        ArrayList<ArrayList<Pair<Pair<Integer, Integer>, Pair<Integer, Integer>>>> overlappingPairs = new ArrayList<>();
        ArrayList<ArrayList<Pair<Integer, Integer>>> lonelyStructures = new ArrayList<>();
        for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
            overlappingPairs.add(getOverlappingPairs(observedSecondaryStructure, predictedSecondaryStructure, STATES.charAt(stateIndex)));
            lonelyStructures.add(getLonelyObserved(observedSecondaryStructure, predictedSecondaryStructure, STATES.charAt(stateIndex)));
        }

        
        double score = 0;
        int N = 0;

        // calculate 1/N term in formula
        for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
            N += sovN_i(overlappingPairs.get(stateIndex), lonelyStructures.get(stateIndex));
        }

        score += 100.0 / (double) N;

        // calculate large sum
        double rearSum = 0;
        for (int stateIndex = 0; stateIndex < 3; stateIndex++) {
            rearSum += sovLargeSum_i(overlappingPairs.get(stateIndex));
        }

        return score * rearSum;
    }

    // needs >id <Q3> <SOV> <QH> <QE> <QC> <SOV_H> <SOV_E> <SOV_C> as header
    public static void write_detailed_file(String filePath, ArrayList<String> ids, ArrayList<String> sequences, ArrayList<String> predictions, ArrayList<String> observed) {

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filePath))) {

            for (int i = 0; i < predictions.size(); i++) {

                String id = ids.get(i);
                String sequence = sequences.get(i);
                String prediction = knabber16(predictions.get(i));
                if (prediction == null) continue; // should skip if seq too short for prediction
                String observe = knabber16(observed.get(i));

                double q3Score = q3(prediction, observe);
                double sovScore = sov(observe, prediction);

                double q3_h_score = q3_i(prediction, observe, "H");
                double q3_e_score = q3_i(prediction, observe, "E");
                double q3_c_score = q3_i(prediction, observe, "C");

                double sov_h_score = sov_i(observe, prediction, 'H');
                double sov_e_score = sov_i(observe, prediction, 'E');
                double sov_c_score = sov_i(observe, prediction, 'C');

                StringBuilder res = new StringBuilder();

                DecimalFormat df = new DecimalFormat("0.0");

                res.append("\n");
                res.append("> " + id + " " + padString(df.format(q3Score)) + padString(df.format(sovScore)) + padString(df.format(q3_h_score)) + padString(df.format(q3_e_score)) + padString(df.format(q3_c_score)) + padString(df.format(sov_h_score)) + padString(df.format(sov_e_score)) + padString(df.format(sov_c_score)) + "\n");
                res.append("AS " + sequence + "\n");
                res.append("PS " + predictions.get(i) + "\n");
                res.append("SS " + observed.get(i) + "\n");

                bw.write((res.toString() + "\n").replaceAll(" NaN  ", "  -   "));

            }

        } catch (IOException e) {
            System.err.println("Error writing detailed file: " + e.getMessage());
        }
    }

    /* write detailed file with the stats */
    public static void write_summary_file(String filePath, ArrayList<String> ids, ArrayList<String> predictions, ArrayList<String> observed) {

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filePath))) {

            ArrayList<ArrayList<Double>> all_scores = new ArrayList<>();
            DecimalFormat df = new DecimalFormat("0.0");

            int numProteins = predictions.size();
            int sumProteinLength = 0;
            int sumPredictedPositions = 0;

            ArrayList<Double> q3Scores = new ArrayList<>();
            ArrayList<Double> q3_h_scores = new ArrayList<>();
            ArrayList<Double> q3_e_scores = new ArrayList<>();
            ArrayList<Double> q3_c_scores = new ArrayList<>();

            ArrayList<Double> sovScores = new ArrayList<>();
            ArrayList<Double> sov_h_scores = new ArrayList<>();
            ArrayList<Double> sov_e_scores = new ArrayList<>();
            ArrayList<Double> sov_c_scores = new ArrayList<>();

            for (int i = 0; i < predictions.size(); i++) {

                String prediction = predictions.get(i);
                if (prediction == null) continue; // should skip if seq too short for prediction

                String observe = observed.get(i);

                q3Scores.add(q3(prediction, observe));
                sovScores.add(sov(observe, prediction));

                q3_h_scores.add(Double.valueOf(df.format(q3_i(prediction, observe, "H"))));
                q3_e_scores.add(Double.valueOf(df.format(q3_i(prediction, observe, "E"))));
                q3_c_scores.add(Double.valueOf(df.format(q3_i(prediction, observe, "C"))));

                sov_h_scores.add(Double.valueOf(df.format(sov_i(observe, prediction, 'H'))));
                sov_e_scores.add(Double.valueOf(df.format(sov_i(observe, prediction, 'E'))));
                sov_c_scores.add(Double.valueOf(df.format(sov_i(observe, prediction, 'C'))));

                sumProteinLength += prediction.length();
                sumPredictedPositions += prediction.replace("-", "").length();
            }

            double meanProteinLength = (double) sumProteinLength / numProteins;

            all_scores.add(q3Scores);
            all_scores.add(q3_h_scores);
            all_scores.add(q3_e_scores);
            all_scores.add(q3_c_scores);
            all_scores.add(sovScores);
            all_scores.add(sov_h_scores);
            all_scores.add(sov_e_scores);
            all_scores.add(sov_c_scores);

            ArrayList<ArrayList<String>> stats = get_stats_from_scores(all_scores);

            StringBuilder res = new StringBuilder();
            res.append("Number of proteins:          " + numProteins + "\n");
            res.append("Mean protein length:         " + df.format(Double.valueOf(meanProteinLength)) + "\n");
            res.append("Sum of Protein Length:       " + sumProteinLength + "\n");
            res.append("Sum of Predicted Positions:  " + sumPredictedPositions + "\n");

            for (ArrayList<String> stat : stats) {
                for (String s : stat) {
                    res.append(s);
                }
                res.append("\n");
            }

            bw.write(res.toString());

            bw.close();

        } catch (IOException e) {

            System.err.println("Error writing summary file" + e.getMessage());

        }
    }

    // needs >id <Q3> <SOV> <QH> <QE> <QC> <SOV_H> <SOV_E> <SOV_C> as header
    public static ArrayList<ArrayList<Double>> parse_detailed_file_for_metrics(String filePath) {

        try {

            ArrayList<Double> q3Scores = new ArrayList<>();
            ArrayList<Double> q3_h_scores = new ArrayList<>();
            ArrayList<Double> q3_e_scores = new ArrayList<>();
            ArrayList<Double> q3_c_scores = new ArrayList<>();
            ArrayList<Double> sovScores = new ArrayList<>();
            ArrayList<Double> sov_h_scores = new ArrayList<>();
            ArrayList<Double> sov_e_scores = new ArrayList<>();
            ArrayList<Double> sov_c_scores = new ArrayList<>();

            BufferedReader br = new BufferedReader(new FileReader(filePath));

            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    line = line.substring(1).trim();
                    String[] parts = line.split("\\s+");
                    if (!parts[1].equals("-")) q3Scores.add(Double.valueOf(parts[1]));
                    if (!parts[2].equals("-")) sovScores.add(Double.valueOf(parts[2]));
                    if (!parts[3].equals("-")) q3_h_scores.add(Double.valueOf(parts[3]));
                    if (!parts[4].equals("-")) q3_e_scores.add(Double.valueOf(parts[4]));
                    if (!parts[5].equals("-")) q3_c_scores.add(Double.valueOf(parts[5]));
                    if (!parts[6].equals("-")) sov_h_scores.add(Double.valueOf(parts[6]));
                    if (!parts[7].equals("-")) sov_e_scores.add(Double.valueOf(parts[7]));
                    if (!parts[8].equals("-")) sov_c_scores.add(Double.valueOf(parts[8]));
                }
            }

            br.close();

            return new ArrayList<ArrayList<Double>>() {{
                add(q3Scores);
                add(q3_h_scores);
                add(q3_e_scores);
                add(q3_c_scores);
                add(sovScores);
                add(sov_h_scores);
                add(sov_e_scores);
                add(sov_c_scores);
            }};

        } catch (IOException e) {
            System.err.println("Error reading detailed file: " + e.getMessage());
            return null;
        }

    }

    public static void write_summary_file_from_detailed_file(String filePath, ArrayList<String> ids, ArrayList<String> predictions, ArrayList<String> observed, String detailedFile) {

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(filePath))) {

            StringBuilder res = new StringBuilder();
            DecimalFormat df = new DecimalFormat("0.0");

            int numProteins = predictions.size();
            int sumProteinLength = 0;
            int sumPredictedPositions = 0;

            for (int i = 0; i < predictions.size(); i++) {
                if (predictions.get(i) == null) continue;
                sumProteinLength += predictions.get(i).length();
                sumPredictedPositions += predictions.get(i).replace("-", "").length();
            }

            double meanProteinLength = (double) sumProteinLength / numProteins;

            res.append("Number of proteins:          " + numProteins + "\n");
            res.append("Mean protein length:         " + df.format(meanProteinLength) + "\n");
            res.append("Sum of Protein Length:       " + sumProteinLength + "\n");
            res.append("Sum of Predicted Positions:  " + sumPredictedPositions + "\n");

            ArrayList<ArrayList<Double>> all_scores = parse_detailed_file_for_metrics(detailedFile);
            ArrayList<ArrayList<String>> stats = get_stats_from_scores(all_scores);

            for (ArrayList<String> stat : stats) {
                for (String s : stat) {
                    res.append(s);
                }
                res.append("\n");
            }

            bw.write(res.toString());

        } catch (IOException e) {
            System.err.println("Error writing summary file: " + e.getMessage());
        }
    }

    public static ArrayList<ArrayList<String>> get_stats_from_scores(ArrayList<ArrayList<Double>> all_scores) {

        ArrayList<String> scoreNames = new ArrayList<>();
        scoreNames.add("q3");
        scoreNames.add("qObs_H");
        scoreNames.add("qObs_E");
        scoreNames.add("qObs_C");
        scoreNames.add("SOV");
        scoreNames.add("SOV_H");
        scoreNames.add("SOV_E");
        scoreNames.add("SOV_C");

        ArrayList<ArrayList<String>> results = new ArrayList<>();
        DecimalFormat df = new DecimalFormat("0.0");

        for (int i = 0; i < all_scores.size(); i++) {

            ArrayList<String> currentResults = new ArrayList<>();

            String currentScoreName = scoreNames.get(i);
            ArrayList<Double> currentScores = all_scores.get(i);

            currentResults.add(currentScoreName + "\t");
            currentResults.add("Mean:" + (df.format(mean(currentScores)) + "\t"));
            currentResults.add("Dev:" + (df.format(standard_deviation(currentScores)) + "\t"));
            currentResults.add("Min:" + (df.format(min(currentScores)) + "\t"));
            currentResults.add("Max:" + (df.format(max(currentScores)) + "\t"));
            currentResults.add("Median:" + (df.format(median(currentScores)) + "\t"));
            currentResults.add("Quantil_25:" + (df.format(quantile(currentScores, 0.25)) + "\t"));
            currentResults.add("Quantil_75:" + (df.format(quantile(currentScores, 0.75)) + "\t"));
            currentResults.add("Quantil_5:" + (df.format(quantile(currentScores, 0.05)) + "\t"));
            currentResults.add("Quantil_95:" + (df.format(quantile(currentScores, 0.95)) + ""));

            results.add(currentResults);

        }
        return results;
    }

    
    public static ArrayList<String> knabber16(ArrayList<String> inputSequences) {
        ArrayList<String> knabberedSequences = new ArrayList<>();
        for (String sequence : inputSequences) {
            knabberedSequences.add(sequence.length() > 16 ? sequence.substring(8, sequence.length() - 8) : null);
        }

        return knabberedSequences;
    }

    public static String knabber16(String sequence) {
        return sequence.length() > 16 ? sequence.substring(8, sequence.length() - 8) : null;
    }

    // pad string helper for the detailed output (terrible format)
    public static String padString(String input) {
        if (input == null) return "      ";
        int len = input.length();

        switch (len) {
            case 1:
                return "  " + input + "   ";
            case 2:
                return "  " + input + "  ";
            case 3:
                return " " + input + "  ";
            case 4:
                return " " + input + " ";
            case 5:
                return "" + input + " ";
            default:
                return input;
        }
    }

    /*
    public static void main(String[] args) {
        String p = "--------HHCECCHEEHHHHCECECECCEHHHHHHHHHHEHHHHHCCHHHHHHHHHHCCCCHHHHEHHHHHEHHHHCCCEECECECHEEEECCCEECCCCCECECCEEEHHHHEHHHHHHHHHHHHHHHCECCCHHHHHHHHHHEHHCHEHHCCECCCCECHCHHHEHHHHHHHHH--------";
        String o = "CCCCCCCHHHCCCCCECHHHHCCCCCCCCEHHHHHHHHHHCHHHHHCCHHHHHHHHHHHCCCHHHHHHHHHHHHHHHCCCECCECCCCCEECCCCEECCCCCCCCCCCCHHHHHHHHHHHHHHHHHHHHHCCCCCHHHHHHHHHHHHHHCHHHCCCCCCCCCCCHHHCHHHHHHHHHHHHHHCCC";
        o = o.substring(8, o.length() - 8);
        p = p.substring(8, p.length() - 8);
        //for (int i = 0; i < 3; i++) {
        //    System.out.println(STATES.charAt(i) + " " + sov_i(o, p, STATES.charAt(i)));
        //}
        System.out.println(sov(o, p));
        //System.out.println(sov("CCEEECCCCCCEEEEEECCC", "CCCCCCCEEEEECCCEECCC"));
        
        String o = "CCCHHHCCHHHCEEEEEEECCCCHHHHHHH";
        String p = "CCHHHHCCCHHHCCCEEECCCCCCCEEEEE";
        char state = 'E';
        System.out.println(getOverlappingPairs(o, p, state));
        System.out.println(getLonelyObserved(o, p, state));
        System.out.println(sov_i(o, p, state));
        System.out.println(sov(o, p));
        
    }
    */
}


