package main.utils;

import java.util.HashMap;
import java.util.Objects;

import static main.utils.Constants.*;

public class output_utils {
	public static final String DBLU = "\033[34m";
	public static final String LBLU = "\033[36m";
	public static final String LGRN = "\033[92m";
	public static final String YEL = "\033[33m";
	public static final String ORG = "\033[38;5;208m";
	public static final String RED = "\033[31m";
	public static final String WHT = "\033[37m";
	public static final String GRY = "\033[90m";
	public static final String RESET = "\033[0m";
	public static final String BOLD = "\033[1m";


	public static String colorSecondaryStructure(String secondaryStructureSequence) {
		return secondaryStructureSequence.replaceAll("E", DBLU + "E" + RESET).replaceAll("H", RED + "H" + RESET).replaceAll("C", YEL + "C" + RESET);
	}

	// output formatting functions
	public static String buildSequenceId(String seqId, boolean coloredOutput) {
		return coloredOutput ? BOLD + ">" + seqId + RESET : ">" + seqId;
	}

	public static String getPredictedStructure(double[][] probabilities) {
		StringBuilder predictedStructure = new StringBuilder();


		for (double[] probability : probabilities) {

			// find highest value at current probabilities position
			int maxIndex = 0;
			for (int stateIndex = 1; stateIndex < STATES.length(); stateIndex++) {
				if (probability[stateIndex] > probability[maxIndex]) {
					maxIndex = stateIndex;
				}
			}

			// build predicted structure
			predictedStructure.append(STATES.charAt(maxIndex));
		}
		return predictedStructure.toString();
	}

	public static String buildPrediction(String sequence, double[][] probabilities, boolean applyPostProcessing, PostProcessingFunction postProcessingFunction, boolean coloredOutput) {
		// catch sequences which are too short
		if (probabilities == null) {
			return "-".repeat(sequence.length());
		}

		// postprocessing
		String rawStructure = getPredictedStructure(probabilities);
		String completeStructure = applyPostProcessing ? postProcessingFunction.apply(rawStructure, Constants.smoothingWindowSize) : rawStructure;

		// add non-predicted residues
		return "-".repeat(8) + (coloredOutput ? colorSecondaryStructure(completeStructure) : completeStructure) + "-".repeat(8);
	}

	public static String[] getProbabilityLines(double[][] probabilities, boolean coloredOutput) {
		// init string builders
		StringBuilder[] probabilityLines = new StringBuilder[3];
		for (int i = 0; i < 3; i++) {
			probabilityLines[i] = new StringBuilder();
		}

		// build each probability line
		for (int sequenceIndex = 0; sequenceIndex < probabilities.length; sequenceIndex++) {
			// if colored output is wished --> find highest index and apply coloring
			int maxIndex = 0;
			if (coloredOutput) {
				for (int stateIndex = 0; stateIndex < STATES.length(); stateIndex++) {
					if (probabilities[sequenceIndex][stateIndex] > probabilities[sequenceIndex][maxIndex]) {
						maxIndex = stateIndex;
					}
				}
			}

			// loop through each state and calculate the exponential values
			double[] expS = new double[3];
			double sumExp = 0.0;
			for (int stateIndex = 0; stateIndex < STATES.length(); stateIndex++) {
				expS[stateIndex] = Math.exp(probabilities[sequenceIndex][stateIndex]);
				sumExp += expS[stateIndex];
			}

			// loop through each state, calcualte the rounded values and apply coloring to max
			for (int stateIndex = 0; stateIndex < STATES.length(); stateIndex++) {
				int roundedProbability = (int) Math.round((expS[stateIndex] / sumExp) * 9);
				probabilityLines[stateIndex].append(coloredOutput && stateIndex != maxIndex ? roundedProbability : LGRN + roundedProbability + RESET);
			}
		}

		return new String[]{probabilityLines[0].toString(), probabilityLines[1].toString(), probabilityLines[2].toString()};
	}

	public static String buildProbabilityLines(double[][] probabilities, boolean coloredOutput, char[] probabilityOutputOrder) {
		// get (colored) probability lines in the order of STATES
		String[] probabilityLines = getProbabilityLines(probabilities, coloredOutput);

		StringBuilder output = new StringBuilder();
		for (char state : probabilityOutputOrder) {
			String probabilityLine = "P" + state + " " + "-".repeat(8) + probabilityLines[STATES.indexOf(state)] + "-".repeat(8);
			output.append(probabilityLine).append("\n");
		}

		return output.toString();
	}

	public static String applyFormat(String output, String format) {
		return Objects.equals(format, "html") ? "<pre>" + output + "</pre>" : output;
	}

	// function which builds output for a single sequence given the specified parameters
	public static String buildSingleSequenceOutput(String seqId, String sequence, double[][] probabilities, boolean outputProbabilities, boolean coloredOutput, String format, boolean applyPostProcessing, PostProcessingFunction postProcessingFunction) {
		StringBuilder output = new StringBuilder();

		output.append(buildSequenceId(seqId, coloredOutput)).append("\n");
		output.append("AS ").append(sequence).append("\n");
		output.append("PS ").append(buildPrediction(sequence, probabilities, applyPostProcessing, postProcessingFunction, coloredOutput)).append("\n");
		output.append(outputProbabilities && probabilities != null ? buildProbabilityLines(probabilities, coloredOutput, new char[]{'H', 'E', 'C'}) : ""); // no new line at end, because buildProbability Lines already does this

		return applyFormat(output.toString(), format);
	}

	public static String buildOutput(HashMap<String, Pair<String, double[][]>> idsSequencesProbabilities, boolean outputProbabilities, boolean coloredOutput, String format, boolean applyPostProcessing, PostProcessingFunction postProcessingFunction) {
		StringBuilder output = new StringBuilder();

		for (String seqId : idsSequencesProbabilities.keySet()) {
			String sequence = idsSequencesProbabilities.get(seqId).getV1();
			double[][] probabilities = idsSequencesProbabilities.get(seqId).getV2();
			output.append(buildSingleSequenceOutput(seqId, sequence, probabilities, outputProbabilities, coloredOutput, format, applyPostProcessing, postProcessingFunction)).append("\n");
		}

		return output.toString();
	}
}
