package main.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.stream.Collectors;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.util.stream.Stream;

public class parse_utils {

    public static HashMap<String, Pair<String, String>> parse_seclib(String seclibPath) {
        try {
            HashMap<String, Pair<String, String>> idsSequencesStructures = new HashMap<>();
            
            Path alignmentFilePath = Paths.get(seclibPath);
            String[] fileContent = Files.readString(alignmentFilePath).split(">");

            String seqId = "";
            StringBuilder aminoAcidSequence = new StringBuilder();
            StringBuilder labelSequence = new StringBuilder();
            boolean appendingToAASequence = false;

            for (String entry : fileContent) {
                String[] lines = entry.split("\n");
                seqId = lines[0].trim();
                for (int i = 1; i < lines.length; i++) {
                    String line = lines[i];
                    if (line.startsWith("AS ")) {
                        appendingToAASequence = true;
                        aminoAcidSequence.append(line.split(" ")[1].trim());
                    } else if (line.startsWith("PS ") || line.startsWith("SS ")) {
                        appendingToAASequence = false;
                        labelSequence.append(line.split(" ")[1].trim());
                    } else if (appendingToAASequence && seqId != "") {
                        aminoAcidSequence.append(line.trim());
                    } else if (seqId != ""){
                        labelSequence.append(line.trim());
                    }
                }

                if (seqId != "" && !aminoAcidSequence.isEmpty() && !labelSequence.isEmpty()) {
                    idsSequencesStructures.put(seqId, new Pair<String,String>(aminoAcidSequence.toString(), labelSequence.toString()));
                    seqId = "";
                    aminoAcidSequence = new StringBuilder();
                    labelSequence = new StringBuilder();
                    appendingToAASequence = false;
                }
            }

            return idsSequencesStructures;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    /*
     * Parses a seclib file and returns list of three lists, for the ids,
     * for the amino acid sequences, for the label sequences.
     */
    /*
    public static HashMap<String, Pair<String, String>> parse_seclib(String seclibPath) {

        try {

            FileReader fr = new FileReader(seclibPath);
            BufferedReader br = new BufferedReader(fr);

            //ArrayList<String> ids = new ArrayList<>();
            //ArrayList<String> aminoAcidSequences = new ArrayList<>();
            //ArrayList<String> labelSequences = new ArrayList<>();
            HashMap<String, Pair<String, String>> idsSequencesStructures = new HashMap<>();

            String line;
            String seqId = "";
            String aminoAcidSequence = "";
            String labelSequence = "";

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    //ids.add(line.substring(1).trim());
                    seqId = line.substring(1).trim();
                } else if (line.startsWith("AS ")) {
                    //aminoAcidSequences.add(line.substring(3).trim());
                    aminoAcidSequence = line.substring(3).trim();
                } else if (line.startsWith("SS ")) {
                    //labelSequences.add(line.substring(3).trim());
                    labelSequence = line.substring(3).trim();
                }

                if (seqId != "" && aminoAcidSequence != "" && labelSequence != "") {
                    idsSequencesStructures.put(seqId, new Pair<String,String>(aminoAcidSequence, labelSequence));
                    seqId = "";
                    aminoAcidSequence = "";
                    labelSequence = "";
                }
            }

            br.close();
            fr.close();

            //ArrayList<ArrayList<String>> result = new ArrayList<>();

            //result.add(ids);
            //result.add(aminoAcidSequences);
            //result.add(labelSequences);

            return idsSequencesStructures;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
    */


    public static HashMap<String, String> parse_fasta(String seclibPath) {
        try {
            HashMap<String, String> idsSequences = new HashMap<>();
            
            Path alignmentFilePath = Paths.get(seclibPath);
            String[] fileContent = Files.readString(alignmentFilePath).split(">");

            String seqId = "";
            StringBuilder aminoAcidSequence = new StringBuilder();

            for (String entry : fileContent) {
                String[] lines = entry.split("\n");
                seqId = lines[0];
                for (int i = 1; i < lines.length; i++) {
                    aminoAcidSequence.append(lines[i].trim());
                }

                if (seqId != "" && !aminoAcidSequence.isEmpty()) {
                    idsSequences.put(seqId, aminoAcidSequence.toString());
                    seqId = "";
                    aminoAcidSequence = new StringBuilder();
                }
            }

            return idsSequences;
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }    

    /*
     * Parses the fasta file and returns a list of two lists for the
     * IDs and for the sequences.
     */
    /*
    public static HashMap<String, String> parse_fasta(String fastaPath) {

        try {

            FileReader fr = new FileReader(fastaPath);
            BufferedReader br = new BufferedReader(fr);

            ArrayList<String> ids = new ArrayList<>();
            ArrayList<String> sequences = new ArrayList<>();

            String line;
            StringBuilder currentSequence = new StringBuilder();

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    if (!currentSequence.isEmpty()) {
                        sequences.add(currentSequence.toString());
                        currentSequence.setLength(0);
                    }
                    ids.add(line.substring(1).trim());
                } else {
                    currentSequence.append(line.trim());
                }
            }

            if (!currentSequence.isEmpty()) {
                sequences.add(currentSequence.toString());
            }

            br.close();
            fr.close();

            //ArrayList<ArrayList<String>> result = new ArrayList<>();
            HashMap<String, String> idsSequences = new HashMap<>();

            for (int i = 0; i < ids.size(); i++) {
                idsSequences.put(ids.get(i), sequences.get(i));
            }

            //result.add(ids);
            //result.add(sequences);

            return idsSequences;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }

    }
    */

    public static HashMap<String, Pair<String, String>> parse_prd(String prdPath) {

        try {

            FileReader fr = new FileReader(prdPath);
            BufferedReader br = new BufferedReader(fr);

            
            HashMap<String, Pair<String, String>> idsSequencesPredictions = new HashMap<>();

            String line;
            String seqId = "";
            String aminoAcidSequence = "";
            String prediction = "";

            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    //ids.add(line.substring(1).trim());
                    seqId = line.substring(1).trim();
                } else if (line.startsWith("AS ")) {
                    //sequences.add(line.substring(3).trim());
                    aminoAcidSequence = line.substring(3).trim();
                } else if (line.startsWith("PS ") || line.startsWith("SS ")) {
                    //predictions.add(line.substring(3).trim());
                    prediction = line.substring(3).trim();
                }

                if (seqId != "" && aminoAcidSequence != "" && prediction != "") {
                    idsSequencesPredictions.put(seqId, new Pair<String,String>(aminoAcidSequence, prediction));
                }
            }

            br.close();
            fr.close();

            //ArrayList<ArrayList<String>> result = new ArrayList<>();

            //result.add(ids);
            //result.add(sequences);
            //result.add(predictions);

            return idsSequencesPredictions;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

    public static Pair<ArrayList<String>, Pair<String, String>> pastaAln(String alignmentFile) {
        try {
            // read file
            Path alignmentFilePath = Paths.get(alignmentFile);
            String[] fileContent = Files.readString(alignmentFilePath).split("\n");

            String originalSequence = "";
            String seqId = "";

            // find AS sequence
            for (String line : fileContent) {
                if (line.startsWith(">")) {
                    seqId = line.substring(1).trim();
                } else if (line.startsWith("AS")) {
                    originalSequence = line.split(" ")[1].trim();
                }
            }

            if (originalSequence == "") {
                throw new IllegalArgumentException("Alignment File does not contain original Sequence! Please make sure your file contains a line starting with AS_ followed by the sequence");
            }

            // read all aligned sequences
            ArrayList<String> alignmentSequences = new ArrayList<>();
            alignmentSequences.add(originalSequence);
            for (String line : fileContent) {
                String[] parts = line.split(" ");
                if (isInteger(parts[0])) {
                    //parts[1] = parts[1].replaceAll("-", "X");
                    parts[1] = parts[1].replaceAll("[UZOB]", "X");
                    parts[1] = parts[1].trim();
                    alignmentSequences.add(parts[1]);
                }
            }
            

            return new Pair<>(alignmentSequences, new Pair<>(seqId, originalSequence));
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static boolean isInteger(String input) {
        for (int i = 0; i < input.length(); i++) {
            if (!Character.isDigit(input.charAt(i))) {
                return false;
            }
        }

        return true;
    }

    public static ArrayList<String> getAlignmentFiles(String directory) {
    try (Stream<Path> stream = Files.list(Paths.get(directory))) {
        return stream
            .filter(path -> path.toString().endsWith(".aln"))
            .map(path -> path.toAbsolutePath().toString())
            .collect(Collectors.toCollection(ArrayList::new));
    } catch (IOException e) {
        e.printStackTrace();
        return null;
    }
}

    public static int identifyModelFile(String modelFile) {
        try {
            Path modelFilePath = Paths.get(modelFile);
            String firstLine = Files.readString(modelFilePath).split("\n")[0];

            return switch (firstLine) {
                case "// Matrix3D" -> 1;
                case "// Matrix4D" -> 3;
                case "// Matrix6D" -> 4;
                default -> -1;
            };
        } catch (IOException e) {
            e.printStackTrace();
            return -1;
        }
    }
}
