package main;


import static main.gor4_utils.prediction_utils.read_matrices_from_file_6d_4d;

import java.util.ArrayList;

import static main.validation_utils.val_utils.write_detailed_file;

public class testing {

    public static void main(String[] args) {

        String id = "test";
        String sequence = "ACDEFGHIKLMNPQRSTVWY";
        String filePath = "/home/jonas/Documents/Uni/ProPra/Data/gor_examples/cb513_gor4_cb513testtestettest.sum";

        ArrayList<String> ids = new ArrayList<>();
        ids.add(id);
        ArrayList<String> sequences = new ArrayList<>();
        sequences.add(sequence);
        ArrayList<String> predictions = new ArrayList<>();
        ArrayList<String> observed = new ArrayList<>();
        observed.add(sequence);
        write_detailed_file(filePath, ids, sequences, predictions, observed);

    }
}
