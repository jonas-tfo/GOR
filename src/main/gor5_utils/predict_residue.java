package main.gor5_utils;

public interface predict_residue<T> {
  double[] predictResidue(String window, T scoreMatrices);
}
