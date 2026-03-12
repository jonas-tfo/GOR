package main.utils;

@FunctionalInterface
public interface PostProcessingFunction {
    String apply(String rawPrediction, int windowSize);
}
