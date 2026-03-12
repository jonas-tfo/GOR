package main.utils;

import java.util.ArrayList;

public class Math2 {

    /* standard deviation */
    public static double standard_deviation(ArrayList<Double> nums) {

        int n = nums.size();
        double standardDeviation = 0.0;

        double sum = 0.0;

        for (int i = 0; i < n; i++) {
            sum = sum + nums.get(i);
        }

        double mean = sum / n;

        for (int i = 0; i < n; i++) {
            standardDeviation += Math.pow((nums.get(i) - mean), 2);
        }
        double sq = standardDeviation / n;
        double res = Math.sqrt(sq);

        if (Double.isNaN(res)) {
            return 0.0;
        } else {
            return res;
        }
    }


    public static double mean(ArrayList<Double> nums) {
        int n = nums.size();
        double sum = 0.0;

        for (int i = 0; i < n; i++) {
            sum = sum + nums.get(i);
        }

        double res = sum / n;
        if (Double.isNaN(res)) {
            return 0.0;
        } else {
            return res;
        }
    }

    public static double max(ArrayList<Double> nums) {
        double max = Double.NEGATIVE_INFINITY;

        for (double num : nums) {
            if (num > max) {
                max = num;
            }
        }

        return max;
    }

    public static double min(ArrayList<Double> nums) {
        double min = Double.POSITIVE_INFINITY;

        for (double num : nums) {
            if (num < min) {
                min = num;
            }
        }

        return min;
    }

    public static double median(ArrayList<Double> nums) {
        int n = nums.size();

        ArrayList<Double> sortedNums = new ArrayList<>(nums);
        sortedNums.sort(Double::compareTo);

        if (n % 2 == 0) {
            return (sortedNums.get(n / 2 - 1) + sortedNums.get(n / 2)) / 2.0;
        } else {
            return sortedNums.get(n / 2);
        }
    }

    public static double quantile(ArrayList<Double> nums, double q) {

        int n = nums.size();

        ArrayList<Double> sortedNums = new ArrayList<>(nums);
        sortedNums.sort(Double::compareTo);

        double pos = q * (n + 1);
        int index = (int) pos;

        if (index < 1) {
            return sortedNums.get(0);
        } else if (index >= n) {
            double res = sortedNums.get(n - 1);
            if (Double.isNaN(res)) {
                return 0.0;
            } else {
                return res;
            }
        } else {
            double fraction = pos - index;
            double res = sortedNums.get(index - 1) + fraction * (sortedNums.get(index) - sortedNums.get(index - 1));
            if (Double.isNaN(res)) {
                return 0.0;
            } else {
                return res;
            }
        }
    }





}
