package edu.umich.andykong.ptmshepherd.iterativelocalization;

import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Stream;

public class MatchedIonTable {
    public ArrayList<MatchedIonRow> rows;
    private double minProjVal;
    private double maxProjVal;

    public MatchedIonTable(){
        this.rows = new ArrayList<>();
    }

    public Stream<MatchedIonRow> stream() {
        return rows.stream();
    }

    public void addRow(float massError, float intensity, boolean isDecoy) {
        this.rows.add(new MatchedIonRow(massError, intensity, isDecoy));
    }

    public void mergeProjectedData(RealMatrix projectedData) {
        // Iterate over each row in the projected data and update the corresponding MatchedIonRow
        for (int i = 0; i < projectedData.getRowDimension(); i++) {
            // Assuming a 1-dimensional projection, get the value from the first column
            Double projectionValue = projectedData.getEntry(i, 0);
            this.rows.get(i).setProjectedDataRow(projectionValue);
        }
    }

    public void sortRowsByProjVal() {
        Collections.sort(this.rows);
        this.maxProjVal = this.rows.get(this.rows.size()-1).projVal;
        this.minProjVal = this.rows.get(0).projVal;
    }

    public void reverseProjVals() {
        for (MatchedIonRow row : this.rows)
            row.projVal *= -1;
    }

    public double getMaxProjVal() {
            return this.maxProjVal;
    }

    public double getMinProjVal() {
        return this.minProjVal;
    }

    public class MatchedIonRow implements Comparable<MatchedIonRow> {
        double massError;
        double intensity;
        boolean isDecoy;
        double projVal;

        public MatchedIonRow(float massError, float intensity, boolean isDecoy) {
            this.massError = Math.log10(massError+0.1);
            this.intensity = Math.log10(intensity);
            this.isDecoy = isDecoy;
        }

        public String toString() {
            return String.format("%f\t%f\t%f\t%b", this.massError, this.intensity, this.projVal, this.isDecoy);
        }

        public void setProjectedDataRow(double projVal) {
            this.projVal = projVal;
        }

        @Override
        public int compareTo(MatchedIonRow other) {
            return Double.compare(this.projVal, other.projVal);
        }
    }

}
