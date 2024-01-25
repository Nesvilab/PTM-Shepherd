package edu.umich.andykong.ptmshepherd.iterativelocalization;

import java.util.PriorityQueue;

public class TwoDimJointQValue {
    double[][] qValGrid;
    boolean[][] gridMask;

    TwoDimJointPMF targetPMF;
    TwoDimJointPMF decoyPMF;

    TwoDimJointQValue(TwoDimJointPMF targetPMF, TwoDimJointPMF decoyPMF) {
        this.qValGrid = new double[targetPMF.getXDim()][targetPMF.getYDim()];
        this.gridMask = new boolean[targetPMF.getXDim()][targetPMF.getYDim()];

        this.targetPMF = targetPMF;
        this.decoyPMF = decoyPMF;

        calculateQValGrid();
        propagateMin();
    }

    private void calculateQValGrid() {
        // First two loops are to calculate this for every value
        for (int i1 = 0; i1 < getXDim(); i1++) {
            for (int j1 = 0; j1 < getYDim(); j1++) {
                int targetCount = 0; // Number of target vals to be integrated
                int decoyCount = 1; // Number of decoy vals to be integrated, unbiased estimator
                // These two loops are the values that will be integrated
                for (int i2 = 0; i2 <= i1; i2++) { // Sum values up to X (intensity)
                    for (int j2 = j1; j2 >= 0; j2--) { // Sum values great than Y (mass error)
                        targetCount += this.targetPMF.grid[i2][j2];
                        decoyCount += this.decoyPMF.grid[i2][j2];
                    }
                }
                this.qValGrid[i1][j1] = (double) decoyCount / (double) Math.max(targetCount, 1);
            }
        }
    }

    class Cell implements Comparable<Cell> {
        int x, y;
        double value;

        Cell(int x, int y, double value) {
            this.x = x;
            this.y = y;
            this.value = value;
        }

        @Override
        public int compareTo(Cell other) {
            return Double.compare(this.value, other.value);
        }
    }

    public void propagateMin() {
        PriorityQueue<Cell> pq = new PriorityQueue<>();
        for (int i = 0; i < getXDim(); i++) {
            for (int j = 0; j < getYDim(); j++) {
                pq.offer(new Cell(i, j, this.qValGrid[i][j]));
            }
        }

        while (!pq.isEmpty()) {
            Cell cell = pq.poll();

            if (this.gridMask[cell.x][cell.y]) {
                continue;
            }

            double minValue = cell.value;

            for (int i = cell.x; i < getXDim(); i++) {
                for (int j = cell.y; j >= 0; j--) {
                    if (!this.gridMask[i][j]) {
                        this.qValGrid[i][j] = Math.min(this.qValGrid[i][j], minValue);
                        this.gridMask[i][j] = true;
                    }
                }
            }
        }
    }

    public double getQVal(int x, int y) {
        return this.qValGrid[x][y];
    }

    public int getXDim() {
        return this.qValGrid.length;
    }

    public int getYDim() {
        return this.qValGrid[0].length;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();

        sb.append("decoy grid\t");
        for (int j = 0; j < getXDim(); j++) {
            sb.append(j + "\t");
        }
        sb.append("\n");
        for (int i = 0; i < getYDim(); i++) {
            sb.append(i + ".\t");
            for (int j = 0; j < getXDim(); j++)
                sb.append(String.format("%d\t", this.decoyPMF.grid[j][i]));
            sb.append("\n");
        }
        sb.append("\t");
        sb.append("\n");

        sb.append("target grid\t");
        for (int j = 0; j < getXDim(); j++) {
            sb.append(j + "\t");
        }
        sb.append("\n");
        for (int i = 0; i < getYDim(); i++) {
            sb.append(i + ".\t");
            for (int j = 0; j < getXDim(); j++)
                sb.append(String.format("%d\t", this.targetPMF.grid[j][i]));
            sb.append("\n");
        }
        sb.append("\t");
        sb.append("\n");

        sb.append("q-val grid\t");
        for (int j = 0; j < getXDim(); j++) {
            sb.append(j + "\t");
        }
        sb.append("\n");
        for (int i = 0; i < getYDim(); i++) {
            sb.append(i + ".\t");
            for (int j = 0; j < getXDim(); j++)
                sb.append(String.format("%.2f\t", this.qValGrid[j][i]));
            sb.append("\n");
        }
        sb.append("\t");
        sb.append("\n");
        return sb.toString();
    }

    public static void main(String args[]) {
        TwoDimJointPMF targetPMF = new TwoDimJointPMF(2,4,true,true);
        TwoDimJointPMF decoyPMF = new TwoDimJointPMF(2,4,true,true);
        targetPMF.addVal(0, 0);
        targetPMF.addVal(0, 1);
        targetPMF.addVal(0, 2);
        targetPMF.addVal(0, 3);
        targetPMF.addVal(1, 0);
        targetPMF.addVal(1, 1);
        targetPMF.addVal(1, 2);
        targetPMF.addVal(1, 3);

        decoyPMF.addVal(0, 3);
        decoyPMF.addVal(0, 2);

        TwoDimJointQValue tmpQValGrid = new TwoDimJointQValue(targetPMF, decoyPMF);

        System.out.println(tmpQValGrid);

    }
}
