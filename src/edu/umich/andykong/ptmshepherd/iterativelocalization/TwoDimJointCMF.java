package edu.umich.andykong.ptmshepherd.iterativelocalization;

public class TwoDimJointCMF {
    double probGrid[][];
    int countGrid[][];
    int nVals;
    TwoDimJointPMF PMF;
    boolean direction1;
    boolean direction2;

    public TwoDimJointCMF(TwoDimJointPMF PMF, boolean direction1, boolean direction2) {
        this.probGrid = new double[PMF.getXDim()][PMF.getYDim()];
        this.countGrid = new int[PMF.getXDim()][PMF.getYDim()];
        this.PMF = PMF;
        this.nVals = PMF.nVals;
        this.direction1 = direction1;
        this.direction2 = direction2;

        calculateCMF();
    }

    private void calculateCMF() {
        // First two loops are to calculate this for every value
        for (int i1 = 0; i1 < getXDim(); i1++) {
            for (int j1 = 0; j1 < getYDim(); j1++) {
                int count = 0; // Number of vals to be integrated
                // These two loops are the values that will be integrated
                for (int i2 = 0; i2 <= i1; i2++) { // Sum values up to X (intensity)
                    for (int j2 = j1; j2 >= 0; j2--) { // Sum values great than Y (mass error)
                        count++;
                    }
                }
                this.countGrid[i1][j1] = count;
                this.probGrid[i1][j1] = (double) count / (double) this.nVals;
            }
        }
    }

    public int getXDim() {
        return this.countGrid.length;
    }

    public int getYDim() {
        return this.countGrid[0].length;
    }

    public String toString() {
        StringBuffer sb = new StringBuffer();

        sb.append("count grid");
        for (int j = 0; j < getXDim(); j++) {
            sb.append(j + ".\t");
        }
        sb.append("\n");
        for (int i = 0; i < getYDim(); i++) {
            sb.append(i + ".\t");
            for (int j = 0; j < getXDim(); j++)
                sb.append(String.format("%d\t", this.countGrid[j][i]));
            sb.append("\n");
        }
        sb.append("\t");
        sb.append("\n");

        sb.append("prob grid\n");
        for (int i = 0; i < getYDim(); i++) {
            sb.append(i + ".\t");
            for (int j = 0; j < getXDim(); j++)
                sb.append(String.format("%.2f\t", this.probGrid[j][i]));
            sb.append("\n");
        }
        sb.append("\t");
        for (int j = 0; j < getXDim(); j++) {
            sb.append(j + ".\t\t");
        }
        sb.append("\n");
        return sb.toString();
    }


    public static void main(String args[]) {
        TwoDimJointPMF tmpPMF = new TwoDimJointPMF(2,4,true,true);
        tmpPMF.addVal(0, 0);
        tmpPMF.addVal(0, 1);
        tmpPMF.addVal(0, 2);
        tmpPMF.addVal(0, 3);
        tmpPMF.addVal(1, 0);
        tmpPMF.addVal(1, 1);
        tmpPMF.addVal(1, 2);
        tmpPMF.addVal(1, 3);

        TwoDimJointCMF tmpCMF = new TwoDimJointCMF(tmpPMF, true, true);

        System.out.println(tmpCMF);

    }
}