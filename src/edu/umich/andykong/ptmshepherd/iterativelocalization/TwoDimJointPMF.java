package edu.umich.andykong.ptmshepherd.iterativelocalization;

public class TwoDimJointPMF {
        int[][] grid;
        int nVals;

        boolean direction1;
        boolean direction2;

        public TwoDimJointPMF(int size1, int size2, boolean direction1, boolean direction2) {
                this.grid = new int[size1][size2];
                this.nVals = 0;
                this.direction1 = direction1;
                this.direction2 = direction2;
        }

        public void addVal(int x, int y) {
                this.grid[x][y]++;
                this.nVals++;
        }

        public void addVals(int[] xs, int[] ys) {
                for (int i = 0; i < xs.length; i++)
                        addVal(xs[i], ys[i]);
        }

        public int getXDim() {
                return this.grid.length;
        }

        public int getYDim() {
                return this.grid[0].length;
        }

        public int getnVals() {
                return this.nVals;
        }
}