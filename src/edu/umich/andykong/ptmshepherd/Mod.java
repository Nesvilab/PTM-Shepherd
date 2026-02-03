package edu.umich.andykong.ptmshepherd;

public class Mod {
    public int position;
    public double mass;

    public Mod(int position, double mass) {
        this.position = position;
        this.mass = mass;
    }
    
    public String toString() {
        return String.format("%d(%.4f)", position, mass);
    }

    public boolean equals(Mod other) {
        return this.position == other.position && this.mass == other.mass;
    }
}
