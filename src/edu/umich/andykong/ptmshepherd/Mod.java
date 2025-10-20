package edu.umich.andykong.ptmshepherd;

public class Mod {
    public int position;
    public float mass;

    public Mod(int position, float mass) {
        this.position = position;
        this.mass = mass;
    }
    
    public String toString() {
        return String.format("%d(%.4f)", position, mass);
    }
}
