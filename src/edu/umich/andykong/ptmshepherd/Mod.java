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
    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (o == null || getClass() != o.getClass())
            return false;
        Mod other = (Mod) o;
        return this.position == other.position && Double.compare(this.mass, other.mass) == 0;
    }

    @Override
    public int hashCode() {
        int result = Integer.hashCode(position);
        result = 31 * result + Double.hashCode(mass);
        return result;
    }
}
