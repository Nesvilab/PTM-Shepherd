package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.common.util.concurrent.AtomicDouble;
import org.checkerframework.checker.units.qual.A;

import java.util.concurrent.atomic.AtomicInteger;

public class DiagnosticProfileRecord {
    double peakApex;
    String type;
    double mass;
    double adjustedMass;
    double q;
    double rbc;
    double propWIonTreatIonLevel;
    double propWIonControlIonLevel;
    double wIonIntTreatIonLevel;
    double wIonIntContIonLevel;

    // For immonium and Y ions
    AtomicInteger nTotal;
    AtomicInteger nWithIon;
    AtomicDouble wIonInt; // running sum
    AtomicDouble nWithIonZb;
    AtomicDouble wIonIntZb;

    // For remainder masses
    AtomicInteger nShiftedIons;
    AtomicInteger nUnshiftedIons;
    AtomicDouble pctCoverage; // running sum of shifted ions
    AtomicDouble pctCoverageUnmod; // running sum of unshifted ions
    //TODO add intensities?

    public DiagnosticProfileRecord (double peakApex, String type, double mass, double adjustedMass, double q,
                                    double rbc, double propWIonTreat, double propWIonCont,
                                    double propWIonIntensityTreat, double propWIonIntensityCont) {
        this.peakApex = peakApex;
        this.type = type;
        this.mass = mass;
        this.adjustedMass = adjustedMass;
        this.q = q;
        this.rbc = rbc;

        this.propWIonTreatIonLevel = propWIonTreat;
        this.propWIonControlIonLevel = propWIonCont;
        this.wIonIntTreatIonLevel = propWIonIntensityTreat;
        this.wIonIntContIonLevel = propWIonIntensityCont;

        this.nTotal = new AtomicInteger();
        this.nWithIon = new AtomicInteger();
        this.wIonInt = new AtomicDouble();

        this.nShiftedIons = new AtomicInteger();
        this.nUnshiftedIons = new AtomicInteger();
        this.pctCoverage = new AtomicDouble();
        this.pctCoverageUnmod = new AtomicDouble();
    }

    public String toString() {
        String newLine;
        if (this.type.equals("imm") || this.type.equals("Y")) {
            float propWIonSpectrumLevel = this.nWithIon.get() / (float) this.nTotal.get();
            float wIonIntensity = (float) (this.wIonInt.get() / this.nWithIon.get());
            newLine = String.format("%.04f\t%s\t%.04f\t%.04f\t" +
                            "%.02f\t%.02f\t\t" +
                            "%.02f\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%e\t%f\n",
                    this.peakApex, this.type, this.mass, this.adjustedMass, //basic stats
                    propWIonSpectrumLevel, wIonIntensity, //spectrum level stats
                    this.propWIonControlIonLevel, this.propWIonTreatIonLevel, //ion level stats for propensity
                    this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                    this.q, this.rbc); //selection stats
        } else {
            float remainderOdds = (float) this.nShiftedIons.get() / (float) this.nUnshiftedIons.get();
            newLine = String.format("%.04f\t%s\t%.04f\t%.04f\t" +
                            "\t\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%e\t%f\n",
                    this.peakApex, this.type, this.mass, this.adjustedMass, //basic stats
                    remainderOdds, //spectrum level stats
                    this.propWIonControlIonLevel, this.propWIonTreatIonLevel, //ion level stats for propensity
                    this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                    this.q, this.rbc);
        }
        return newLine;
    }
}
