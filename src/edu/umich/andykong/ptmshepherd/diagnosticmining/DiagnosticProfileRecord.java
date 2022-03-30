package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.common.util.concurrent.AtomicDouble;
import org.checkerframework.checker.units.qual.A;

import java.util.concurrent.atomic.AtomicInteger;

public class DiagnosticProfileRecord {
    double peakApex;
    String modName;
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

    public DiagnosticProfileRecord (double peakApex, String modName, String type, double mass, double adjustedMass, double q,
                                    double rbc, double propWIonTreat, double propWIonCont,
                                    double propWIonIntensityTreat, double propWIonIntensityCont) {
        this.peakApex = peakApex;
        this.modName = modName;
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
        if (this.type.equals("diagnostic")) {
            float propWIonSpectrumLevel = this.nWithIon.get() / (float) this.nTotal.get();
            float wIonIntensity = (float) (this.wIonInt.get() / this.nWithIon.get());
            newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                            "\t\t" +
                            "%.02f\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%f\n",
                    this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                    //propWIonSpectrumLevel, wIonIntensity, //spectrum level stats
                    this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                    this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                    //this.q, this.rbc); //selection stats
                    this.rbc);
        } else if (this.type.equals("peptide")) {
            float propWIonSpectrumLevel = this.nWithIon.get() / (float) this.nTotal.get();
            float wIonIntensity = (float) (this.wIonInt.get() / this.nWithIon.get());
            float remainderDelta = (float) this.peakApex - (float) this.adjustedMass;
            newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                            "\t%.04f\t" +
                            "%.02f\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%f\n",
                    this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                    remainderDelta, //lost mass stats
                    //propWIonSpectrumLevel, wIonIntensity, //spectrum level stats
                    this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                    this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                    //this.q, this.rbc); //selection stats
                    this.rbc);
        } else {
            float remainderOdds = (float) this.pctCoverage.get() / (float) this.nTotal.get();
            double remainderDelta = this.peakApex - this.adjustedMass;
            newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                            "%.02f\t%.04f\t" +
                            "%.02f\t%.02f\t" +
                            "%.02f\t%.02f\t" +
                            "%f\n",
                    this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                    remainderOdds, remainderDelta, //spectrum level stats
                    this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                    this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                    //this.q, this.rbc);
                    this.rbc);
        }
        return newLine;
    }
}
