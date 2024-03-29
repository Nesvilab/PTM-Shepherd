/*
 *    Copyright 2022 University of Michigan
 *
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *
 *        http://www.apache.org/licenses/LICENSE-2.0
 *
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

package edu.umich.andykong.ptmshepherd.diagnosticmining;

import com.google.common.util.concurrent.AtomicDouble;
import edu.umich.andykong.ptmshepherd.PTMShepherd;

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
    long nMod;
    long nUnmod;

    // For immonium and Y ions
    AtomicInteger nTotal;
    AtomicInteger nWithIon;
    AtomicDouble wIonInt; // running sum
    AtomicDouble nWithIonZb;
    AtomicDouble wIonIntZb;
    AtomicInteger zSum; // running sum of charge state to be averaged later

    // For remainder masses
    AtomicInteger nShiftedIons;
    AtomicInteger nUnshiftedIons;
    AtomicDouble pctCoverage; // running sum of shifted ions
    AtomicDouble pctCoverageUnmod; // running sum of unshifted ions

    public DiagnosticProfileRecord(double peakApex, String modName, String type, double mass, double adjustedMass, double q,
                                    double rbc, double propWIonTreat, double propWIonCont,
                                    double propWIonIntensityTreat, double propWIonIntensityCont,
                                   long nMod, long nUnmod) {
        this.peakApex = peakApex;
        this.modName = modName;
        this.type = type;
        this.mass = mass;
        this.adjustedMass = adjustedMass;
        this.q = q;
        this.rbc = rbc;
        this.nMod = nMod;
        this.nUnmod = nUnmod;

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

        this.zSum = new AtomicInteger();
    }

    public String toString(int printAuc) {
        String newLine;

        // Optionally remove isotopic peaks
        boolean printIsos = Boolean.parseBoolean(PTMShepherd.getParam("diagmine_printIsotopes"));
        if (!printIsos) {
            if (this.modName.contains("isotopic peak"))
                return "";
        }

        if (printAuc == 1) {
            if (this.type.equals("diagnostic")) {
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "\t\t" +
                                "%.02f\t%.02f\t" +
                                "%.04f\t" +
                                "%.02f\t%.02f\t" +
                                "%f\t%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        //propWIonSpectrumLevel, wIonIntensity, //spectrum level stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        (float) this.zSum.get() / (float) this.nWithIon.get(), //average charge
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        //this.q, this.rbc); //selection stats
                        foldChange, this.rbc);
            } else if (this.type.equals("peptide")) {
                double remainderDelta = -1.0 * (this.peakApex - this.adjustedMass);
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "%.04f\t\t" +
                                "%.02f\t%.02f\t" +
                                "\t" +
                                "%.02f\t%.02f\t" +
                                "%f\t%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        remainderDelta, //lost mass stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        foldChange, this.rbc);
            } else {
                float remainderOdds = (float) this.pctCoverage.get() / (float) this.nTotal.get();
                if (remainderOdds < Double.parseDouble(PTMShepherd.getParam("diagmine_fragMinPropensity")) / 100.0) //todo this should be calculated and filtered elsewhere
                    return "";
                double remainderDelta = -1.0 * (this.peakApex - this.adjustedMass);
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "%.04f\t%.02f\t" +
                                "%.02f\t%.02f\t" +
                                "\t" +
                                "%.02f\t%.02f\t" +
                                "%f\t%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        remainderDelta, remainderOdds * 100.0, //ion level stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        foldChange, this.rbc);
            }
            return newLine;
        }
        else {
            if (this.type.equals("diagnostic")) {
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "\t\t" +
                                "%.02f\t%.02f\t" +
                                "%.04f\t" +
                                "%.02f\t%.02f\t" +
                                "%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        //propWIonSpectrumLevel, wIonIntensity, //spectrum level stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        (float) this.zSum.get() / (float) this.nWithIon.get(), //average charge
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        //this.q, this.rbc); //selection stats
                        foldChange);
            } else if (this.type.equals("peptide")) {
                double remainderDelta = -1.0 * (this.peakApex - this.adjustedMass);
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "%.04f\t\t" +
                                "%.02f\t%.02f\t" +
                                "\t" +
                                "%.02f\t%.02f\t" +
                                "%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        remainderDelta, //lost mass stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        foldChange);
            } else {
                float remainderOdds = (float) this.pctCoverage.get() / (float) this.nTotal.get();
                if (remainderOdds < Double.parseDouble(PTMShepherd.getParam("diagmine_fragMinPropensity")) / 100.0) //todo this should be calculated and filtered elsewhere
                    return "";
                double remainderDelta = -1.0 * (this.peakApex - this.adjustedMass);
                double foldChange = (this.propWIonTreatIonLevel * this.wIonIntTreatIonLevel) /
                        (this.propWIonControlIonLevel * this.wIonIntContIonLevel);
                foldChange = (foldChange > 100.0) ? 100.0 : foldChange;
                newLine = String.format("%.04f\t%s\t%s\t%.04f\t" +
                                "%.04f\t%.02f\t" +
                                "%.02f\t%.02f\t" +
                                "\t" +
                                "%.02f\t%.02f\t" +
                                "%f\n",
                        this.peakApex, this.modName, this.type, this.adjustedMass, //basic stats
                        remainderDelta, remainderOdds * 100.0, //ion level stats
                        this.propWIonTreatIonLevel * 100.0, this.propWIonControlIonLevel * 100.0, //ion level stats for propensity
                        this.wIonIntTreatIonLevel, this.wIonIntContIonLevel, //ion level stats for intensity
                        foldChange);
            }
            return newLine;
        }
    }
}
