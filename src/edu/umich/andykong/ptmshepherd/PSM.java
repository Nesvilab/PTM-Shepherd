package edu.umich.andykong.ptmshepherd;

import edu.umich.andykong.ptmshepherd.core.Spectrum;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.reNormName;

/**
 * PSM class to hold parsed line info.
 */
public class PSM {
    public final int lineNum; // 0 indexed starting from header, 1 indexed starting from data
    public ArrayList<String> spLine;
    private final String spec;
    private String fileName;
    private final int specNum;
    private final String pep;
    private float [] modArr;
    private Float dMass;
    private float calcPepMass;
    private final int charge;
    private float originalDeltaMass;	// what was listed in the PSM table before analysis
    private TreeMap<Integer, Float> originalAssignedMods;
    private TreeMap<Integer, Float> assignedMods;		// position -> mass

    PSM(int lineNum, String line, int massdiffToVarmod, int specCol, int pepCol, int chargeCol, int calcMassCol, int dMassCol, int assignedModCol, int msfraggerLocalizationCol) {
        this.lineNum = lineNum;
        this.spLine = new ArrayList<>(Arrays.asList(line.replace("\n","").split("\t")));
        this.fileName = null;
        this.spec = reNormName(spLine.get(specCol));
        String[] spSpec = spec.split("\\.", -1);
        specNum = Integer.parseInt(spSpec[spSpec.length-2]);
        this.pep = spLine.get(pepCol);
        this.charge = Integer.parseInt(spLine.get(chargeCol));
        this.calcPepMass = Float.parseFloat(spLine.get(calcMassCol));
        initializeMods(massdiffToVarmod, dMassCol, assignedModCol, msfraggerLocalizationCol);
    }

    public String printLine() {
        return String.join("\t", spLine);
    }

    public void addValAtColumn(int colIdx, String val){
        spLine.add(colIdx, val);
    }

    public void replaceValAtColumn(int colIdx, String val){
        spLine.set(colIdx, val);
    }

    public String getSpec() {
        return spec;
    }

    public int getScanNum() {
        return specNum;
    }

    public String getFileName() {
        if (fileName == null)
            fileName = spec.substring(0, spec.indexOf("."));
        return fileName;
    }

    public String getPep() {
        return pep;
    }

    public int getCharge() {
        return charge;
    }

    public Float getCalcPepmass() {
        return calcPepMass;
    }

    /**
     * Initialize delta mass and assigned mods, accounting for mass-diff-to-varmod setting
     * from MSFragger and any previous modifications to the PSM table (e.g., if this is a re-run)
     * @param massdiffToVarmod MSFragger setting for mass diff to varmod (0 = none, 1 = remove, 2 = keep)
     */
    public void initializeMods(int massdiffToVarmod, int dMassCol, int assignedModCol, int msfraggerLocalizationCol) {
        originalDeltaMass = Float.parseFloat(spLine.get(dMassCol));
        originalAssignedMods = initAssignedMods(assignedModCol);

        if (massdiffToVarmod == 0) {
            // no mass diff to varmod, use original delta mass and assigned mods
            dMass = originalDeltaMass;
            assignedMods = originalAssignedMods;
        } else {
            String msfraggerLocStr = spLine.get(msfraggerLocalizationCol);
            if (msfraggerLocStr.isEmpty()) {
                // unmodified PSM
                dMass = originalDeltaMass;
                assignedMods = originalAssignedMods;
            } else {
                // find the variable mod that was the delta mass using the MSFragger localization info
                int deltaMassPos = 0;
                for (int i = 0; i < msfraggerLocStr.length(); i++) {
                    if (Character.isLowerCase(msfraggerLocStr.charAt(i))) {
                        deltaMassPos = i + 1;
                        break;
                    }
                }

                // remove the delta mass from the assigned mods and add its mass to the dMass for analysis
                assignedMods = new TreeMap<>();
                boolean foundDeltaMod = false;
                for (Map.Entry<Integer, Float> mod : originalAssignedMods.entrySet()) {
                    if (mod.getKey() != deltaMassPos) {
                        assignedMods.put(mod.getKey(), mod.getValue());
                    } else {
                        foundDeltaMod = true;
                    }
                }

                // if delta mass was removed, add it back (but do not if it was kept)
                if (massdiffToVarmod == 1 && foundDeltaMod) {
                    dMass = originalDeltaMass + originalAssignedMods.get(deltaMassPos);
                } else {
                    dMass = originalDeltaMass;
                }
            }
        }
    }

    public TreeMap<Integer, Float> initAssignedMods(int assignedModCol) {
        TreeMap<Integer, Float> mods = new TreeMap<>();
        String strMods = spLine.get(assignedModCol);
        if (!strMods.isEmpty()) {
            String[] spMods = strMods.split(",", -1);
            for (String spMod : spMods) {
                int p = spMod.indexOf("(");
                int q = spMod.indexOf(")");
                String spos = spMod.substring(0, p).trim();
                float mass = Float.parseFloat(spMod.substring(p + 1, q).trim());
                int pos;
                if (spos.equals("N-term"))
                    pos = 0;
                else if (spos.equals("c"))
                    pos = this.getPep().length();
                else
                    pos = Integer.parseInt(spos.substring(0, spos.length() - 1));
                mods.put(pos, mass);
            }
        }
        return mods;
    }

    public float [] getModsAsArray() {
        if (modArr == null) {
            modArr = new float[getPep().length()];
            Arrays.fill(modArr, 0.0f);
            for (Map.Entry<Integer, Float> mod : assignedMods.entrySet()) {
                if (mod.getKey() == 0)
                    modArr[0] = mod.getValue();
                else
                    modArr[mod.getKey()-1] = mod.getValue();
            }
        }
        return modArr;
    }

    // use this to get the delta mass for actual analyses
    public float getDMass() {
        return dMass;
    }

    public TreeMap<Integer, Float> getAssignedMods() {
        return assignedMods;
    }

    public String printAssignedMods() {
        ArrayList<String> modStrs = new ArrayList<>();
        for (Map.Entry<Integer, Float> mod : assignedMods.entrySet()) {
            StringBuilder sb = new StringBuilder();
            if (mod.getKey() == 0) {
                sb.append("N-term");
            } else if (mod.getKey() == getPep().length()) {
                sb.append("C-term");
            } else {
                sb.append(mod.getKey()).append(getPep().charAt(mod.getKey()-1));
            }
            sb.append(String.format("(%.4f)", mod.getValue()));
            modStrs.add(sb.toString());
        }
        return String.join(",", modStrs);
    }

    public String getColumnValue(int colIndex) {
        return spLine.get(colIndex);
    }

    public ArrayList<String> getSpLine() {
        return spLine;
    }

    /**
     * Update the delta mass for this PSM (e.g., after glycan analysis), including updating the
     * calculated peptide mz and mass columns. If the delta mass was previously placed as a variable mod,
     * reset the delta mass to the original value before updating with the new delta mass.
     * @param newDeltaMass new delta mass
     */
    public void updateDeltaMass(float newDeltaMass, int massdiffToVarmod, float prevTheoreticalMass, int peptideCalcMassCol, int calcMZcol, int dmassCol, int assignedModCol) {
        double prevCalcPeptideMass = Double.parseDouble(spLine.get(peptideCalcMassCol));
        if (massdiffToVarmod == 1) {
            float originalDeltaMass = dMass + prevTheoreticalMass;
            dMass = originalDeltaMass - newDeltaMass;		// update delta mass in case the glycan composition has changed
            // Remove previous mass prior to subtracting the new delta mass (in case the new delta mass is different)
            double correctedMass = prevCalcPeptideMass - prevTheoreticalMass;
            spLine.set(peptideCalcMassCol, String.format("%.4f", correctedMass + newDeltaMass));
            spLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (correctedMass + newDeltaMass), getCharge())));
        } else {
            // original delta mass was left intact, simply subtract the glycan mass
            dMass = dMass - newDeltaMass;
            spLine.set(peptideCalcMassCol, String.format("%.4f", prevCalcPeptideMass + newDeltaMass));
            spLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ((float) (prevCalcPeptideMass + newDeltaMass), getCharge())));
        }
        updatePSMline(dmassCol, assignedModCol);
    }

    // Because delta mass and assigned mods can change (e.g., in glyco), this method updates the PSM line
    public void updatePSMline(int dMassCol, int assignedModCol) {
        spLine.set(dMassCol, String.format("%.4f", dMass));
        spLine.set(assignedModCol, printAssignedMods());
    }

    public String toString() {
        return String.join("\t", this.spLine);
    }

}