package edu.umich.andykong.ptmshepherd;

import edu.umich.andykong.ptmshepherd.core.AAMasses;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.glyco.GlycanAssignmentResult;

import java.util.*;
import java.util.stream.Collectors;

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
    private final String peptide;
    private float [] modArr;
    private Float dMass;
    private float calcPepMass;
    private float originalCalcPepMass;
    private final int charge;
    private float originalDeltaMass;	// what was listed in the PSM table before analysis
    private ArrayList<Mod> originalAssignedMods;
    private ArrayList<Mod> assignedMods;		// position -> mass
    private String modifiedPeptide;
    private String originalModifiedPeptide;
    public GlycanAssignmentResult glycanAssignmentResult = null; // result of glycan assignment, if any

    PSM(int lineNum, String line, int massdiffToVarmod, int specCol, int pepCol, int modPepCol, int chargeCol, int calcMassCol, int dMassCol, int assignedModCol, int msfraggerLocalizationCol) {
        this.lineNum = lineNum;
        this.spLine = new ArrayList<>(Arrays.asList(line.replace("\n","").split("\t")));
        this.fileName = null;
        this.spec = reNormName(spLine.get(specCol));
        String[] spSpec = spec.split("\\.", -1);
        specNum = Integer.parseInt(spSpec[spSpec.length-2]);
        this.peptide = spLine.get(pepCol);
        this.charge = Integer.parseInt(spLine.get(chargeCol));
        initializeMods(massdiffToVarmod, dMassCol, assignedModCol, modPepCol, msfraggerLocalizationCol, calcMassCol);
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

    public String getPeptide() {
        return peptide;
    }

    public int getCharge() {
        return charge;
    }

    public Float getCalcPepmass() {
        return calcPepMass;
    }

    public float getOriginalDeltaMass() {
        return originalDeltaMass;
    }
    public String getModifiedPeptide() {
        return modifiedPeptide;
    }
    public String getOriginalModifiedPeptide() {
        return originalModifiedPeptide;
    }
    public float getCalcPepMass() {
        return calcPepMass;
    }
    public float getOriginalCalcPepMass() {
        return originalCalcPepMass;
    }

    /**
     * Initialize delta mass and assigned mods, accounting for mass-diff-to-varmod setting
     * from MSFragger and any previous modifications to the PSM table (e.g., if this is a re-run)
     * @param massdiffToVarmod MSFragger setting for mass diff to varmod (0 = none, 1 = remove, 2 = keep)
     */
    public void initializeMods(int massdiffToVarmod, int dMassCol, int assignedModCol, int modPepCol, int msfraggerLocalizationCol, int calcMassCol) {
        originalDeltaMass = Float.parseFloat(spLine.get(dMassCol));
        originalAssignedMods = initAssignedMods(assignedModCol);
        originalModifiedPeptide = spLine.get(modPepCol).isEmpty() ? peptide : spLine.get(modPepCol);    // use peptide if no modified peptide
        originalCalcPepMass = Float.parseFloat(spLine.get(calcMassCol));

        if (massdiffToVarmod == 0) {
            // no mass diff to varmod, use original delta mass and assigned mods
            dMass = originalDeltaMass;
            assignedMods = originalAssignedMods;
            modifiedPeptide = originalModifiedPeptide;
            calcPepMass = originalCalcPepMass;
        } else {
            String msfraggerLocStr = spLine.get(msfraggerLocalizationCol);
            if (msfraggerLocStr.isEmpty()) {
                // unmodified PSM
                dMass = originalDeltaMass;
                assignedMods = originalAssignedMods;
                modifiedPeptide = originalModifiedPeptide;
                calcPepMass = originalCalcPepMass;
            } else {
                // find the variable mod that was the delta mass using the MSFragger localization info
                int deltaMassPos = 0;
                for (int i = 0; i < msfraggerLocStr.length(); i++) {
                    if (Character.isLowerCase(msfraggerLocStr.charAt(i))) {
                        deltaMassPos = i + 1;
                        break;
                    }
                }
                deltaMassPos = handleTerminalDeltaMass(deltaMassPos);

                // remove the delta mass from the assigned mods and add its mass to the dMass for analysis
                assignedMods = new ArrayList<>();
                boolean foundDeltaMod = false;
                float modMassAtDeltaPos = 0.0f;
                for (Mod mod : originalAssignedMods) {
                    if (mod.position != deltaMassPos) {
                        assignedMods.add(mod);
                    } else {
                        // this is the delta mod. Do not include it in the new assigned mods
                        foundDeltaMod = true;
                        modMassAtDeltaPos += mod.mass;  // in case there are multiple (e.g., fixed + offset at same site)
                    }
                }

                // remove the delta mass from the modified peptide if it is present as an assigned mod
                if (foundDeltaMod) {
                    modifiedPeptide = removeModFromModifiedPeptide(originalModifiedPeptide, deltaMassPos);
                } else {
                    modifiedPeptide = originalModifiedPeptide;
                }

                // if delta mass was removed, add it back (but do not if it was kept)
                if (massdiffToVarmod == 1 && foundDeltaMod) {
                    dMass = originalDeltaMass + modMassAtDeltaPos;
                    calcPepMass = originalCalcPepMass - modMassAtDeltaPos;
                } else {
                    dMass = originalDeltaMass;
                    calcPepMass = originalCalcPepMass;
                }
            }
        }
    }

    /**
     * Handle the case where the delta mass is at the N- or C-terminus and has been placed by MSFragger. This means
     * the MSFragger localization string will show it as the first or last side chain, but the assigned modifications
     * string will show it as a terminal modification.
     * Note: if there is both a terminal mod and a side chain mod at the first or last AA, this assumes the terminal
     * mod is the delta mass rather than the side chain mod.
     * todo: (ideally) handle case if there is both an offset & a variable/fixed mod at first or last AA (one terminal and one on side chain). Requires reading pepxml or propagating more info to psm.tsv
     * @param deltaMassPos initial deltaMassPos from MSFragger localization string
     * @return updated deltaMassPos
     */
    private int handleTerminalDeltaMass(int deltaMassPos) {
        if (deltaMassPos == 1) {
            for (Mod mod : originalAssignedMods) {
                if (mod.position == 0) {
                    deltaMassPos = 0;
                    break;
                }
            }
        }
        if (deltaMassPos == peptide.length()) {
            for (Mod mod : originalAssignedMods) {
                if (mod.position == peptide.length() + 1) {
                    deltaMassPos = peptide.length() + 1;
                    break;
                }
            }
        }
        return deltaMassPos;
    }

    public ArrayList<Mod> initAssignedMods(int assignedModCol) {
        ArrayList<Mod> mods = new ArrayList<>();
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
                else if (spos.equals("C-term"))
                    pos = this.getPeptide().length() + 1;   // record C-term as last residue + 1
                else
                    pos = Integer.parseInt(spos.substring(0, spos.length() - 1));
                mods.add(new Mod(pos, mass));
            }
        }
        return mods;
    }

    public float [] getModsAsArray() {
        if (modArr == null) {
            modArr = new float[getPeptide().length()];
            Arrays.fill(modArr, 0.0f);
            for (Mod mod : assignedMods) {
                if (mod.position == 0)
                    modArr[0] = mod.mass;
                else if (mod.position == getPeptide().length() + 1)
                    modArr[modArr.length - 1] = mod.mass; // C-term is last residue + 1, so last index
                else
                    modArr[mod.position-1] = mod.mass;
            }
        }
        return modArr;
    }

    // use this to get the delta mass for actual analyses
    public float getDMass() {
        return dMass;
    }

    public ArrayList<Mod> getAssignedMods() {
        return assignedMods;
    }

    public ArrayList<Mod> getOriginalAssignedMods() {
        return originalAssignedMods;
    }

    public String printAssignedMods() {
        ArrayList<String> modStrs = new ArrayList<>();
        for (Mod mod : assignedMods) {
            StringBuilder sb = new StringBuilder();
            if (mod.position == 0) {
                sb.append("N-term");
            } else if (mod.position == getPeptide().length() + 1) {
                sb.append("C-term");
            } else {
                sb.append(mod.position).append(getPeptide().charAt(mod.position-1));
            }
            sb.append(String.format("(%.4f)", mod.mass));
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
        double prevCalcPeptideMass = calcPepMass;
        if (massdiffToVarmod == 1) {
            float originalDeltaMass = dMass + prevTheoreticalMass;
            dMass = originalDeltaMass - newDeltaMass;		// update delta mass in case the glycan composition has changed
            // Remove previous mass prior to subtracting the new delta mass (in case the new delta mass is different)
            double correctedMass = prevCalcPeptideMass - prevTheoreticalMass;
            calcPepMass = (float) (correctedMass + newDeltaMass);
        } else {
            // original delta mass was left intact, simply subtract the glycan mass
            dMass = dMass - newDeltaMass;
            calcPepMass = (float) (prevCalcPeptideMass + newDeltaMass);
        }

        // update the spLine
        spLine.set(peptideCalcMassCol, String.format("%.4f", calcPepMass));
        spLine.set(calcMZcol, String.format("%.4f", Spectrum.neutralMassToMZ(calcPepMass, getCharge())));
        spLine.set(dmassCol, String.format("%.4f", dMass));
        spLine.set(assignedModCol, printAssignedMods());
    }

    /**
     * Add a modification to a modified peptide string. If massdiff to varmod was enabled, check if a previous
     * mod needs to be removed first and do so if necessary.
     */
    public void editModifiedPeptide(int modLocation, double modMass, int massdiffToVarmod, int modPeptideCol) {
        if (modLocation == 0 || modLocation == peptide.length() + 1) {
            return;   // do not edit terminal mods because glycan assignment will never adjust them
        }
        if (massdiffToVarmod != 0) {
            int adjustedModIndex = getModPepIndex(modifiedPeptide, modLocation);
            if (modifiedPeptide.charAt(adjustedModIndex + 1) == '[') {
                // previous mod present: remove it
                modifiedPeptide = removeModFromModifiedPeptide(modifiedPeptide, modLocation);
            }
        }
        modifiedPeptide = addModToModifiedPeptide(modifiedPeptide, modLocation, modMass);
        spLine.set(modPeptideCol, modifiedPeptide);
    }

    // add a modification to a modified peptide string
    public static String addModToModifiedPeptide(String previousModPep, int modLocation, double modMass) {
        if (modLocation == 0 || modLocation == previousModPep.length() + 1) {
            return previousModPep;   // do not edit terminal mods because glycan assignment will never adjust them
        }
        int residueIndex = getModPepIndex(previousModPep, modLocation);
        int AAindex = previousModPep.charAt(residueIndex) - 65;		// capital alphabet starts at 65
        int roundedModMass = (int) Math.round(modMass + AAMasses.monoisotopic_masses[AAindex]);     // round to match existing pepxml format
        return previousModPep.substring(0, residueIndex + 1) + String.format("[%d]", roundedModMass) + previousModPep.substring(residueIndex + 1);
    }

    // delete a modification from a modified peptide string
    public String removeModFromModifiedPeptide(String previousModPep, int modLocation) {
        if (modLocation == 0) {
            return previousModPep;   // do not edit terminal mods because glycan assignment will never adjust them
            // N-term mod. Remove n[##] from start (split on first ']')
//            int endBracket = previousModPep.indexOf(']') + 1;
//            return previousModPep.substring(endBracket);
        } else if (modLocation == peptide.length() + 1) {
            return previousModPep;   // do not edit terminal mods because glycan assignment will never adjust them
            // C-term mod. Remove c[##] from end (split on last '[')
//            int startBracket = previousModPep.lastIndexOf('[');
//            return previousModPep.substring(0, startBracket - 1);   // need to also remove preceding 'c'
        } else {
            // non-terminal mods
            int adjustedIndex = getModPepIndex(previousModPep, modLocation);

            // Find the start and end of the bracketed number
            int startBracket = previousModPep.indexOf('[', adjustedIndex);
            int endBracket = previousModPep.indexOf(']', startBracket) + 1;

            // Remove the bracketed number
            return previousModPep.substring(0, startBracket) + previousModPep.substring(endBracket);
        }
    }

    // get the 0-indexed site in a modified peptide String for the mod location (including other mods, if present)
    // Note: does NOT handle terminal mods (since glycan assignment will never adjust them)
    private static int getModPepIndex(String previousModPep, int modLocation) {
        int residueIndex = -1;
        int letterCount = 0;
        for (int i = 0; i <= previousModPep.length(); i++) {
            if (letterCount == modLocation) {
                residueIndex = i;
                break;
            }
            if (Character.isUpperCase(previousModPep.charAt(i))) {
                letterCount++;
            }
        }
        return residueIndex - 1;
    }

    public String toString() {
        return String.join("\t", this.spLine);
    }

}