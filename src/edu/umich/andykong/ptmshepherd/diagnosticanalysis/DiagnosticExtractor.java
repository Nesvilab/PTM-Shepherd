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

package edu.umich.andykong.ptmshepherd.diagnosticanalysis;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PSM;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import edu.umich.andykong.ptmshepherd.core.MXMLReader;
import edu.umich.andykong.ptmshepherd.core.Spectrum;
import edu.umich.andykong.ptmshepherd.glyco.GlycoProfile;
import edu.umich.andykong.ptmshepherd.localization.SiteLocalization;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;


public class DiagnosticExtractor {
    String dsName;
    File rawDiagnosticFile;             // .diagnosticIon.tsv file
    MXMLReader mr;
    ArrayList<String> lineWithoutSpectra = new ArrayList<>();
    int totalLines;
    HashMap<String, MXMLReader> multiMr;
    float ppmTol, peakTol;
    int condPeaks;
    int specCol, pepCol, modpepCol, chargecol, deltaCol, rtCol, intCol, pmassCol, modCol;
    double condRatio;
    double[] capYShifts;
    double[] oxoniumIons;
    double[] remainderMasses;

    public DiagnosticExtractor(String dsName) {
        this.dsName = dsName;
        this.rawDiagnosticFile = new File(PTMShepherd.normFName(dsName + PTMShepherd.diagIonsExtractName));
    }


    public void extractDiagPSMs(PSMFile pf, HashMap<String, File> mzMappings, ExecutorService executorService, int numThreads) throws Exception {
        //open up output file
        HashMap<String, ArrayList<Integer>> mappings = new HashMap<>();
        PrintWriter diagnosticOut = new PrintWriter(new FileWriter(rawDiagnosticFile));
        ArrayList<String> linesWithoutSpectra = new ArrayList<>();

        //get necessary params
        ppmTol = Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"));
        condPeaks = Integer.parseInt(PTMShepherd.getParam("spectra_condPeaks"));
        condRatio = Double.parseDouble(PTMShepherd.getParam("spectra_condRatio"));
        //cap y ions
        String[] capYstrs;
        if (!PTMShepherd.getParam("cap_y_ions").isEmpty())
            capYstrs = PTMShepherd.getParam("cap_y_ions").split("[,/\\s]+");
        else
            capYstrs = new String[0];
        capYShifts = new double[capYstrs.length];
        for (int i = 0; i < capYstrs.length; i++)
            capYShifts[i] = Double.parseDouble(capYstrs[i]);
        //oxonium ions
        String[] oxStrs;
        if (!PTMShepherd.getParam("diag_ions").isEmpty())
            oxStrs = PTMShepherd.getParam("diag_ions").split("[,/\\s]+");
        else
            oxStrs = new String[0];
        oxoniumIons = new double[oxStrs.length];
        for (int i = 0; i < oxStrs.length; i++)
            oxoniumIons[i] = Double.parseDouble(oxStrs[i]);
        //remainder masses
        String[] remainderStrs;
        if (!PTMShepherd.getParam("remainder_masses").isEmpty())
            remainderStrs = PTMShepherd.getParam("remainder_masses").split("[,/\\s]+");
        else
            remainderStrs = new String[0];
        remainderMasses = new double[remainderStrs.length];
        for (int i = 0; i < remainderStrs.length; i++)
            remainderMasses[i] = Double.parseDouble(remainderStrs[i]);

        //write header
        StringBuilder diagnosticHeader = new StringBuilder(String.format("%s\t%s\t%s\t%s\t%s", "Spectrum", "Peptide", "Mods", "Pep Mass", "Mass Shift"));
        String ionTypes = PTMShepherd.concatIonTypes();

        for (double capYShift : capYShifts) diagnosticHeader.append(String.format("\tY_%.4f_intensity", capYShift));
        for (double oxoniumIon : oxoniumIons) diagnosticHeader.append(String.format("\tox_%.4f_intensity", oxoniumIon));
        for (double remainderMass : remainderMasses) {
            diagnosticHeader.append(String.format("\tdeltascore_%.4f\tlocalization_%.4f", remainderMass, remainderMass));
            for (int k=0; k < ionTypes.length(); k++) {
                diagnosticHeader.append(String.format("\tCount_%s_%.4f", ionTypes.charAt(k), remainderMass));
                diagnosticHeader.append(String.format("\tInt_%s_%.4f", ionTypes.charAt(k), remainderMass));
            }
        }
        diagnosticOut.println(diagnosticHeader);

        //map PSMs to file
        SiteLocalization.initSpectrumMappings(pf, mappings);

        /* Loop through spectral files -> indexed lines in PSM -> process each line */
        for (String cf : mappings.keySet()) { //for file in relevant spectral files
            long t1 = System.currentTimeMillis();
            //System.out.println(cf);
            mr = new MXMLReader(mzMappings.get(cf), Integer.parseInt(PTMShepherd.getParam("threads")));
            mr.readFully();
            long t2 = System.currentTimeMillis();
            ArrayList<Integer> clines = mappings.get(cf); //lines corr to curr spec file

            /* set up parallelization blocks */
            final int BLOCKSIZE = 100; //number of scans to be parsed per thread (to cut down on thread creation overhead)
            int nBlocks = clines.size() / (BLOCKSIZE); //number of jobs submitted to queue
            if (clines.size() % BLOCKSIZE != 0) //if there are missing scans, add one more block
                nBlocks++;

            ArrayList<Future> futureList = new ArrayList<>(nBlocks);
            /* Process PSM chunks */
            for (int i = 0; i < nBlocks; i++) {
                int startInd = i * BLOCKSIZE;
                int endInd = Math.min((i + 1) * BLOCKSIZE, clines.size());
                ArrayList<PSM> cBlock = new ArrayList<>();
                for (int j = startInd; j < endInd; j++)
                    cBlock.add(pf.psms.get(clines.get(j)));
                futureList.add(executorService.submit(() -> processLinesBlock(cBlock, diagnosticOut)));
            }
            /* Wait for all processes to finish */
            for (Future future : futureList)
                future.get();

            long t3 = System.currentTimeMillis();
            PTMShepherd.print(String.format("\t%s - %d (%d ms, %d ms)", cf, clines.size(), t2 - t1, t3 - t2));
        }
        diagnosticOut.close();

        if (!linesWithoutSpectra.isEmpty()) {
            PTMShepherd.print(String.format("\tCould not find %d/%d (%.1f%%) spectra.\n", linesWithoutSpectra.size(), this.totalLines,
                    100.0*((double)linesWithoutSpectra.size()/this.totalLines)));
            int previewSize = Math.min(linesWithoutSpectra.size(), 5);
            PTMShepherd.print(String.format("\tShowing first %d of %d spectra IDs that could not be found: \n\t%s\n", previewSize, linesWithoutSpectra.size(),
                    String.join("\n\t\t", linesWithoutSpectra.subList(0, previewSize))));
        }
    }

    public void processLinesBlock(ArrayList<PSM> cBlock, PrintWriter out) {
        StringBuilder newBlock  = new StringBuilder();
        for (PSM psm : cBlock) {
            newBlock.append(processLine(psm)).append("\n");
        }
        printLines(out, newBlock.toString());
    }

    private synchronized void printLines(PrintWriter out, String linesBlock) {
        out.print(linesBlock);
    }

    public String processLine(PSM psm) {
        StringBuilder diagnosticResultString = new StringBuilder();

        diagnosticResultString.append(String.format("%s\t%s\t%s\t%.4f\t%.4f", psm.getSpec(), psm.getPep(), psm.printAssignedMods(), psm.getCalcPepmass(), psm.getDMass()));

        Spectrum spec = mr.getSpectrum(psm.getSpec());
        if (spec == null) {
            this.lineWithoutSpectra.add(psm.getSpec());
            return "ERROR";
        }
        spec.conditionOptNorm(condPeaks, condRatio, false);

        //System.out.println("got spec");
        double[] capYIonIntensities;
        double[] oxoniumIonIntensities;
        capYIonIntensities = findCapitalYIonMasses(spec, psm.getCalcPepmass());
        oxoniumIonIntensities = findOxoniumIonMasses(spec, psm.getCalcPepmass());

        for (double capYIonIntensity : capYIonIntensities)
            diagnosticResultString.append(String.format("\t%.2f", capYIonIntensity));
        for (double oxoniumIonIntensity : oxoniumIonIntensities)
            diagnosticResultString.append(String.format("\t%.2f", oxoniumIonIntensity));
        float[] deltaScores = new float[remainderMasses.length];
        String ionTypes = PTMShepherd.concatIonTypes();
        float[][] remainderIntensities = new float[remainderMasses.length][ionTypes.length()];
        int[][] remainderCounts = new int[remainderMasses.length][ionTypes.length()];
        boolean[][] isMaxScores = localizeRemainderFragments(spec, psm.getPep(), psm.getAssignedMods(), deltaScores, remainderIntensities, remainderCounts);

        for (int i = 0; i < remainderMasses.length; i++) {
            diagnosticResultString.append(String.format("\t%.1f", deltaScores[i]));
            StringBuilder locSb = new StringBuilder("\t");
            for (int j = 0; j < psm.getPep().length(); j++) {
                if (isMaxScores[i][j]) {
                    locSb.append(String.format("%d%c", j + 1, psm.getPep().charAt(j))); //position (1 indexed), character
                }
            }
            // add remainder intensities
            for (int k=0; k < ionTypes.length(); k++) {
                locSb.append(String.format("\t%d", remainderCounts[i][k]));
                locSb.append(String.format("\t%.1f", remainderIntensities[i][k]));
            }
            diagnosticResultString.append(locSb);
        }
        return diagnosticResultString.toString();
    }

    public double[] findCapitalYIonMasses(Spectrum spec, double pepMass) {
        //implement charge states //todo

        int normToBasePeak = Integer.parseInt(PTMShepherd.getParam("glyco_cap_y_ions_normalize"));
        //System.out.println(normToBasePeak);

        //initialize final capYion masses
        double[] capYIons = new double[capYShifts.length];
        double [] capYIonIntensities = new double[capYShifts.length];
        for (int i = 0; i < capYIons.length; i++)
            capYIons[i] = capYShifts[i] + pepMass;
        //find capital Y ion intensities
        for (int i = 0; i < capYIons.length; i++) {
            //System.out.println(capYIons[i]);
            capYIonIntensities[i] = spec.findIonNeutral(capYIons[i],
                    Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol")),
                    Integer.parseInt(PTMShepherd.getParam("spectra_maxPrecursorCharge"))); //todo simplify parameter calling
            if (normToBasePeak == 1) {
                //System.out.print(capYIonIntensities[i]);
                //System.out.println(" 1");
                //System.out.print(spec.findBasePeakInt());
                //System.out.println(" 2");
                capYIonIntensities[i] /= spec.basePeakInt;
                capYIonIntensities[i] *= 100.0;
                //System.out.print(capYIonIntensities[i]);
                //System.out.println(" 3");
            }
        }
        return capYIonIntensities;
    }

    public double[] findOxoniumIonMasses(Spectrum spec, double pepMass) {
        //initialize oxonium masses //todo
        //initialize capYion masses /todo
        //implement charge states //todo
        //initialize oxonium ion intensities
        int normToBasePeak = Integer.parseInt(PTMShepherd.getParam("glyco_diag_ions_normalize"));
        double[] oxoniumIonIntensities = new double[oxoniumIons.length];
        //for ion in oxonium masses/capYions
        for (int i = 0; i < oxoniumIons.length; i++) {
            oxoniumIonIntensities[i] = spec.findIon(oxoniumIons[i], Float.parseFloat(PTMShepherd.getParam("spectra_ppmtol"))); //todo simplify parameter calling
            if (normToBasePeak == 1) {
                //System.out.print(oxoniumIonIntensities[i]);
                //System.out.println(" 1");
                //System.out.print(spec.findBasePeakInt());
                //System.out.println(" 2");
                oxoniumIonIntensities[i] /= spec.basePeakInt;
                oxoniumIonIntensities[i] *= 100.0;
                //System.out.print(oxoniumIonIntensities[i]);
                //System.out.println(" 3");
            }
        }
        return oxoniumIonIntensities;
    }

    public boolean[][] localizeRemainderFragments(Spectrum spec, String seq, TreeMap<Integer, Float> smods, float[] deltaScores, float[][] remainderInts, int[][] remainderCounts) {
        //initialize allowed positions
        boolean [] allowedPoses = SiteLocalization.parseAllowedPositions(seq, PTMShepherd.getParam("localization_allowed_res"));
        //initialize remainder delta scores
        //double[] remainderDscores = new double[remainderMasses.length];
        //add variable and fixed mods to frag masses for peptide
        float [] mods = new float[seq.length()];
        Arrays.fill(mods, 0f);
        SiteLocalization.localizeMods(smods, mods);
        //iterate through remainder masses
        //these 3 variables store values for each remainder mass
        float [] maxScores = new float[remainderMasses.length];
        int [] maxFrags = new int[remainderMasses.length];
        boolean [][] isMaxScores = new boolean[remainderMasses.length][seq.length()]; //1 if localized AND = max score
        //these 3 variables store values that are constant for the PSM

        String ionTypes = PTMShepherd.concatIonTypes();
        float baseScore = spec.getHyper(seq, mods, ppmTol);
        int baseFrags = spec.getFrags(seq, mods, ppmTol);
        //these 3 variables need to be reinitialized every remainder mass
        float [] scores;
        int [] frags;
        double dmass;
        //begin iterate through remainder masses
        for (int i = 0; i < remainderMasses.length; i++){
            //reinit for each remainder mass
            dmass = remainderMasses[i];
            scores = new float[seq.length()];
            frags = new int[seq.length()];
            maxScores[i] = baseScore;
            maxFrags[i] = baseFrags;

            //localize at each position
            for(int j = 0; j < seq.length(); j++) {
                if (allowedPoses[j])
                    mods[j] += dmass;
                scores[j] = spec.getHyper(seq, mods, ppmTol);
                //System.out.println(scores[j] + "score");
                frags[j] = spec.getFrags(seq, mods, ppmTol);
                if(frags[j] > maxFrags[i]) {
                    maxFrags[i] = frags[j];
                }
                if(scores[j] > maxScores[i]) {
                    maxScores[i] = scores[j];
                }
                if (allowedPoses[j])
                    mods[j] -= dmass;
            }

            //System.out.println(maxScores[i]+"maxscore");
            //determine if localized and record max positions
            if (maxScores[i] > baseScore) {
                deltaScores[i] = maxScores[i] - baseScore;
                for (int j = 0; j < seq.length(); j++) {
                    if (scores[j] == maxScores[i]) {
                        isMaxScores[i][j] = true;
                        // save remainder fragment intensities for best position only
                        TreeMap<Character, Float> ionIntensities = initIonIntensities();
                        TreeMap<Character, Integer> ionCounts = initIonCounts();
                        mods[j] += (float) dmass;
                        spec.getRemainderFrags(seq, mods, ppmTol, ionIntensities, ionCounts, j);
                        mods[j] -= (float) dmass;
                        for (int k=0; k < ionIntensities.size(); k++) {
                            remainderInts[i][k] = ionIntensities.get(ionTypes.charAt(k));
                            remainderCounts[i][k] = ionCounts.get(ionTypes.charAt(k));
                        }
                    } else {
                        isMaxScores[i][j] = false;
                    }
                }
            } else {
                for (int j = 0; j < seq.length(); j++) {
                    isMaxScores[i][j] = false;
                }
            }
        }
        return isMaxScores;
    }


    public boolean isDiagnosticComplete() throws Exception {
        if(rawDiagnosticFile.exists()) {
            RandomAccessFile raf = new RandomAccessFile(rawDiagnosticFile, "r");
            raf.seek(Math.max(0, rawDiagnosticFile.length() - 20));
            String cline;
            while((cline = raf.readLine())!=null)
                if(cline.equals("COMPLETE")) {
                    raf.close();
                    return true;
                }
            raf.close();
            rawDiagnosticFile.delete();
        }
        return false;
    }

    private TreeMap<Character, Float> initIonIntensities() {
        TreeMap<Character, Float> ionIntensities = new TreeMap<>();
        String ionTypes = PTMShepherd.concatIonTypes();
        for (int i=0; i < ionTypes.length(); i++) {
            ionIntensities.put(ionTypes.charAt(i), 0.0F);
        }
        return ionIntensities;
    }

    private TreeMap<Character, Integer> initIonCounts() {
        TreeMap<Character, Integer> ionCounts = new TreeMap<>();
        String ionTypes = PTMShepherd.concatIonTypes();
        for (int i=0; i < ionTypes.length(); i++) {
            ionCounts.put(ionTypes.charAt(i), 0);
        }
        return ionCounts;
    }

    public void completeDiagnostic() throws Exception {
        PrintWriter out = new PrintWriter(new FileWriter(rawDiagnosticFile,true));
        out.println("COMPLETE");
        out.close();
    }
    public void updateGlycoProfiles(GlycoProfile[] profiles) throws Exception {
        BufferedReader in = new BufferedReader(new FileReader(rawDiagnosticFile));
        String cline;
        in.readLine();
        while ((cline = in.readLine()) != null) {
            if (cline.equals("COMPLETE"))
                break;
            if (cline.startsWith("Spectrum"))
                continue;
            if (cline.startsWith("ERROR"))
                continue;
            String[] sp = cline.split("\\t");
            double md = Double.parseDouble(sp[4]);
            for (int i = 0; i < profiles.length; i++) {
                int cind = profiles[i].locate.getIndex(md);
                if (cind != -1) {
                    profiles[i].records[cind].updateWithLine(sp);
                }
            }
        }
        in.close();
    }
}
