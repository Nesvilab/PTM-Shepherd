PTM-Shepherd automates characterization of PTM profiles detected in open searches based on attributes such as amino acid localization, fragmentation spectra similarity, retention time shifts, and relative modification rates. PTM-Shepherd can also perform multi-experiment comparisons for studying changes in modification profiles, e.g. in data generated in different laboratories or under different conditions.

![PTM-Shepherd Workflow](./Workflow.png)

Data processing begins by aggregating the mass shifts across all datasets into a common histogram. Peaks are determined based on their prominence. The 500 most intense peaks in aggregate are then quantified for each dataset and normalized to size. Peptides with each mass shift are iteratively rescored with the peptide at each position, producing localization scores for each peptide and an aggregate localization enrichment for each mass shift. Finally, modified peptides and their unmodified counterparts are analyzed to have their pairwise cosine spectral similarity and change in retention time calculated.

### [Documentation](https://github.com/Nesvilab/PTM-Shepherd/wiki)

### Usage

We strongly recommend running PTM-Shepherd via [FragPipe](http://fragpipe.nesvilab.org/) to simplify open search analysis. PTM-Shepherd can also be run from the command line as a standalone JAR file. 
* You can download the latest Fragpipe release [here](https://github.com/Nesvilab/FragPipe/releases).
* You can download the latest standlone release of PTM-Shepherd [here](https://github.com/Nesvilab/PTM-Shepherd/releases).

#### Command line usage
If running via command line, parameters should be passed to PTM-Shepherd as a text configuration file. PTM-Shepherd can be executed by the command line by running
```
java -jar ptm-shepherd.jar path/to/config.txt
```
Datasets are passed to PTM-Shepherd in the form of [Philosopher](https://philosopher.nesvilab.org/) psm.tsv files. If experiments (parameter: dataset) are given different names, they will be analyzed as separate experiments. If experiments share the same name, they will be pooled together during analysis. To pass datasets to PTM-Shepherd, datasets should be presented in the configuration file (one line per psm.tsv file) as:
```
dataset = $DATASETNAME01 path/to/psm01.tsv path/to/mzML/directory
dataset = $DATASETNAME02 path/to/psm02.tsv path/to/mzML/directory
```

### Optional parameters
```

threads = 8                 #number of threads used for processing. Default is either 8 or the number of available threads, whichever is lower.

#MS1 delta mass histogram parameters
histo_bindivs = 5000        #takes integer values > 0. Number of bins per dalton to be used for mass shift binning. The default is 5000 bins, or 0.0002 Da bins.
histo_smoothbins = 3        #takes integer values >= 0. Number of bins on each side of a bin that the weight of the bin is smoothed across. This smoothing traces a normal distribution. A value of 1, e.g., will smooth the weight of the bin across 1 bin to either side (3 bins total) using the weights 0.23 (bin to left), 0.49 (same bin), 0.23 (bin to right). Changing this parameter is not recommended for non-advanced users.
histo_normalizeTo = psms    #takes "psms" and "scans". Normalize output to total PSMs in dataset. 

#MS1 delta mass peakpicking parameters
peakpicking_promRatio = 0.3 #takes values between 0 and 1. Ratio of peak apex to peak height used to call peaks. Peak height is determined by checking the shorter of the two peak shoulders. Changing this parameter is not recommended for non-advanced users.
peakpicking_mass_units = 0  #units in which peakpicking is done, 0  = Da and 1 = ppm. PPM peakpicking widths are calculated based on an average peptide mass plus the mass shift. Default is 0. Changing this parameter is not recommended for non-advanced users.
peakpicking_width = 0.002   #takes values > 0. Width of signal used in signal-to-noise calculation for peakpicking. This parameter is applied symmetrically, so the true width is 2\*peakpicking_width. Changing this parameter is not recommended for non-advanced users.
peakpicking_topN = 500      #takes values > 0. Number of peaks reported in output. Larger numbers have very small effects on processing time, so it is rarely beneficial to reduce the number of peaks reported.
peakpicking_minPsm = 10     #takes values > 0. Number of PSMs required for a peak to be called.
mass_offsets = 0            #takes a / separated list of values, at which modification peaks will be checked for (e.g. 0/79.9663). Default is unused.
isotope_error = 0           #takes a / separated list of isotope states that modify mass offsets to check for combinations of a mass shift and isotopic peaks (e.g. 0/1/2). Default is unused.

#PSM parameters
precursor_mass_units = 0    #units in which peak widths are selected, 0  = Da and 1 = ppm. PPM peak widths are calculated based on an average peptide mass plus the mass shift. Default is 0. Changing this parameter is not recommended for non-advanced users.
precursor_tol = 0.02        #takes values > 0. Width of annotated peaks. Will assign final peak widths to either precursor_tol or detected width, whichever is more precise. Default is 0.01. Also determines which spectra are “zero bin” when doing modified to unmodified comparisons.
spectra_ppmtol = 20         #takes number > 0. This is the tolerance applied when matching peaks for localization and similarity scoring
spectra_condPeaks = 150     #takes integer > 0. This is the number of peaks that will be used when performing localization and similarity scoring.
spectra_condRatio = 0.01    #takes number between 0 and 1. Minimum intensity used for peaks in localization and similarity scoring, expressed as a ratio of a peak’s intensity to the most intense peak.

#Annotation parameters
varmod_masses = None:0.0    #takes a series of mass shifts and names expressed like Oxidation:15.9949,Acetylation:42.0106. These mass shifts are used and prioritized during peak annotation. The default is a failed carbamidomethylation event, so that modifications occurring on cysteine can be expressed as a combination of those two modifications. This is useful when working with enriched datasets or datasets with PTMs not enumerated in Unimod.
annotation_file =           #takes path to custom annotation file. Must be of tab-separated lines of format $MODNAME\t$MODMASS.
annotation_tol = 0.01       #annotation tolerance (in daltons) for unimod matching to MS1 peak apex from histogram

#Localization
localization_background = 3 #takes integer values 1-4. Background of peptides against which localization enrichment scores are calculated.
                              #1 = unique peptides with localizable mass shifts within the same bin are used as background
                              #2 = PSMs with localizable mass shifts within the same bin are used as background
                              #3 = unique peptides with localizable mass shifts within the entire dataset are used as background
                              #4 = PSMs with localizable mass shifts within the entire dataset are used as background
localization_allowed_res = all  #takes all or ABCDEF. Restrcts localization to these amino acids.
#Use iontype_Q in localiztions. Takes 0 or 1.
iontype_a = 0              
iontype_b = 1
iontype_c = 0
iontype_x = 0
iontype_y = 1
iontype_z = 0

#Extract known diagnostic features NOTE:parameters will be migrated to different names in version 3.0
diagextract_mode = false     #enable diagnostic feature extractions
#extract peptide remainder masses from spectra (comma, space, slash separation accepted)
cap_y_ions = 0,203.07937,349.137279,406.15874,568.21156,730.26438,892.3172
#extract diagnostic ions from spectra (comma, space, slash separation accepted)
diag_ions = 144.0656,138.055,168.065526,186.076086,204.086646,243.026426,274.0921325,292.1026925,308.09761,366.139466,405.079246,485.045576,512.197375,657.2349
#extract fragment remainder masses from spectra and localize (comma, space, slash separation accepted)
remainder_masses = 203.07937

#Diagnostic feature mining
diagmine_mode = false      #turns on diagnostic feature mining mode

diagmine_minSignal = 0.001 #controls the minimum signal required per MS2 feature to be send to statistical filtering. Accepts values >0.
diagmine_filterIonTypes = aby #defines which fragment ion series are removed from spectra before feature calculation. Accepts abcxyz.
diagmine_ionTypes = by      #defines which fragment ion series will have features calculated. Accepts abcxyz.
diagmine_maxP = 0.05        #maximum E-value cutoff for spectral features modified and unmodified peptides. Set to large number (e.g. 10000) to not filter by E-value.
diagmine_minPeps = 25 #minimum number of unique peptide ions that must be present for a mass shift to be processed
diagmine_maxPsms = 1000 #maximum number of PSMs that will be processed per mass shift
diagmine_printIsotopes = false #retains collapsed isotopic peaks instead of filtering them out. Takes true/false.

#The following parameters control cutoffs for peptide remainder masses
diagmine_pepMinSpecDiff = 25.0 #minimum percentage of spectra containing feature
diagmine_pepMinFoldChange = 3.0 #minimum fold-change difference in average intensity between modified and unmodified peptides

#The following parameters control cutoffs for diagnostic ions masses
diagmine_diagMinSpecDiff = 25.0 #minimum percentage of spectra containing feature
diagmine_diagMinFoldChange = 3.0 #minimum fold-change difference in average intensity between modified and unmodified peptides

#The following parameters control cutoffs for fragment remainder masses
diagmine_minIonsPerSpec = 2 #minumum number of shifted ions that must be present in an ion series for it to be considered a hit
                            #decrease to increase sensitivty and noise
diagmine_fragMinSpecDiff = 15.0 #minimum percentage of spectra containing feature
diagmine_fragMinPropensity = 12.5 #minimum value for n_shifted_ions / n_shifted_ions * 100 that must be present for an ion series
diagmine_fragMinFoldChange = 3.0  #minimum fold-change difference in average intensity between

#The following parameters control output
output_extended = false #takes true/false. Prints additional files related to the analysis. Useful if you are interested in spectrum-level analysis rather than aggregate analysis. Default is false.
output_path =   #direct output to directory
```

### How to cite
#### For all PTM-Shepherd uses, please cite this manuscript:
Daniel J. Geiszler, Andy T. Kong, Dmitry M. Avtonomov, Fengchao Yu, Felipe V. Leprevost, Alexey I. Nesvizhskii. *PTM-Shepherd: analysis and summarization of post-translational and chemical modifications from open search results*. doi: (https://doi.org/10.1074/mcp.TIR120.002216)[https://doi.org/10.1074/mcp.TIR120.002216].

#### For glycan assignment, please cite this manuscript:
Daniel A. Polasky, Daniel J. Geiszler, Fengchao Yu, Alexey I. Nesvizhskii. *Multiattribute Glycan Identification and FDR Control for Glycoproteomics*. doi: (https://doi.org/10.1016/j.mcpro.2022.100205)[https://doi.org/10.1016/j.mcpro.2022.100205].

##### For diagnostic feature mining, please cite this manuscript:
 Daniel J. Geiszler, Daniel A. Polasky, Fengchao Yu, Alexey I. Nesvizhskii. *Mining for ions: diagnostic feature detection in MS/MS spectra of post-translationally modified peptides*. doi: (https://doi.org/10.1101/2022.09.12.507594)[https://doi.org/10.1101/2022.09.12.507594].

### Check out our other tools!

* [FragPipe](https://fragpipe.nesvilab.org/)
* [MSFragger](https://msfragger.nesvilab.org)
* [Philosopher](https://philosopher.nesvilab.org)
* [IonQuant](https://ionquant.nesvilab.org)

