# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
./gradlew build                   # Full build with tests
./gradlew test                    # Run all tests
./gradlew test --tests ClassName  # Run a single test class
./gradlew packageWithDeps         # Build fat JAR with dependencies
./gradlew shadowJar               # Alternative fat JAR target
```

Tests are in `test/` with resources in `test-resources/`. The CI workflow (`.github/workflows/gradle.yml`) runs `gradle build` and `gradle test` separately on Java 17.

Running the tool:
```bash
java -jar ptm-shepherd.jar path/to/config.txt
java -jar ptm-shepherd.jar --config    # Print default config
```

## Architecture Overview

PTM-Shepherd is a proteomics pipeline for characterizing post-translational modifications (PTMs) from mass spectrometry data. It consumes Philosopher `psm.tsv` output files and raw spectral data (mzML/mzXML).

### Entry Point & Global State

`PTMShepherd.java` is the main class and acts as a static singleton orchestrating the entire pipeline. It holds all global state: `params` (config map), `datasets` (TreeMap of experiment → PSM file paths), `mzMap` (spectral file lookup), and a shared `ExecutorService` for parallelism. Utility methods `getParam()`, `die()`, `print()`, and `normFName()` are used throughout the codebase.

### Pipeline Stages

The `main()` method runs these stages sequentially:
1. Load PSM files (`PSMFile`) and map spectral data files
2. Cache spectra to MZBin format for fast random access
3. Generate MS1 mass-shift histograms → detect peaks → annotate peaks with mod IDs
4. (Optional) Iterative localization refinement
5. Site localization scoring
6. Spectral similarity & retention time shift analysis
7. (Optional) Diagnostic feature mining
8. Combine results into summary tables

### Key Packages

| Package | Purpose |
|---------|---------|
| `peakpicker` | MS1 histogram generation, peak detection, Unimod annotation |
| `localization` | PTM site localization scoring against FASTA sequences |
| `iterativelocalization` | PTMiner-style iterative refinement with Poisson-binomial statistics and LDA |
| `specsimilarity` | Spectral similarity scoring and RT shift analysis |
| `diagnosticmining` | MS/MS diagnostic ion detection (Mann-Whitney U test) |
| `cleaner` | Aggregation of per-experiment results into combined output tables |
| `utils` | Stateless helpers: `Peptide`, `StringParsingUtils`, `Variance` |

### Core Data Models

- **`PSM.java`** — one peptide-spectrum match; holds delta mass, peptide sequence, modifications
- **`PSMFile.java`** — parses Philosopher `psm.tsv`; dynamically maps column headers; groups spectra by raw file
- **`Spectrum.java`** — MS/MS spectrum with peaks, precursor info, RT; multiple constructors for BatMass, MZBin, and MGF inputs
- **`Mod.java`** — lightweight (position, mass) tuple for a modification site

### Parallelism

All CPU-intensive analysis tasks submit work to a shared `ExecutorService` (thread count from `numThreads` param). Each `PSMFile` is processed independently, making the per-experiment loop the natural parallelism boundary.

