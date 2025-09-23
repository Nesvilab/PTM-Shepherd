package glyco;

import edu.umich.andykong.ptmshepherd.PSMFile;
import edu.umich.andykong.ptmshepherd.PTMShepherd;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import utils.TestUtilities;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static edu.umich.andykong.ptmshepherd.PTMShepherd.makeOutputDir;

public class GlycoPsmFilterTest {
    protected File outputDir;

    @TempDir
    Path tempOutputDir;

    @Test
    public void GlycoPSMFilterTests() {
        String singleConfigPath = new File("test-resources/shepherd.config").getAbsolutePath();
        String multiConfigPath = new File("test-resources/2expts/shepherd.config").getAbsolutePath();

        runGlycoAnalysisTest(singleConfigPath);
        runGlycoAnalysisTest(multiConfigPath);
    }

    public void runGlycoAnalysisTest(String configPath) {
        HashMap<String, Integer> firstGlycoPSMCounts = new HashMap<>();
        HashMap<String, Integer> initalPSMCounts = new HashMap<>();

        // Set up temp output directory
        outputDir = tempOutputDir.toFile();
        makeOutputDir(outputDir.getAbsolutePath());

        // initialize PTM-Shepherd
        PTMShepherd.init(new String[] {configPath});
        PTMShepherd.params.put("output_path", outputDir.getAbsolutePath());
        PTMShepherd.outputPath = outputDir.getAbsolutePath();

        // copy the psm files to temp dir (to avoid overwriting the originals)
        try {
            copyPsmsToTempDir();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        // do PTM-S setup
        PTMShepherd.loadPSMFiles();
        PTMShepherd.print("Finding spectral data");
        PTMShepherd.getMzDataMapping();
        PTMShepherd.print("Done finding spectral data\n");
        PTMShepherd.deletePreviousFiles();
        PTMShepherd.print("Caching spectral data");
        PTMShepherd.rewriteDataToMzBin();
        PTMShepherd.print("Done caching spectral data\n");

        // run glycan assignment test twice, counting PSMs before and after each run to confirm that numbers remain consistent
        for (Map.Entry<String, ArrayList<PSMFile>> entry : PTMShepherd.psmFiles.entrySet()) {
            for (PSMFile psmFile : entry.getValue()) {
                TestUtilities.countGlycoPSMs(psmFile, 1);
                initalPSMCounts.put(entry.getKey(), psmFile.psms.size());
            }
        }
        PTMShepherd.runGlycanAssignment();
        for (Map.Entry<String, ArrayList<PSMFile>> entry : PTMShepherd.psmFiles.entrySet()) {
            for (PSMFile psmFile : entry.getValue()) {
                int gPsms = TestUtilities.countGlycoPSMs(psmFile, PTMShepherd.glycoParams.glycoFDR);
                firstGlycoPSMCounts.put(entry.getKey(), gPsms);
                assert psmFile.psms.size() == initalPSMCounts.get(entry.getKey());
            }
        }

        // 2nd run
        PTMShepherd.deletePreviousFiles();
        PTMShepherd.loadPSMFiles();     // re-load PSMs to make sure loading the unfiltered tables works correctly
        PTMShepherd.runGlycanAssignment();

        for (Map.Entry<String, ArrayList<PSMFile>> entry : PTMShepherd.psmFiles.entrySet()) {
            for (PSMFile psmFile : entry.getValue()) {
                int SecondGPSMs = TestUtilities.countGlycoPSMs(psmFile, PTMShepherd.glycoParams.glycoFDR);
                assert SecondGPSMs == firstGlycoPSMCounts.get(entry.getKey());
                assert psmFile.psms.size() == initalPSMCounts.get(entry.getKey());
            }
        }
    }

    public void copyPsmsToTempDir() throws IOException {
        // copy PSM files to temp dir
        for (Map.Entry<String, ArrayList<String[]>> datasetEntry: PTMShepherd.datasets.entrySet()) {
            String expName = datasetEntry.getKey();
            for (String[] input : datasetEntry.getValue()) {
                String psmPath = input[0];
                Path psmPathObj = Paths.get(psmPath);
                Path newExpDir = Paths.get(outputDir.getAbsolutePath(), expName);
                if (!Files.exists(newExpDir)) {
                    Files.createDirectories(newExpDir);
                }
                Path destPath = Paths.get(newExpDir.toString(), psmPathObj.getFileName().toString());
                Files.copy(psmPathObj, destPath);
                // update the path in PTMShepherd
                input[0] = destPath.toString();
            }
        }

    }
}
