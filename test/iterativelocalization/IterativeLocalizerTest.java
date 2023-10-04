package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.IterativeLocalizer;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class IterativeLocalizerTest {

    @Test
    // TODO write tests here
     void computePoissonBinomialLikelihood() {
        String pep = "PASGAGAGAGAGKR";
        float[] mods = new float[]{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        float dMass = 100.0188f;
        boolean[] allowedPoses = new boolean[]{true, true, true, true, true, true, true, true, true, true, true, true,
                true, true, true, true};
        float[] peakMzs = new float[]{169.09694f, 175.11827f, 256.12845f, 313.1513f, 384.1857f, 512.244f, 588.312f,
                716.3641f, 844.4323f, 972.4885f, 1059.5171f, 1130.5404f};
        float[] peakInts = new float[]{19.48378f, 3.9259293f, 6.327259f, 5.4565396f, 6.5460205f, 4.2508264f, 7.016416f,
                23.633015f, 29.399632f, 7.5406113f, 30.278028f, 4.4198275f};

    }


    @Test // Can use reflection to test private methods
    void findMatchedIons() { // todo test
        // Add test case where there are some overlapping ions GGNFGGRGGYGGGGGGSR / 02330a_GB1_3990_02_PTM_TrainKit_Rmod_Methyl_200fmol_3xHCD_R1.8156.8156"
        // Another one 02330a_GB1_3990_02_PTM_TrainKit_Rmod_Methyl_200fmol_3xHCD_R1.34399.34399
    }
}
