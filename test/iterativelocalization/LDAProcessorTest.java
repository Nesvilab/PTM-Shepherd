package iterativelocalization;

import edu.umich.andykong.ptmshepherd.iterativelocalization.LDAProcessor;
import edu.umich.andykong.ptmshepherd.iterativelocalization.MatchedIonTable;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import static org.junit.jupiter.api.Assertions.*;


public class LDAProcessorTest {
    private MatchedIonTable data;

    @BeforeEach
    public void setUp() {
        // Initialize test data
        data = new MatchedIonTable();
        // Adding sample data points
        data.rows.add(data.new MatchedIonRow(2.0f, 1.1f, false));
        data.rows.add(data.new MatchedIonRow(2.5f, 1.2f, false));
        data.rows.add(data.new MatchedIonRow(2.4f, 1.3f, false));
        data.rows.add(data.new MatchedIonRow(2.2f, 2.2f, true));
        data.rows.add(data.new MatchedIonRow(2.3f, 2.4f, true));
        data.rows.add(data.new MatchedIonRow(2.2f, 2.3f, true));
    }


    @Test
    public void testSolveLDA() {
        LDAProcessor ldaProcessor = new LDAProcessor(data);

        assertDoesNotThrow(() -> ldaProcessor.solveLDA(Executors.newFixedThreadPool(2)));

        // Check that values are equivalent when project onto learned EigenMatrix
        assertEquals(ldaProcessor.projectedData.getEntry(0,0),
                ldaProcessor.projectData(Math.log10(0.1+2.0f), Math.log10(1.1f)).getEntry(0,0));

    }
}

