import edu.umich.andykong.ptmshepherd.PTMShepherd;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class PTMShepherdTest {

    @Test
    void reNormName() {
        String specName = "filename.0001.0001.3";
        String expectedSpecName = "filename.1.1";
        String reNormedName = PTMShepherd.reNormName(specName);
        System.out.println(reNormedName);
        assertEquals(expectedSpecName, reNormedName);
    }
}
