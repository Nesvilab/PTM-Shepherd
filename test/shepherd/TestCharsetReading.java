package shepherd;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

public class TestCharsetReading {

//  @Test
  public void readFileWithBadSymbolTest() throws Exception {
    String fn = "/peaksummary.annotated.tsv";
    System.out.println("Opening stream");
    try (InputStream is = this.getClass().getResourceAsStream(fn)) {
      if (is == null) {
        System.out.println("Stream didn't open, exiting");
        return;
      }
      InputStreamReader isr = new InputStreamReader(is);
      BufferedReader br = new BufferedReader(isr);
      String line;
      System.out.println("Start reading");
      while ((line = br.readLine()) != null) {
        System.out.println("read line: " + line);
      }
      System.out.println("Done reading");
    }
  }

//  @Test
  public void readAllLines() throws Exception {
    String path = "C:\\code\\msfragger\\ptmshepherd\\test-resources\\peaksummary.annotated-from-shepherd.tsv";
    List<String> strings = Files.readAllLines(Paths.get(path), StandardCharsets.UTF_8);

    for (String string : strings) {
      System.out.println("Line: " + string);
    }
  }


//  @Test
  public void readFile() throws Exception {
    String path = "C:\\code\\msfragger\\ptmshepherd\\test-resources\\peaksummary.annotated-from-shepherd.tsv";

    try (InputStream is = Files.newInputStream(Paths.get(path))) {
      if (is == null) {
        System.out.println("Stream didn't open, exiting");
        return;
      }
      InputStreamReader isr = new InputStreamReader(is, "UTF-8");
      BufferedReader br = new BufferedReader(isr);
      String line;
      System.out.println("Start reading");
      while ((line = br.readLine()) != null) {
        System.out.println("read line: " + line);
      }
      System.out.println("Done reading");
    }
  }

//  @Test
  public void newBufferedReader() throws Exception {
    String path = "C:\\code\\msfragger\\ptmshepherd\\test-resources\\peaksummary.annotated-from-shepherd.tsv";

    try (BufferedReader br = Files.newBufferedReader(Paths.get(path), StandardCharsets.UTF_8)) {
      if (br == null) {
        System.out.println("Stream didn't open, exiting");
        return;
      }
      String line;
      System.out.println("Start reading");
      while ((line = br.readLine()) != null) {
        System.out.println("read line: " + line);
      }
      System.out.println("Done reading");
    }
  }

}
